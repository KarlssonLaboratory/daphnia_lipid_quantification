# =============================================================================
# Daphnia Lipid Droplet Analysis Pipeline
# =============================================================================
# Loads per-well imaging data and droplet-level segmentation data from
# multiple experimental replicates, fits GLMs / LMMs per response variable,
# and writes diagnostic plots plus pairwise-comparison tables to disk.
# =============================================================================

# ── Install and load dependencies ─────────────────────────────────────────────

packages <- c(
  "tidyverse",
  "readxl",
  "multcomp",
  "lme4",
  "emmeans",
  "DHARMa",
  "MASS",
  "patchwork"
)

install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

invisible(lapply(packages, install_and_load))


# =============================================================================
# Response-variable specification table
# =============================================================================
# Each row names one response variable and the GLM family to use for it.
# Add or remove rows here to change what the pipeline analyses.
# "negbin" is a sentinel string handled specially inside fit_glm() for negative binomial models.
# =============================================================================
response <- tibble::tibble(
  response = c(
    "droplet_intensity_total",  # total lipid-droplet fluorescence intensity
    "droplet_area_total",       # summed area of all droplets per animal
    "daphnia_size",             # body-size proxy (e.g. pixel area)
    "num_droplets",             # raw count of detected droplets
    "num_droplets_per_area"     # count normalised by Daphnia body area
  ),
  family = list(
    Gamma(link = "log"),  # strictly positive, right-skewed intensity data
    gaussian(),           # area totals – approximately normal after logging
    gaussian(),           # size – approximately normal
    "negbin",             # count data with overdispersion
    Gamma(link = "log")   # ratio – strictly positive, right-skewed
  )
)


# =============================================================================
# load_experiment()
# =============================================================================
# Reads all data files for one experiment folder and returns a single
# tidy data frame enriched with treatment labels and a normalised droplet
# density column.  Also writes the individual-droplet data frame into the
# global environment as `droplet_data` (used later by analyze_experiment).
#
# Arguments:
#   experiment  – character; name of the sub-folder under "data/"
#
# Returns:
#   A data frame with one row per well × frame observation.
# =============================================================================
load_experiment <- function(experiment) {

  data_dir <- file.path("data/TBT-CL chronic mothers", experiment)

  # ── 1. Load well-treatment labels ──────────────────────────────────────────
  # Each replicate has an Excel file "*_wells.xlsx" mapping well letter to
  # food_treatment and `ug/L`.  We stack all replicates into one frame
  # and retain the replicate ID as a factor.
  well_files <- list.files(
    data_dir,
    pattern    = ".*_wells.*\\.xlsx",
    full.names = TRUE
  )

  treatment <- map_dfr(well_files, function(f) {
    rep_id <- sub(".*measures_export_(.*)\\_wells.*\\.xlsx", "\\1", basename(f))
    df <- read_excel(f, col_names = TRUE)
    df$replicate <- rep_id
    df
  })
  treatment$replicate <- as.factor(treatment$replicate)

  # ── 2. Load per-well summary CSVs ──────────────────────────────────────────
  # One CSV per replicate; filename encodes the replicate ID.
  data_files <- list.files(
    data_dir,
    pattern    = ".*\\.csv",
    full.names = TRUE
  )

  data <- map_dfr(data_files, function(f) {
    rep_id <- sub("measures_export_(\\d+).*\\.csv", "\\1", basename(f))
    df <- read.csv(f)
    df$replicate  <- rep_id
    df$experiment <- experiment
    df
  })
  data$replicate <- as.factor(data$replicate)

  # ── 3. Attach treatment labels ─────────────────────────────────────────────
  # Join on (well letter, replicate) so each well row gets its treatment group.
  data <- data |>
    mutate(well_letter = substr(well_name, 1, 1)) |>
    inner_join(treatment, by = c("well_letter" = "well", "replicate"))

  # Enforce meaningful factor ordering for food treatment
  data$food_treatment <- factor(
    data$food_treatment,
    levels = c("no food", "low food", "high food")
  )
  data$`ug/L` <- factor(data$`ug/L`)

  # ── 4. Derived variable: droplets per unit body area ───────────────────────
  data <- data |>
    mutate(num_droplets_per_area = num_droplets / daphnia_area)

  # ── 5. Quality control: remove wells with ≥ 5 zero-valued numeric columns ──
  # Wells with this many zeros are likely imaging failures rather than
  # biologically meaningful observations.
  data <- data |>
    rowwise() |>
    mutate(zero_count = sum(c_across(where(is.numeric)) == 0, na.rm = TRUE)) |>
    ungroup()

  removed_rows <- data |> filter(zero_count >= 5)

  if (nrow(removed_rows) > 0) {
    warning(
      "Rows removed due to >=5 zeros:\n",
      paste0(
        "  Well ", removed_rows$well_name,
        " (", removed_rows$zero_count, " zeros)",
        collapse = "\n"
      )
    )
  }

  data <- data |>
    dplyr::filter(zero_count < 5) |>
    dplyr::select(-zero_count)

  # ── 6. Load individual-droplet segmentation data ───────────────────────────
  # The masks_<rep_id>/ sub-folders contain one CSV per well with one row per
  # detected droplet.  We stack them, attach treatment labels, and store in the
  # global environment so analyze_experiment() can access them.
  masks_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  masks_dirs <- masks_dirs[grepl("masks_", basename(masks_dirs))]

  droplet_data <- map_dfr(masks_dirs, function(masks_dir) {
    rep_id     <- sub(".*masks_(\\d+)", "\\1", basename(masks_dir))
    mask_files <- list.files(
      masks_dir,
      pattern    = "_droplets\\.csv$",
      full.names = TRUE
    )

    map_dfr(mask_files, function(f) {
      df <- tryCatch(read.csv(f), error = function(e) data.frame())
      if (nrow(df) == 0) return(data.frame())        # skip empty files
      df$well_name  <- sub("^([A-Z]\\d+)_.*", "\\1", basename(f))
      df$replicate  <- rep_id
      df$experiment <- experiment
      df
    })
  })

  if (nrow(droplet_data) > 0) {
    droplet_data <- droplet_data |>
      mutate(
        well_letter = substr(well_name, 1, 1),
        replicate   = as.factor(replicate)
      ) |>
      inner_join(treatment, by = c("well_letter" = "well", "replicate"))

    droplet_data$food_treatment <- factor(
      droplet_data$food_treatment,
      levels = c("no food", "low food", "high food")
    )
    droplet_data$`ug/L` <- factor(droplet_data$`ug/L`)

    # Expose to the calling environment so analyze_experiment() can use it
    assign("droplet_data", droplet_data, envir = .GlobalEnv)
  }

  return(data)
}


# =============================================================================
# Helpers used inside analyze_experiment()
# =============================================================================

# ── resolve_family() ──────────────────────────────────────────────────────────
# Converts the raw family specification (family object, function, or character
# string) from the `response` table into a usable object.  Returns the string
# "negbin" unchanged so fit_glm() can route to glm.nb().
resolve_family <- function(fam) {
  if (is.character(fam) && fam == "negbin") return("negbin")
  if (inherits(fam, "family"))              return(fam)
  if (is.function(fam))                     return(fam())
  # last resort: parse a character expression such as "Gamma(link='log')"
  tryCatch(
    match.fun(fam)(),
    error = function(e) eval(parse(text = fam))
  )
}

# ── fit_glm() ─────────────────────────────────────────────────────────────────
# Thin wrapper that dispatches to glm.nb() for negative-binomial responses
# and to glm() for everything else.
fit_glm <- function(formula, family_obj, data, is_negbin = FALSE) {
  if (is_negbin) {
    MASS::glm.nb(formula, data = data)
  } else {
    glm(formula, family = family_obj, data = data)
  }
}

# ── build_formulas() ──────────────────────────────────────────────────────────
# Returns a named list of formulas for the main model and the replicate-
# adjusted model, plus the matching emmeans specs.  Automatically adjusts
# whether to include `ug/L` interaction and/or replicate offset.
build_formulas <- function(response, data) {
  has_chem <- length(unique(data$`ug/L`)) > 1
  has_reps <- length(unique(data$replicate))      > 1

  base_terms <- if (has_chem) "`ug/L`" else "food_treatment"
  rep_terms  <- if (has_reps) paste(base_terms, "+ replicate")  else base_terms

  list(
    main          = as.formula(paste(response, "~", base_terms)),
    with_reps     = as.formula(paste(response, "~", rep_terms)),
    emmeans_main  = as.formula(paste("~", base_terms)),
    emmeans_reps  = as.formula(paste("~", rep_terms))
  )
}

# ── append_dharma_diagnostics() ───────────────────────────────────────────────
# Adds one page of DHARMa residual plots to the currently open PDF device.
# Falls back to base-R plot.glm() if DHARMa simulation fails.
append_dharma_diagnostics <- function(model, label) {
  res_sim <- try(simulateResiduals(fittedModel = model, plot = FALSE), silent = TRUE)

  op <- par(no.readonly = TRUE)
  on.exit(par(op))     # always restore graphics parameters

  if (!inherits(res_sim, "try-error")) {
    par(oma = c(0, 0, 3, 0))
    plot(res_sim)
    mtext(label, outer = TRUE, cex = 1.1)
  } else {
    # Fallback: standard GLM diagnostic quad-plot
    par(mfrow = c(2, 2))
    plot(model)
  }
}
# ── make_histogram() ──────────────────────────────────────────────────────────
# Produces a multi-panel histogram figure for a single response variable,
# with density curves overlaid on each panel.  The number of panels depends
# on whether a chemical treatment is present in the data:
#
#   has_chem = FALSE  →  2 panels: overall | by food treatment
#   has_chem = TRUE   →  4 panels: overall | by food | by chem | food × chem
#
# Panels are stacked vertically and assembled with patchwork.
# The function is called inside the histogram loop of analyze_experiment()
# and the returned plot is printed to the currently open PDF device.
#
# Binwidth is chosen automatically via the Freedman-Diaconis rule:
#   bw = 2 * IQR(x) / n^(1/3)
# which adapts to each variable's spread and sample size.  A fallback of
# range/30 is used when IQR = 0 (e.g. zero-inflated variables).
#
# Arguments:
#   data        – data frame returned by load_experiment(); must contain
#                 columns `food_treatment`, ``ug/L``, and `resp`
#   resp        – character; name of the numeric response column to plot
#   experiment  – character; experiment name used in the plot title
#
# Returns:
#   A patchwork object, or invisible(NULL) if `resp` is absent or non-numeric.
# ─────────────────────────────────────────────────────────────────────────────
make_histogram <- function(data, resp, experiment) {
  
  if (!resp %in% names(data) || !is.numeric(data[[resp]])) return(invisible(NULL))
  
  df <- data |> filter(!is.na(.data[[resp]]))
  
  bw <- 2 * IQR(df[[resp]], na.rm = TRUE) / (nrow(df)^(1/3))
  if (bw == 0) bw <- diff(range(df[[resp]], na.rm = TRUE)) / 30
  
  has_chem <- length(unique(df$`ug/L`)) > 1
  
  

  # ── Panel A: overall ──────────────────────────────────────────────────────
  p_overall <- ggplot(df, aes(x = .data[[resp]])) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth  = bw, fill = "grey70", colour = "white", linewidth = 0.3
    ) +
    geom_density(colour = "black", linewidth = 0.7) +
    labs(x = resp, y = "density", subtitle = "all groups combined") +
    theme_bw()

  # ── Panel B: by food treatment ────────────────────────────────────────────
  p_food <- ggplot(df, aes(x = .data[[resp]], fill = food_treatment)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = bw, colour = "white", linewidth = 0.3, alpha = 0.8
    ) +
    geom_density(aes(colour = food_treatment), linewidth = 0.7, show.legend = FALSE) +
    facet_wrap(~ food_treatment, ncol = 3, scales = "free_y") +
    scale_fill_manual(
      values = c("no food" = "#F4A460", "low food" = "#4169E1", "high food" = "#2E8B57"),
      guide  = "none"
    ) +
    scale_colour_manual(
      values = c("no food" = "#c47a20", "low food" = "#1a3a8f", "high food" = "#1a5c2e")
    ) +
    labs(x = resp, y = "density", subtitle = "by food treatment") +
    theme_bw()

  # ── Early return if no chem treatment ────────────────────────────────────
  if (!has_chem) {
    return(
      (p_overall / p_food) +
        plot_annotation(
          title = paste(experiment, "-", resp),
          theme = theme(plot.title = element_text(face = "bold"))
        )
    )
  }

  # ── Panel C: by chem treatment ────────────────────────────────────────────
  p_chem <- ggplot(df, aes(x = .data[[resp]], fill = `ug/L`)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = bw, colour = "white", linewidth = 0.3, alpha = 0.8
    ) +
    geom_density(aes(colour = `ug/L`), linewidth = 0.7, show.legend = FALSE) +
    facet_wrap(~ `ug/L`, scales = "free_y") +
    scale_fill_viridis_d(option = "plasma", end = 0.8, guide = "none") +
    scale_colour_viridis_d(option = "plasma", end = 0.8, guide = "none") +
    labs(x = resp, y = "density", subtitle = "by chemical treatment") +
    theme_bw()

  # ── Panel D: interaction (food × chem) ───────────────────────────────────
  p_interaction <- ggplot(df, aes(x = .data[[resp]], fill = food_treatment)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = bw, colour = "white", linewidth = 0.3, alpha = 0.8
    ) +
    geom_density(aes(colour = food_treatment), linewidth = 0.7, show.legend = FALSE) +
    facet_grid(
      `ug/L` ~ food_treatment,
      scales   = "free_y",
      labeller = labeller(
        food_treatment = c("no food" = "NF", "low food" = "LF", "high food" = "HF")
      )
    ) +
    scale_fill_manual(
      values = c("no food" = "#F4A460", "low food" = "#4169E1", "high food" = "#2E8B57"),
      guide  = "none"
    ) +
    scale_colour_manual(
      values = c("no food" = "#c47a20", "low food" = "#1a3a8f", "high food" = "#1a5c2e")
    ) +
    labs(x = resp, y = "density", subtitle = "food × chemical treatment") +
    theme_bw()

  # ── Combine all four panels ───────────────────────────────────────────────
  (p_overall / p_food / p_chem / p_interaction) +
    plot_annotation(
      title = paste(experiment, "-", resp),
      theme = theme(plot.title = element_text(face = "bold"))
    )
}



# =============================================================================
# analyze_experiment()
# =============================================================================
# For each response variable in `responses`:
#   • fits a GLM (with and without replicate offset)
#   • runs pairwise contrasts (BH-adjusted) via emmeans
#   • writes a boxplot with compact-letter-display significance labels
#   • appends DHARMa diagnostics to a combined PDF
#
# Then, if individual-droplet data exist, fits a mixed model (LMM with
# well_name as random intercept) for droplet area and produces a violin plot.
#
# Arguments:
#   data       – data frame returned by load_experiment()
#   responses  – tibble with columns `response` and `family` (default: global)
# =============================================================================
analyze_experiment <- function(data, responses = response) {

  experiment <- unique(data$experiment)
  out_dir    <- file.path("output", experiment)

  # Create output sub-directories (silently if they already exist)
  dir.create(file.path(out_dir, "graphs"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

  # Open a single multi-page PDF that will collect all diagnostic plots
  pdf(file.path(out_dir, "graphs", "histograms.pdf"), width = 12, height = 12, onefile = TRUE)
  # ── Histograms for all responses (one page each) ──────────────────────────
  for (r in responses$response) {
    if (!r %in% names(data)) next

  p_hist <- make_histogram(data, r, experiment)
  if (!is.null(p_hist)) print(p_hist)   # printed to the open PDF device

  }

  dev.off()  # close the histogram PDF before starting the diagnostics PDF
  # ── Per-response modelling loop ───────────────────────────────────────────
  pdf(file.path(out_dir, "graphs", "model_diagnostics.pdf"), width = 12, height = 6, onefile = TRUE)

  for (i in seq_len(nrow(responses))) {

    resp       <- responses$response[i]
    fam        <- responses$family[[i]]
    is_negbin  <- is.character(fam) && fam == "negbin"
    family_obj <- resolve_family(fam)

    message(sprintf("[%d/%d] Analyzing response: %s", i, nrow(responses), resp))

    # Build the four formulas needed for this response
    fmls <- build_formulas(resp, data)

    # ── Model with replicate as covariate (used for replicate-level contrasts) 
    glm_rep  <- fit_glm(fmls$with_reps, family_obj, data, is_negbin)
    anodev_rep <- car::Anova(glm_rep, type = 2)
    emm_rep  <- emmeans(glm_rep, fmls$emmeans_reps)
    pairs_rep_df <- as.data.frame(pairs(emm_rep, adjust = "BH"))

    # ── Main model (food ± chem treatment only, no replicate term) ───────────
    glm_main <- fit_glm(fmls$main, family_obj, data, is_negbin)
    anodev_main <- car::Anova(glm_main, type = 2)
    emm_main <- emmeans(glm_main, fmls$emmeans_main)
    pairs_df <- as.data.frame(pairs(emm_main, adjust = "BH"))

    # Compact letter display for annotating the plot
    cld_df <- cld(emm_main, Letters = letters, adjust = "BH") |>
      as.data.frame() |>
      filter(!is.na(.group))

    # ── Boxplot: response by food treatment, coloured by chem treatment ──────
    y_max <- max(data[[resp]], na.rm = TRUE)

    p_main <- ggplot(data, aes(x = `ug/L`, y = .data[[resp]],
                               fill = food_treatment, color = replicate)) +
      geom_boxplot(
        outlier.shape = NA,
        position      = position_dodge(width = 0.6),
        width         = 0.5,
        color         = "black"
      ) +
      geom_jitter(
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6),
        size     = 2, alpha = 0.6
      ) +
      stat_summary(
        fun      = median,
        geom     = "crossbar",
        position = position_dodge(width = 0.6),
        width    = 0.4,
        linewidth = 0.3,
        color    = "black"
      ) +
      # Compact letter display: groups that share a letter are not significantly
      # different after BH correction
      geom_text(
        data     = cld_df,
        aes(x = `ug/L`, y = y_max, label = .group),
        inherit.aes = FALSE,
        colour      = "black",
        show.legend = FALSE,
        hjust    = 0.5,
        vjust    = -0.5,
        position = position_dodge(width = 0.6),
        fontface = "bold",
        size     = 4
      ) +
      scale_x_discrete(
        labels = c("no food" = "NF", "low food" = "LF", "high food" = "HF")
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      labs(
        y     = resp,
        x     = NULL,
        title = paste(experiment, "-", resp),
        color  = "Replicate",
        subtitle = paste(paste(deparse(fmls$main), collapse = ""), "| family:", if (is_negbin) "negbin | log" else paste(family_obj$family, "|", family_obj$link))
      ) +
      scale_color_brewer(palette = "Dark2") +
      theme_bw()

    ggsave(
      file.path(out_dir, "graphs", paste0(resp, "_by_treatment.pdf")),
      p_main
    )

    # ── Save pairwise-comparison tables ──────────────────────────────────────
    write.csv(
      pairs_rep_df,
      file.path(out_dir, "results", paste0(resp, "_pairs_replicates.csv")),
      row.names = FALSE
    )
    write.csv(
      pairs_df,
      file.path(out_dir, "results", paste0(resp, "_pairs_treatment.csv")),
      row.names = FALSE
    )

    write.csv(
      anodev_rep,
      file.path(out_dir, "results", paste0(resp, "_anova_replicates.csv")),
      row.names = TRUE
    )

    write.csv(
      anodev_main,
      file.path(out_dir, "results", paste0(resp, "_anova_treatment.csv")),
      row.names = TRUE
    )

    # ── Append DHARMa diagnostics to the combined PDF ────────────────────────
    fam_label <- if (is_negbin) "negbin | log" else
      paste(family_obj$family, "|", family_obj$link)

    append_dharma_diagnostics(
      model = glm_main,
      label = paste("Response:", resp, "| Family:", fam_label)
    )

  } # end response loop
  
  dev.off()  # close the model diagnostics PDF

  # ── Individual-droplet area analysis (optional) ───────────────────────────
  # Only runs if load_experiment() populated the global `droplet_data` object
  # and it contains rows for the current experiment.
  if (!exists("droplet_data") || nrow(droplet_data) == 0) {
    cat("\n---", experiment, "done (no droplet data) ---\n")
    return(invisible(NULL))
  }

  droplet_current <- droplet_data |>
    filter(experiment == !!experiment)

  if (nrow(droplet_current) == 0) {
    cat("\n---", experiment, "done (no droplet data for this experiment) ---\n")
    return(invisible(NULL))
  }

  # Build formulas for the droplet mixed model
  fmls_drop <- build_formulas("area", droplet_current)

  # LMM: fixed effects = food (± chem) + replicate; random intercept per well
  # The random effect accounts for the non-independence of droplets from the
  # same well (pseudo-replication at the droplet level).
  lmm_drop <- glmer(
    update(fmls_drop$with_reps, . ~ . + (1 | well_name)),
    data   = droplet_current,
    family = Gamma(link = "log")
  )

  emm_drop    <- emmeans(lmm_drop, fmls_drop$emmeans_main)
  pairs_drop  <- as.data.frame(pairs(emm_drop, adjust = "BH"))
  cld_drop_df <- cld(emm_drop, Letters = letters, adjust = "BH") |>
    as.data.frame() |>
    filter(!is.na(.group))

  # ── LMM diagnostic PDF ────────────────────────────────────────────────────
  diag_path <- file.path(out_dir, "graphs", "individual_droplet_area_diagnostics.pdf")
  pdf(diag_path, width = 12, height = 6)
  
  hist(droplet_current$area, main = paste(experiment, "- droplet area"),
       xlab = "area", col = "grey80", breaks = 30)
  append_dharma_diagnostics(
    model = lmm_drop,
    label = "Response: droplet area | LMM (random intercept: well_name)"
  )
  dev.off()

  # ── Violin plot: droplet area by food treatment ───────────────────────────
  y_max_drop <- max(droplet_current$area, na.rm = TRUE)

  p_drop <- ggplot(
    droplet_current,
    aes(x = `ug/L`, y = area, fill = food_treatment)
  ) +
    geom_violin(position = position_dodge(width = 0.6), trim = FALSE) +
    stat_summary(
      fun      = median,
      geom     = "crossbar",
      position = position_dodge(width = 0.6),
      width    = 0.4,
      linewidth = 0.3
    ) +
    geom_text(
      data     = cld_drop_df,
      aes(x = `ug/L`, y = y_max_drop, label = .group),
      inherit.aes = FALSE,
      hjust    = 0.5,
      vjust    = -0.5,
      position = position_dodge(width = 0.6),
      fontface = "bold",
      size     = 4
    ) +
    scale_x_discrete(
      labels = c("no food" = "NF", "low food" = "LF", "high food" = "HF")
    ) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      y     = "area (log10)",
      x     = NULL,
      colour   = "Replicate",
      title = paste(experiment, "- droplet area by treatment"),
      subtitle = paste(paste(deparse(fmls$main), collapse = ""), "| family:", if (is_negbin) "negbin | log" else paste(family_obj$family, "|", family_obj$link))
    ) +
    theme_bw()

  ggsave(
    file.path(out_dir, "graphs", "individual_droplets_by_treatment.pdf"),
    p_drop
  )

  write.csv(
    pairs_drop,
    file.path(out_dir, "results", "individual_droplet_area_pairs_treatment.csv"),
    row.names = FALSE
  )

  cat("\n---", experiment, "done ---\n")
}


# =============================================================================
# Entry point
# =============================================================================
# Prompt the user to choose one experiment or process all at once.
# =============================================================================
all_experiments <- list.dirs("data/TBT-CL chronic mothers/", full.names = FALSE, recursive = FALSE)
all_experiments <- all_experiments[!grepl("pilot", all_experiments, ignore.case = TRUE)] # Exclude pilot folders if present

run_experiment <- select.list(
  c("ALL", all_experiments),
  title = "Select experiment to run"
)

experiments <- if (run_experiment == "ALL") all_experiments else run_experiment


# Load all chosen experiments, then analyse each one
all_data <- map(experiments, load_experiment)
walk(all_data, analyze_experiment)