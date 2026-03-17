library(tidyverse)
library(readxl)
library(multcomp)
library(lme4)

if (!require(emmeans)) {
  install.packages("emmeans")
  library(emmeans)
}

load_experiment <- function(experiment) {
  
  data_dir <- file.path("data", experiment)
  
  # Load well labels
  well_files <- list.files(data_dir, pattern = paste0(".*_wells\\.xlsx"), full.names = TRUE)
  treatment <- map_dfr(well_files, function(f) {
    rep_id <- sub(".*measures_export_(.*)\\_wells\\.xlsx", "\\1", basename(f))
    df <- read_excel(f, col_names = TRUE)
    df$replicate <- rep_id  
    df
  })
  treatment$replicate <- as.factor(treatment$replicate)
  
  # Load data
  data_files <- list.files(data_dir, pattern = paste0(".*\\.csv"), full.names = TRUE)
  data <- map_dfr(data_files, function(f) {
    rep_id <- sub("measures_export_(\\d+)\\.csv", "\\1", basename(f))
    df <- read.csv(f)
    df$replicate <- rep_id  
    df$experiment <- experiment
    df
  })
  data$replicate <- as.factor(data$replicate)
  
  # Join treatment labels
  data <- data |>
    mutate(well_letter = substr(well_name, 1, 1)) |>
    inner_join(treatment, by = c("well_letter" = "well", "replicate"))
  
  data$food_treatment <- factor(data$food_treatment, levels = c("no food", "low food", "high food"))
  data$chem_treatment <- factor(data$chem_treatment)

  # Normalize Number of droplets by Daphnia area
  data <- data |>
    mutate(num_droplets_per_area = num_droplets / daphnia_area)


  data <- data |>
    rowwise() |>
    mutate(zero_count = sum(c_across(where(is.numeric)) == 0, na.rm = TRUE)) |>
    ungroup()

  # Identify rows to be removed
  removed_rows <- data |>
    filter(zero_count >= 5)

  # Issue warning if any rows are removed
  if (nrow(removed_rows) > 0) {
    warning(
      paste0(
        "Rows removed due to >=5 zeros:\n",
        paste0("Well ",removed_rows$well_name, " (", removed_rows$zero_count, " zeros)", collapse = "\n")
      )
    )
  }

  data <- data |>
  dplyr::filter(zero_count <= 4) |>
  dplyr::select(-zero_count)
  
# Apply filtering and drop helper column
  
  
  # Load and enrich droplet data
  masks_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  masks_dirs <- masks_dirs[grepl("masks_", basename(masks_dirs))]
  
  droplet_data <- map_dfr(masks_dirs, function(masks_dir) {
    rep_id <- sub(".*masks_(\\d+)", "\\1", basename(masks_dir))
    mask_files <- list.files(masks_dir, pattern = "_droplets\\.csv$", full.names = TRUE)
    
    map_dfr(mask_files, function(f) {
      df <- tryCatch(
        read.csv(f),
        error = function(e) data.frame()
      )
      
      # Skip empty files
      if (nrow(df) == 0) {
        return(data.frame())
      }
      
      # Extract well name from filename
      well_name <- sub("^([A-Z]\\d+)_.*", "\\1", basename(f))
      df$well_name <- well_name
      df$replicate <- rep_id
      df$experiment <- experiment
      df
    })
  })
  
  # Join treatment information to droplet data
  if (nrow(droplet_data) > 0) {
    droplet_data <- droplet_data |>
      mutate(well_letter = substr(well_name, 1, 1), replicate = as.factor(replicate)) |>
      inner_join(treatment, by = c("well_letter" = "well", "replicate"))
    
    droplet_data$food_treatment <- factor(droplet_data$food_treatment, levels = c("no food", "low food", "high food"))
    droplet_data$chem_treatment <- factor(droplet_data$chem_treatment)
    
    # Attach droplet_data to the environment for access after load_experiment
    assign("droplet_data", droplet_data, envir = .GlobalEnv)
  }

  return(data)
}



response <- tibble::tibble(
  response = c("droplet_intensity_total", "droplet_area_total", "daphnia_size", "num_droplets", "num_droplets_per_area"),
  family   = list(Gamma(link = "log"), gaussian(), gaussian(), "negbin", Gamma(link = "log"))
)

analyze_experiment <- function(data, responses = response) {

  # prepare experiment output dirs and a combined diagnostics PDF
  experiment <- unique(data$experiment)
  out_dir <- file.path("output", experiment)
  dir.create(file.path(out_dir, "graphs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

  diagnostics_pdf <- file.path(out_dir, "graphs", paste0("responses_and_diagnostics.pdf"))
  pdf(diagnostics_pdf, width = 12, height = 6, onefile = TRUE)

  # histograms for all requested responses (one page each)
  for (r in responses$response) {
    if (r %in% names(data) && is.numeric(data[[r]])) {
      hist(data[[r]], main = paste(experiment, "-", r), xlab = r, col = "grey", breaks = 30)
    } else if (r %in% names(data)) {
      counts <- table(data[[r]])
      barplot(counts, main = paste(experiment, "-", r), xlab = r)
    }
  }

  for (i in seq_len(nrow(responses))) {
    resp <- responses$response[i]
    fam  <- responses$family[[i]]   # use [[ ]] to get the actual family object

    # Check for negbin before any resolution attempt
  is_negbin <- (is.character(fam) && fam == "negbin")

  # derive family_obj: accept family objects, functions, or character expressions
  if (is_negbin) {
    family_obj <- "negbin"  # keep as sentinel; fit_glm will handle it
  } else if (is.character(fam)) {
    family_obj <- tryCatch(
      match.fun(fam)(),
      error = function(e) {
        tryCatch(eval(parse(text = fam)), error = function(e2) stop("Could not resolve family: ", fam))
      }
    )
  } else if (inherits(fam, "family")) {
    family_obj <- fam
  } else if (is.function(fam)) {
    family_obj <- fam()
  } else {
    family_obj <- fam
  }
      
    # build formulas depending on chem_treatment presence
    has_chem  <- length(unique(data$chem_treatment)) > 1
    has_reps  <- length(unique(data$replicate)) > 1

    base_terms    <- if (has_chem) "food_treatment * chem_treatment" else "food_treatment"
    rep_terms     <- if (has_reps) paste(base_terms, "+ replicate")  else base_terms

    formula                  <- as.formula(paste(resp, "~", base_terms))
    emmeans_formula          <- as.formula(paste("~", base_terms))
    formula_replicates       <- as.formula(paste(resp, "~", rep_terms))
    emmeans_formula_replicates <- as.formula(paste("~", rep_terms))

    fit_glm <- function(formula, family_obj, data, is_negbin = FALSE) {
      if (is_negbin) {
        if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
        MASS::glm.nb(formula, data = data)
      } else {
        glm(formula, family = family_obj, data = data)
      }
    }


    glm_replicates <- fit_glm(formula_replicates, family_obj, data, is_negbin)
    anodev_replicates <- anova(glm_replicates, test = "F")
    emmeans_replicates <- emmeans(glm_replicates, emmeans_formula_replicates)
    pairs_comparisons_replicates <- pairs(emmeans_replicates, adjust = "BH")
    pairs_comparisons_replicates_df <- as.data.frame(pairs_comparisons_replicates)

    experiment <- unique(data$experiment)
    out_dir <- file.path("output", experiment)
    dir.create(file.path(out_dir, "graphs"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
    
    # Plot by replicate: compute ranges from the chosen response
    y_max <- max(data[[resp]], na.rm = TRUE)
    y_range <- diff(range(data[[resp]], na.rm = TRUE))

    # sig_data_rep <- pairs_comparisons_replicates_df %>%
    #   separate(contrast, into = c("group1_full", "group2_full"), sep = " - ", remove = FALSE) %>%
    #   mutate(
    #     treatment1 = str_extract(group1_full, "^(.+?)(?= replicate)", group = 1),
    #     replicate1  = str_extract(group1_full, "replicate(\\d+)$", group = 1),
    #     treatment2 = str_extract(group2_full, "^(.+?)(?= replicate)", group = 1),
    #     replicate2  = str_extract(group2_full, "replicate(\\d+)$", group = 1)
    #   ) %>%
    #   filter(treatment1 == treatment2) %>%
    #   rename(food_treatment = treatment1) %>%
    #   mutate(stars = case_when(
    #     p.value < 0.001 ~ "***",
    #     p.value < 0.01  ~ "**",
    #     p.value < 0.05  ~ "*",
    #     TRUE            ~ "n.s."
    #   )) %>%
    #   mutate(
    #     replicate1 = factor(replicate1, levels = levels(data$replicate)),
    #     replicate2 = factor(replicate2, levels = levels(data$replicate)),
    #     x1_num = as.numeric(replicate1),
    #     x2_num = as.numeric(replicate2)
    #   ) %>%
    #   group_by(food_treatment) %>%
    #   mutate(y_position = y_max + y_range * 0.08 * row_number()) %>%
    #   ungroup()

    # p1 <- ggplot(data, aes(x = replicate, y = .data[[resp]], color = food_treatment)) +
    #   geom_jitter(width = 0.1, size = 2) +
    #   stat_summary(fun = median, geom = "crossbar", width = 0.5,
    #                color = "black", linewidth = 0.3) +
    #   geom_segment(
    #     data = sig_data_rep,
    #     aes(x = x1_num, xend = x2_num, y = y_position, yend = y_position),
    #     inherit.aes = FALSE, color = "black"
    #   ) +
    #   geom_segment(
    #     data = sig_data_rep,
    #     aes(x = x1_num, xend = x1_num, y = y_position, yend = y_position - y_range * 0.02),
    #     inherit.aes = FALSE, color = "black"
    #   ) +
    #   geom_segment(
    #     data = sig_data_rep,
    #     aes(x = x2_num, xend = x2_num, y = y_position, yend = y_position - y_range * 0.02),
    #     inherit.aes = FALSE, color = "black"
    #   ) +
    #   geom_text(
    #     data = sig_data_rep,
    #     aes(x = (x1_num + x2_num) / 2, y = y_position + y_range * 0.03, label = stars),
    #     inherit.aes = FALSE, color = "black", size = 4
    #   ) +
    #   facet_wrap(~ factor(food_treatment, levels = c("no food", "low food", "high food"))) +
    #   scale_color_manual(values = c("no food" = "#F4A460",
    #                                 "low food" = "#4169E1",
    #                                 "high food" = "#2E8B57")) +
    #   scale_x_discrete(labels = c("1" = "I", "2" = "II", "3" = "III")) +
    #   scale_y_continuous(expand = expansion(mult = c(0.05, 0.02))) +
    #   labs(y = resp, x = NULL, title = paste(experiment, "-", resp)) +
    #   theme_bw()

    # ggsave(file.path(out_dir, "graphs", paste0(resp, "_by_replicate.pdf")), p1)
    
    glm <- fit_glm(formula, family_obj, data, is_negbin)
    anodev <- anova(glm, test = "F")
    emmeans_obj <- emmeans(glm, emmeans_formula)
    pairs_comparisons <- as.data.frame(pairs(emmeans_obj, adjust = "BH"))
    cld_df <- cld(emmeans_obj, Letters = letters, adjust = "BH") |>
      as.data.frame() |>
      filter(!is.na(.group))

    y_max2 <- max(data[[resp]], na.rm = TRUE)
    y_range2 <- diff(range(data[[resp]], na.rm = TRUE))

    treatment_levels <- levels(data$food_treatment)

    # sig_data2 <- pairs_comparisons %>%
    #   separate(contrast, into = c("group1_full", "group2_full"), sep = " - ", remove = FALSE) %>%
    #   mutate(
    #     stars = case_when(
    #       p.value < 0.001 ~ "***",
    #       p.value < 0.01  ~ "**",
    #       p.value < 0.05  ~ "*",
    #       TRUE            ~ "n.s."
    #     ),
    #     group1_full = factor(group1_full, levels = treatment_levels),
    #     group2_full = factor(group2_full, levels = treatment_levels),
    #     x1_num = as.numeric(group1_full),
    #     x2_num = as.numeric(group2_full)
    #   ) %>%
    #   group_by() %>%
    #   mutate(y_position = y_max2 + y_range2 * 0.08 * row_number()) %>%
    #   ungroup()

    # bracket_step2 <- y_range2 * 0.08
    # tick_length2  <- bracket_step2 * 0.25

    p2 <- ggplot(data, aes(x = food_treatment, y = .data[[resp]], fill = chem_treatment)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.6), width = 0.5) +
      geom_jitter(
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6),
        size = 2
      ) +
      stat_summary(
        fun = median,
        geom = "crossbar",
        position = position_dodge(width = 0.6),
        width = 0.4,
        linewidth = 0.3
      ) +    
      # geom_segment(
      #   data = sig_data2,
      #   aes(x = x1_num, xend = x2_num, y = y_position, yend = y_position),
      #   inherit.aes = FALSE, color = "black"
      # ) +
      # geom_segment(
      #   data = sig_data2,
      #   aes(x = x1_num, xend = x1_num, y = y_position, yend = y_position - tick_length2),
      #   inherit.aes = FALSE, color = "black"
      # ) +
      # geom_segment(
      #   data = sig_data2,
      #   aes(x = x2_num, xend = x2_num, y = y_position, yend = y_position - tick_length2),
      #   inherit.aes = FALSE, color = "black"
      # ) +
      # geom_text(
      #   data = sig_data2,
      #   aes(x = (x1_num + x2_num) / 2, y = y_position + y_range2 * 0.03, label = stars),
      #   inherit.aes = FALSE, color = "black", size = 4
      # ) +
      geom_text(
        data = cld_df,
        aes(x = food_treatment, y = y_max2, label = .group),
        hjust = 0.5,  # Center over box
        vjust = -0.5, # Slight upward offset
        position = position_dodge(width = 0.6),  # Must match boxplot
        fontface = "bold",
        size = 4
    ) +
      scale_x_discrete(labels = c("no food" = "NF", "low food" = "LF", "high food" = "HF")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      labs(y = resp, x = NULL, title = paste(experiment, "-", resp)) +
      theme_bw()


    ggsave(file.path(out_dir, "graphs", paste0(resp, "_by_treatment.pdf")), p2)

    # save emmeans pairwise results for this response
    write.csv(pairs_comparisons_replicates_df, file.path(out_dir, "results", paste0(resp, "_pairs_replicates.csv")), row.names = FALSE)
    write.csv(pairs_comparisons, file.path(out_dir, "results", paste0(resp, "_pairs_treatment.csv")), row.names = FALSE)

    # append GLM diagnostic plots to the combined PDF using DHARMa
    if (exists("glm") && inherits(glm, "glm")) {
      if (!require(DHARMa, quietly = TRUE)) {
        install.packages("DHARMa")
        library(DHARMa)
      }
      res_sim <- try(simulateResiduals(fittedModel = glm, plot = FALSE), silent = TRUE)
      if (!inherits(res_sim, "try-error")) {
        op <- par(no.readonly = TRUE)
        par(oma = c(0, 0, 3, 0))
        try(plot(res_sim))
        fam_label <- if (is.character(family_obj) && family_obj == "negbin") {
          "negbin | log"
        } else {
          paste(family_obj$family, "|", family_obj$link)
        }
        mtext(paste("Response:", resp, "| Family:", fam_label), outer = TRUE, cex = 1.1)
        par(op)
      } else {
        # fallback to basic diagnostic plots if DHARMa simulation fails
        op <- par(no.readonly = TRUE)
        par(mfrow = c(2, 2))
        try(plot(glm))
        par(op)
      }
    }

  } # end for resp

  # close diagnostics PDF
  dev.off()

  # Analyze droplet data if available
  if (exists("droplet_data") && nrow(droplet_data) > 0) {
    
    # Filter droplet data to current experiment
    droplet_current <- droplet_data |>
      filter(experiment == unique(data$experiment))
    
    if (nrow(droplet_current) > 0) {
      
      # GLM analysis for droplet area by replicate
      
      has_chem  <- length(unique(droplet_data$chem_treatment)) > 1
      has_reps  <- length(unique(droplet_data$replicate)) > 1

      base_terms    <- if (has_chem) "food_treatment * chem_treatment" else "food_treatment"
      rep_terms     <- if (has_reps) paste(base_terms, "+ replicate")  else base_terms

      formula_droplets <- as.formula(paste("area ~", base_terms))
      emmeans_formula_droplets <- as.formula(paste("~", base_terms))
      
      formula_droplets_rep <- as.formula(paste("area ~", rep_terms))
      emmeans_formula_droplets_rep <- as.formula(paste("~", rep_terms))
      
      
      
#       glm_droplets_rep <- lmer(
#   update(formula_droplets_rep, . ~ . + (1 | well_name)),
#   data = droplet_current
# )
#       emmeans_droplets_rep <- emmeans(glm_droplets_rep, emmeans_formula_droplets_rep)
#       pairs_droplets_rep <- pairs(emmeans_droplets_rep, adjust = "BH")
#       pairs_droplets_rep_df <- as.data.frame(pairs_droplets_rep)
#       cld_droplets_rep_df <- cld(emmeans_droplets_rep, Letters = letters, adjust = "BH") |>
#         as.data.frame() |>
#         filter(!is.na(.group))
      
#       # Plot droplet area by replicate with significance
#       y_max_drop <- max(droplet_current$area, na.rm = TRUE)
#       y_range_drop <- diff(range(droplet_current$area, na.rm = TRUE))
      
      # sig_data_drop_rep <- pairs_droplets_rep_df %>%
      #   separate(contrast, into = c("group1_full", "group2_full"), sep = " - ", remove = FALSE) %>%
      #   mutate(
      #     treatment1 = str_extract(group1_full, "^(.+?)(?= replicate)", group = 1),
      #     replicate1  = str_extract(group1_full, "replicate(\\d+)$", group = 1),
      #     treatment2 = str_extract(group2_full, "^(.+?)(?= replicate)", group = 1),
      #     replicate2  = str_extract(group2_full, "replicate(\\d+)$", group = 1)
      #   ) %>%
      #   filter(treatment1 == treatment2) %>%
      #   rename(food_treatment = treatment1) %>%
      #   mutate(stars = case_when(
      #     p.value < 0.001 ~ "***",
      #     p.value < 0.01  ~ "**",
      #     p.value < 0.05  ~ "*",
      #     TRUE            ~ "n.s."
      #   )) %>%
      #   mutate(
      #     replicate1 = factor(replicate1, levels = levels(droplet_current$replicate)),
      #     replicate2 = factor(replicate2, levels = levels(droplet_current$replicate)),
      #     x1_num = as.numeric(replicate1),
      #     x2_num = as.numeric(replicate2)
      #   ) %>%
      #   group_by(food_treatment) %>%
      #   mutate(y_position = y_max_drop * (1.5 ^ row_number())) %>%
      #   ungroup()
      
      # p_droplet_rep <- ggplot(droplet_current, aes(x = replicate, y = area, color = food_treatment)) +
      #   geom_violin() +
      #   stat_summary(fun = median, geom = "crossbar", width = 0.5,
      #                color = "black", linewidth = 0.3) +
      #   geom_segment(
      #     data = sig_data_drop_rep,
      #     aes(x = x1_num, xend = x2_num, y = y_position, yend = y_position),
      #     inherit.aes = FALSE, color = "black"
      #   ) +
      #   geom_segment(
      #     data = sig_data_drop_rep,
      #     aes(x = x1_num, xend = x1_num, y = y_position, yend = y_position - y_range_drop * 0.02),
      #     inherit.aes = FALSE, color = "black"
      #   ) +
      #   geom_segment(
      #     data = sig_data_drop_rep,
      #     aes(x = x2_num, xend = x2_num, y = y_position, yend = y_position - y_range_drop * 0.02),
      #     inherit.aes = FALSE, color = "black"
      #   ) +
      #   geom_text(
      #     data = sig_data_drop_rep,
      #     aes(x = (x1_num + x2_num) / 2, y = y_position + y_range_drop * 0.3, label = stars),
      #     inherit.aes = FALSE, color = "black", size = 4
      #   ) +
      #   facet_wrap(~ factor(food_treatment, levels = c("no food", "low food", "high food"))) +
      #   scale_color_manual(values = c("no food" = "#F4A460",
      #                                 "low food" = "#4169E1",
      #                                 "high food" = "#2E8B57")) +
      #   scale_x_discrete(labels = c("1" = "I", "2" = "II", "3" = "III")) +
      #   scale_y_log10(expand = expansion(mult = c(0.05, 0.02))) +
      #   labs(y = "area (log10)", x = NULL, title = paste(unique(data$experiment), "- droplet area by replicate")) +
      #   theme_bw()
      
      # ggsave(file.path(out_dir, "graphs", "individual_droplets_by_replicate.pdf"), p_droplet_rep)
      
      # GLM analysis for droplet area by treatment
      glm_droplets <- lmer(
                      update(formula_droplets_rep, . ~ . + (1 | well_name)),
                      data = droplet_current
                      )
      emmeans_droplets <- emmeans(glm_droplets, emmeans_formula_droplets)
      pairs_droplets <- as.data.frame(pairs(emmeans_droplets, adjust = "BH"))
      cld_droplets_df <- cld(emmeans_droplets, Letters = letters, adjust = "BH") |>
        as.data.frame() |>
        filter(!is.na(.group))
      
      y_max_drop2 <- max(droplet_current$area, na.rm = TRUE)
      y_range_drop2 <- diff(range(droplet_current$area, na.rm = TRUE))
      
      pdf(file.path(out_dir, "graphs", "Individual droplet_area_diagnostics.pdf"), width = 12, height = 6)
      # append GLM diagnostic plots to the combined PDF using DHARMa
      if (exists("glm_droplets") && inherits(glm_droplets, "lmerMod")) {
        if (!require(DHARMa, quietly = TRUE)) {
          install.packages("DHARMa")
          library(DHARMa)
        }
        res_sim <- try(simulateResiduals(fittedModel = glm_droplets, plot = FALSE), silent = TRUE)
        if (!inherits(res_sim, "try-error")) {
          op <- par(no.readonly = TRUE)
          par(oma = c(0, 0, 3, 0))
          try(plot(res_sim))
          mtext(paste("Response: droplet area | Mixed Model Diagnostics"), outer = TRUE, cex = 1.1)
          par(op)
        } else {
          # fallback to basic diagnostic plots if DHARMa simulation fails
          op <- par(no.readonly = TRUE)
          par(mfrow = c(2, 2))
          try(plot(glm_droplets))
          par(op)
        }
      }
      dev.off()


      # treatment_levels_drop <- levels(droplet_current$food_treatment)
      
      # sig_data_drop2 <- pairs_droplets %>%
      #   separate(contrast, into = c("group1_full", "group2_full"), sep = " - ", remove = FALSE) %>%
      #   mutate(
      #     stars = case_when(
      #       p.value < 0.001 ~ "***",
      #       p.value < 0.01  ~ "**",
      #       p.value < 0.05  ~ "*",
      #       TRUE            ~ "n.s."
      #     ),
      #     group1_full = factor(group1_full, levels = treatment_levels_drop),
      #     group2_full = factor(group2_full, levels = treatment_levels_drop),
      #     x1_num = as.numeric(group1_full),
      #     x2_num = as.numeric(group2_full)
      #   ) %>%
      #   group_by() %>%
      #   mutate(y_position = y_max_drop2 * (1.5 ^ row_number())) %>%
      #   ungroup()
      
      # bracket_step_drop <- y_range_drop2 * 0.08
      # tick_length_drop  <- bracket_step_drop * 0.25
      
      summary_stats_drop <- droplet_current |>
        group_by(food_treatment) |>
        summarize(
          n = n(),
          median_area = median(area),
          mean_area = mean(area),
          .groups = "drop"
        )
      
      p_droplet_treat <- ggplot(droplet_current, aes(x = food_treatment, y = area, fill = chem_treatment)) +
        geom_violin(position = position_dodge(width = 0.6), trim = FALSE) +
        stat_summary(
        fun = median,
        geom = "crossbar",
        position = position_dodge(width = 0.6),
        width = 0.4,
        linewidth = 0.3
        ) +
        # geom_segment(
        #   data = sig_data_drop2,
        #   aes(x = x1_num, xend = x2_num, y = y_position, yend = y_position),
        #   inherit.aes = FALSE, color = "black"
        # ) +
        # geom_segment(
        #   data = sig_data_drop2,
        #   aes(x = x1_num, xend = x1_num, y = y_position, yend = y_position - tick_length_drop),
        #   inherit.aes = FALSE, color = "black"
        # ) +
        # geom_segment(
        #   data = sig_data_drop2,
        #   aes(x = x2_num, xend = x2_num, y = y_position, yend = y_position - tick_length_drop),
        #   inherit.aes = FALSE, color = "black"
        # ) +
        # geom_text(
        #   data = sig_data_drop2,
        #   aes(x = (x1_num + x2_num) / 2, y = y_position * 1.2, label = stars),
        #   inherit.aes = FALSE, color = "black", size = 4
        # ) +
        # geom_text(
        #   data = summary_stats_drop,
        #   aes(
        #     x = food_treatment,
        #     y = Inf,
        #     label = paste0(
        #       "n = ", n, "\n",
        #       "Med = ", round(median_area, 1), "\n",
        #       "Mean = ", round(mean_area, 1)
        #     )
        #   ),
        #   vjust = 1.2,
        #   size = 3
        # ) +
        geom_text(
          data = cld_droplets_df,
          aes(x = food_treatment, y = y_max_drop2, label = .group),
          hjust = 0.5,  # Center over box
          vjust = -0.5, # Slight upward offset
          position = position_dodge(width = 0.6),  # Must match violin
          fontface = "bold",
          size = 4
        ) +
        scale_x_discrete(labels = c("no food" = "NF", "low food" = "LF", "high food" = "HF")) +
        scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) +
        labs(y = "area (log10)", x = NULL, title = paste(unique(data$experiment), "- droplet area by treatment")) +
        theme_bw()
      
      ggsave(file.path(out_dir, "graphs", "individual_droplets_by_treatment.pdf"), p_droplet_treat)
      
      # Save results
      # write.csv(pairs_droplets_rep_df, file.path(out_dir, "results", "individual_droplet_area_pairs_replicates.csv"), row.names = FALSE)
      write.csv(pairs_droplets, file.path(out_dir, "results", "individual_droplet_area_pairs_treatment.csv"), row.names = FALSE)
    }
  }

  cat("\n---", unique(data$experiment), "done ---\n")

}

# Run 

run_experiment <- select.list(
  c("ALL", list.dirs("data", full.names = FALSE, recursive = FALSE)),
  title = "Select experiment to run"
)

experiments <- if (run_experiment == "ALL") {
  list.dirs("data", full.names = FALSE, recursive = FALSE)
} else {
  run_experiment
}

# experiments <- list.dirs("data", full.names = FALSE, recursive = FALSE)

all_data <- map(experiments, load_experiment)
walk(all_data, analyze_experiment)
