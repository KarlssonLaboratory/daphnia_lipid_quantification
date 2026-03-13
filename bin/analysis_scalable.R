library(tidyverse)
library(readxl)

if (!require(emmeans)) {
  install.packages("emmeans")
  library(emmeans)
}

load_experiment <- function(experiment) {
  
  data_dir <- file.path("data", experiment)
  
  # Load well labels
  well_files <- list.files(data_dir, pattern = paste0(experiment, ".*_wells\\.xlsx"), full.names = TRUE)
  treatment <- map_dfr(well_files, function(f) {
    rep_id <- sub(paste0(".*", experiment, "_(.*)\\_wells\\.xlsx"), "\\1", basename(f))
    df <- read_excel(f, col_names = TRUE)
    df$replicate <- rep_id  
    df
  })
  treatment$replicate <- as.factor(treatment$replicate)
  
  # Load data
  data_files <- list.files(data_dir, pattern = paste0(experiment, ".*\\.csv"), full.names = TRUE)
  data <- map_dfr(data_files, function(f) {
    rep_id <- sub(paste0(".*", experiment, "_(.*)\\.csv"), "\\1", basename(f))
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


  return(data)
}

response <- tibble::tibble(
  response = c("droplet_intensity_total", "droplet_area_total", "daphnia_size"),
  family   = list(gaussian(), gaussian(), gaussian(link = "log"))
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

    # derive family_obj: accept family objects, functions, or character expressions
    if (is.character(fam)) {
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
    if (length(unique(data$chem_treatment)) > 1) {
      formula <- as.formula(paste(resp, "~ food_treatment * chem_treatment"))
      emmeans_formula <- as.formula("~ food_treatment * chem_treatment")

      formula_replicates <- as.formula(paste(resp, "~ food_treatment * chem_treatment + replicate"))
      emmeans_formula_replicates <- as.formula("~ food_treatment * chem_treatment + replicate")
    } else {
      formula <- as.formula(paste(resp, "~ food_treatment"))
      emmeans_formula <- as.formula("~ food_treatment")

      formula_replicates <- as.formula(paste(resp, "~ food_treatment + replicate"))
      emmeans_formula_replicates <- as.formula("~ food_treatment + replicate")
    }

    glm_replicates <- glm(formula_replicates, family = family_obj, data = data)
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

    sig_data_rep <- pairs_comparisons_replicates_df %>%
      separate(contrast, into = c("group1_full", "group2_full"), sep = " - ", remove = FALSE) %>%
      mutate(
        treatment1 = str_extract(group1_full, "^(.+?)(?= replicate)", group = 1),
        replicate1  = str_extract(group1_full, "replicate(\\d+)$", group = 1),
        treatment2 = str_extract(group2_full, "^(.+?)(?= replicate)", group = 1),
        replicate2  = str_extract(group2_full, "replicate(\\d+)$", group = 1)
      ) %>%
      filter(treatment1 == treatment2) %>%
      rename(food_treatment = treatment1) %>%
      mutate(stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "n.s."
      )) %>%
      mutate(
        replicate1 = factor(replicate1, levels = levels(data$replicate)),
        replicate2 = factor(replicate2, levels = levels(data$replicate)),
        x1_num = as.numeric(replicate1),
        x2_num = as.numeric(replicate2)
      ) %>%
      group_by(food_treatment) %>%
      mutate(y_position = y_max + y_range * 0.08 * row_number()) %>%
      ungroup()

    p1 <- ggplot(data, aes(x = replicate, y = .data[[resp]], color = food_treatment)) +
      geom_jitter(width = 0.1, size = 2) +
      stat_summary(fun = median, geom = "crossbar", width = 0.5,
                   color = "black", linewidth = 0.3) +
      geom_segment(
        data = sig_data_rep,
        aes(x = x1_num, xend = x2_num, y = y_position, yend = y_position),
        inherit.aes = FALSE, color = "black"
      ) +
      geom_segment(
        data = sig_data_rep,
        aes(x = x1_num, xend = x1_num, y = y_position, yend = y_position - y_range * 0.02),
        inherit.aes = FALSE, color = "black"
      ) +
      geom_segment(
        data = sig_data_rep,
        aes(x = x2_num, xend = x2_num, y = y_position, yend = y_position - y_range * 0.02),
        inherit.aes = FALSE, color = "black"
      ) +
      geom_text(
        data = sig_data_rep,
        aes(x = (x1_num + x2_num) / 2, y = y_position + y_range * 0.03, label = stars),
        inherit.aes = FALSE, color = "black", size = 4
      ) +
      facet_wrap(~ food_treatment) +
      scale_color_manual(values = c("no food" = "#F4A460",
                                    "low food" = "#4169E1",
                                    "high food" = "#2E8B57")) +
      scale_x_discrete(labels = c("1" = "I", "2" = "II", "3" = "III")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.02))) +
      labs(y = resp, x = NULL, title = paste(experiment, "-", resp)) +
      theme_bw()

    ggsave(file.path(out_dir, "graphs", paste0(resp, "_by_replicate.pdf")), p1)
    
    glm <- glm(formula, family = family_obj, data = data)
    anodev <- anova(glm, test = "F")
    emmeans_obj <- emmeans(glm, emmeans_formula)
    pairs_comparisons <- as.data.frame(pairs(emmeans_obj, adjust = "BH"))

    y_max2 <- max(data[[resp]], na.rm = TRUE)
    y_range2 <- diff(range(data[[resp]], na.rm = TRUE))

    treatment_levels <- levels(data$food_treatment)

    sig_data2 <- pairs_comparisons %>%
      separate(contrast, into = c("group1_full", "group2_full"), sep = " - ", remove = FALSE) %>%
      mutate(
        stars = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01  ~ "**",
          p.value < 0.05  ~ "*",
          TRUE            ~ "n.s."
        ),
        group1_full = factor(group1_full, levels = treatment_levels),
        group2_full = factor(group2_full, levels = treatment_levels),
        x1_num = as.numeric(group1_full),
        x2_num = as.numeric(group2_full)
      ) %>%
      group_by() %>%
      mutate(y_position = y_max2 + y_range2 * 0.08 * row_number()) %>%
      ungroup()

    bracket_step2 <- y_range2 * 0.08
    tick_length2  <- bracket_step2 * 0.25

    p2 <- ggplot(data, aes(x = food_treatment, y = .data[[resp]])) +
      geom_jitter(width = 0.1, size = 2) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.3) +
      geom_segment(
        data = sig_data2,
        aes(x = x1_num, xend = x2_num, y = y_position, yend = y_position),
        inherit.aes = FALSE, color = "black"
      ) +
      geom_segment(
        data = sig_data2,
        aes(x = x1_num, xend = x1_num, y = y_position, yend = y_position - tick_length2),
        inherit.aes = FALSE, color = "black"
      ) +
      geom_segment(
        data = sig_data2,
        aes(x = x2_num, xend = x2_num, y = y_position, yend = y_position - tick_length2),
        inherit.aes = FALSE, color = "black"
      ) +
      geom_text(
        data = sig_data2,
        aes(x = (x1_num + x2_num) / 2, y = y_position + y_range2 * 0.03, label = stars),
        inherit.aes = FALSE, color = "black", size = 4
      ) +
      scale_x_discrete(labels = c("no food" = "NF", "low food" = "LF", "high food" = "HF")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.02))) +
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
        mtext(paste("Response:", resp, "| Family:", family_obj$family, "| Link:", family_obj$link), outer = TRUE, cex = 1.1)
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
