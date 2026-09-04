# Manuscript figures for the completed heterogeneous-information simulation.
# This script is visualization-only: it reads the frozen performance summary
# and does not alter or recompute simulation estimates, truths, or core metrics.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

CODE_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
STUDY_ROOT <- normalizePath(file.path(CODE_DIR, ".."), winslash = "/", mustWork = TRUE)
SUMMARY_DIR <- file.path(STUDY_ROOT, "03_results_summary")
FIGURE_DIR <- file.path(STUDY_ROOT, "05_figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

RUN_TAG <- Sys.getenv("EACC_RUN_TAG", "FORMAL_V3_R200_B200_7x3_20260819")
INPUT_FILE <- Sys.getenv(
  "EACC_PERFORMANCE_SUMMARY",
  file.path(SUMMARY_DIR, paste0(RUN_TAG, "_performance_summary.csv"))
)
if (!file.exists(INPUT_FILE)) stop("Performance summary not found: ", INPUT_FILE)

raw <- readr::read_csv(INPUT_FILE, show_col_types = FALSE)
required <- c(
  "scenario_id", "outcome", "setting", "architecture", "target_label",
  "target_delta", "method", "absolute_bias", "rmse", "coverage"
)
if (!all(required %in% names(raw))) stop("Performance summary is missing required columns")

architecture_labels <- c(
  A0_Full = "Complete information",
  A1_Common75 = "Common 75% availability",
  A2_Balanced75 = "Heterogeneous 75% availability",
  A3_Structured75 = "Structured 75% availability",
  A4_Common50 = "Common 50% availability",
  A5_Balanced50 = "Heterogeneous 50% availability",
  A6_Structured50 = "Structured 50% availability"
)
architecture_order <- unname(architecture_labels)
target_labels <- c(D0 = "No shift", D0p5 = "0.5 SD", D0p8 = "0.8 SD", D1p2 = "1.2 SD")
target_order <- unname(target_labels)

plot_base <- raw %>%
  mutate(
    synthesis = recode(setting, pairwise = "Pairwise meta-analysis", nma = "Network meta-analysis"),
    information_architecture = factor(unname(architecture_labels[architecture]), levels = architecture_order),
    target_population_shift = factor(unname(target_labels[target_label]), levels = target_order),
    method_group = if_else(method %in% c("Meta", "NMA"), "Standard", "EACC")
  )

# One row for every outcome x synthesis x architecture x target x method.
stopifnot(nrow(plot_base) == 224L)
cell_counts <- plot_base %>% count(outcome, setting, architecture, target_label)
stopifnot(nrow(cell_counts) == 112L, all(cell_counts$n == 2L))
stopifnot(setequal(unique(plot_base$method_group), c("Standard", "EACC")))

figure_long_internal <- plot_base %>%
  select(
    outcome, setting, synthesis, information_architecture, target_population_shift,
    method = method_group, absolute_bias, rmse, coverage
  ) %>%
  pivot_longer(
    cols = c(absolute_bias, rmse, coverage),
    names_to = "metric", values_to = "value"
  ) %>%
  arrange(metric, outcome, synthesis, information_architecture, target_population_shift, method)
figure_long <- figure_long_internal %>% select(-setting)
stopifnot(nrow(figure_long) == 672L)
readr::write_csv(
  figure_long,
  file.path(SUMMARY_DIR, paste0(RUN_TAG, "_figure_data_long.csv"))
)

ratio_data <- plot_base %>%
  select(
    scenario_id, outcome, synthesis, setting, information_architecture,
    target_population_shift, method_group, rmse
  ) %>%
  pivot_wider(names_from = method_group, values_from = rmse) %>%
  transmute(
    outcome, synthesis, setting, information_architecture,
    target_population_shift,
    rmse_ratio = EACC / Standard
  ) %>%
  arrange(outcome, synthesis, information_architecture, target_population_shift)
stopifnot(nrow(ratio_data) == 112L, all(is.finite(ratio_data$rmse_ratio)))
readr::write_csv(
  ratio_data,
  file.path(SUMMARY_DIR, paste0(RUN_TAG, "_FigureS1_rmse_ratio_data.csv"))
)

FONT_FAMILY <- "Times New Roman"
# Restrained two-colour palette inspired by Nature Portfolio figure styling:
# deep slate blue for the reference method and warm vermilion for EACC.
series_colors <- c(Standard = "#3C5488", EACC = "#E64B35")

panel_specs <- tibble::tribble(
  ~panel, ~outcome,      ~setting,    ~title,
  "A",    "binary",     "pairwise", "A  Binary outcome\nPairwise meta-analysis",
  "B",    "binary",     "nma",      "B  Binary outcome\nNetwork meta-analysis",
  "C",    "continuous", "pairwise", "C  Continuous outcome\nPairwise meta-analysis",
  "D",    "continuous", "nma",      "D  Continuous outcome\nNetwork meta-analysis"
)

method_legend_labels <- function(setting) {
  if (setting == "pairwise") {
    c(Standard = "Meta-analysis", EACC = "EACC-adjusted meta-analysis")
  } else {
    c(Standard = "Network meta-analysis", EACC = "EACC-adjusted network meta-analysis")
  }
}

manuscript_theme <- function() {
  theme_minimal(base_size = 10.5, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY, color = "#202020"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, lineheight = 1.05),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#E3E0DD", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 9.5, color = "#303030"),
      axis.text.y = element_text(size = 9.0, hjust = 1, color = "#303030"),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      legend.direction = "vertical",
      legend.text = element_text(size = 9.3),
      legend.key.width = grid::unit(0.78, "cm"),
      legend.key.height = grid::unit(0.32, "cm"),
      plot.margin = margin(7, 7, 6, 4)
    )
}

scenario_layout <- tidyr::expand_grid(
  information_architecture = factor(architecture_order, levels = architecture_order),
  target_population_shift = factor(target_order, levels = target_order)
) %>%
  mutate(
    architecture_index = match(information_architecture, architecture_order),
    shift_index = match(target_population_shift, target_order),
    y_pos = 36 - ((architecture_index - 1L) * 5L + shift_index + 1L)
  )

architecture_headers <- tibble::tibble(
  information_architecture = architecture_order,
  architecture_index = seq_along(architecture_order),
  y_pos = 36 - ((architecture_index - 1L) * 5L + 1L)
)

scenario_breaks <- scenario_layout %>%
  arrange(y_pos) %>%
  transmute(y_pos, label = paste0("   ", as.character(target_population_shift)))
scenario_label_lookup <- stats::setNames(scenario_breaks$label, as.character(scenario_breaks$y_pos))

add_scenario_structure <- function(plot, x_limits) {
  header_x <- x_limits[1] + 0.015 * diff(x_limits)
  plot +
    geom_rect(
      data = architecture_headers,
      aes(ymin = y_pos - 0.42, ymax = y_pos + 0.42),
      xmin = -Inf, xmax = Inf,
      inherit.aes = FALSE, fill = "#F1EFED", color = "#CEC8C3", linewidth = 0.3
    ) +
    geom_text(
      data = architecture_headers,
      aes(x = header_x, y = y_pos, label = information_architecture),
      inherit.aes = FALSE, hjust = 0, fontface = "bold", size = 3.25,
      family = FONT_FAMILY, color = "#202020"
    )
}

metric_number <- function(metric_name, value) {
  if (metric_name == "coverage") sprintf("%.3f", value) else sprintf("%.3f", value)
}

make_metric_panel <- function(metric_name, outcome_name, setting_name, panel_title, x_label, x_limits) {
  dat <- figure_long_internal %>%
    filter(metric == metric_name, outcome == outcome_name, setting == setting_name) %>%
    left_join(
      scenario_layout %>% select(information_architecture, target_population_shift, y_pos),
      by = c("information_architecture", "target_population_shift")
    ) %>%
    mutate(
      method_y = y_pos + if_else(method == "Standard", 0.16, -0.16),
      value_label = metric_number(metric_name, value)
    )
  labels <- method_legend_labels(setting_name)
  label_offset <- 0.012 * diff(x_limits)
  p <- ggplot()
  if (metric_name == "coverage") {
    p <- p + geom_vline(xintercept = 0.95, color = "grey40", linetype = "dashed", linewidth = 0.45)
  }
  p <- add_scenario_structure(p, x_limits)
  p +
    geom_segment(
      data = dat,
      aes(x = 0, xend = value, y = method_y, yend = method_y, color = method),
      inherit.aes = FALSE, linewidth = 4.3, lineend = "butt"
    ) +
    geom_text(
      data = dat,
      aes(x = value + label_offset, y = method_y, label = value_label),
      inherit.aes = FALSE, hjust = 0, size = 2.85, family = FONT_FAMILY, color = "#202020"
    ) +
    scale_color_manual(values = series_colors, breaks = c("Standard", "EACC"), labels = labels, name = NULL) +
    scale_x_continuous(
      limits = x_limits,
      labels = if (metric_name == "coverage") scales::label_percent(accuracy = 1) else scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0.5, 35.5), breaks = scenario_breaks$y_pos,
      labels = function(x) unname(scenario_label_lookup[as.character(x)]),
      expand = expansion(mult = 0)
    ) +
    labs(title = panel_title, x = x_label, y = NULL) +
    guides(
      color = guide_legend(ncol = 1, override.aes = list(linewidth = 4.3))
    ) +
    manuscript_theme()
}

make_metric_figure <- function(metric_name, figure_title, x_label, x_limits) {
  panels <- lapply(seq_len(nrow(panel_specs)), function(i) {
    make_metric_panel(
      metric_name = metric_name,
      outcome_name = panel_specs$outcome[i],
      setting_name = panel_specs$setting[i],
      panel_title = panel_specs$title[i],
      x_label = x_label,
      x_limits = x_limits
    )
  })
  wrap_plots(panels, ncol = 4, guides = "keep") +
    plot_annotation(
      title = figure_title,
      theme = theme(
        text = element_text(family = FONT_FAMILY, color = "#202020"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10))
      )
    )
}

make_ratio_panel <- function(outcome_name, setting_name, panel_title, x_limits) {
  dat <- ratio_data %>%
    filter(outcome == outcome_name, setting == setting_name) %>%
    left_join(
      scenario_layout %>% select(information_architecture, target_population_shift, y_pos),
      by = c("information_architecture", "target_population_shift")
    ) %>%
    mutate(value_label = sprintf("%.2f", rmse_ratio))
  label_offset <- 0.012 * diff(x_limits)
  add_scenario_structure(
    ggplot() + geom_vline(xintercept = 1, color = "grey40", linetype = "dashed", linewidth = 0.5),
    x_limits
  ) +
    geom_segment(
      data = dat,
      aes(x = 0, xend = rmse_ratio, y = y_pos, yend = y_pos),
      inherit.aes = FALSE, color = series_colors[["EACC"]], linewidth = 4.3, lineend = "butt"
    ) +
    geom_text(
      data = dat, aes(x = rmse_ratio + label_offset, y = y_pos, label = value_label),
      inherit.aes = FALSE, hjust = 0, size = 2.85, family = FONT_FAMILY, color = "#202020"
    ) +
    scale_x_continuous(
      limits = x_limits, labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0.5, 35.5), breaks = scenario_breaks$y_pos,
      labels = function(x) unname(scenario_label_lookup[as.character(x)]),
      expand = expansion(mult = 0)
    ) +
    labs(
      title = panel_title,
      x = "RMSE ratio (EACC / Standard)",
      y = NULL
    ) +
    manuscript_theme() +
    theme(legend.position = "none")
}

make_ratio_figure <- function(x_limits) {
  panels <- lapply(seq_len(nrow(panel_specs)), function(i) {
    make_ratio_panel(
      outcome_name = panel_specs$outcome[i],
      setting_name = panel_specs$setting[i],
      panel_title = panel_specs$title[i],
      x_limits = x_limits
    )
  })
  wrap_plots(panels, ncol = 4) +
    plot_annotation(
      title = "Relative RMSE of EACC versus standard synthesis",
      subtitle = "Values below 1 favor EACC; values above 1 favor standard synthesis",
      theme = theme(
        text = element_text(family = FONT_FAMILY, color = "#202020"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10))
      )
    )
}

metric_limits <- list(
  absolute_bias = c(0, max(figure_long$value[figure_long$metric == "absolute_bias"], na.rm = TRUE) * 1.18),
  rmse = c(0, max(figure_long$value[figure_long$metric == "rmse"], na.rm = TRUE) * 1.18),
  coverage = c(0, 1.12),
  ratio = c(0, max(1.08, max(ratio_data$rmse_ratio, na.rm = TRUE) * 1.12))
)

figures <- list(
  Figure1_absolute_bias = make_metric_figure(
    "absolute_bias", "Absolute bias by target population shift and information architecture",
    "Absolute bias", metric_limits$absolute_bias
  ),
  Figure2_RMSE = make_metric_figure(
    "rmse", "RMSE by target population shift and information architecture",
    "RMSE", metric_limits$rmse
  ),
  Figure3_Coverage = make_metric_figure(
    "coverage", "Coverage by target population shift and information architecture",
    "95% CI coverage", metric_limits$coverage
  ),
  FigureS1_RMSE_ratio = make_ratio_figure(metric_limits$ratio)
)

save_figure <- function(plot, stem) {
  png_path <- file.path(FIGURE_DIR, paste0(stem, ".png"))
  pdf_path <- file.path(FIGURE_DIR, paste0(stem, ".pdf"))
  # Freeze the complete patchwork layout once so PNG and PDF have identical
  # panels, architecture headings, scenario rows, axes, and legends.
  plot_grob <- patchwork::patchworkGrob(plot)
  ggsave(pdf_path, plot = plot_grob, width = 24, height = 16, units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(
    png_path, plot = plot_grob, width = 24, height = 16, units = "in", dpi = 300,
    device = grDevices::png, type = "cairo", bg = "white", limitsize = FALSE
  )
  c(png = png_path, pdf = pdf_path)
}

paths <- Map(save_figure, figures, names(figures))

qc <- tibble::tibble(
  figure = c("Figure 1: Absolute bias", "Figure 2: RMSE", "Figure 3: Coverage", "Figure S1: RMSE ratio"),
  panels = 4L,
  panel_layout = "one row x four columns",
  architecture_sections_per_panel = 7L,
  scenario_rows_per_architecture = 4L,
  bars_per_scenario = c(2L, 2L, 2L, 1L),
  source_rows = c(224L, 224L, 224L, 112L),
  expected_combinations_complete = TRUE
)
readr::write_csv(qc, file.path(FIGURE_DIR, paste0(RUN_TAG, "_figure_QC.csv")))

message("Figure data rows: ", nrow(figure_long))
message("RMSE ratio rows: ", nrow(ratio_data))
message("Generated figures:")
print(unlist(paths))
print(qc)
