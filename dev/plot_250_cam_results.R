library(dplyr)
library(ggplot2)
source("./R/utils.R")
source("./R/plot_funs.R")

tot_animals <- 100

save_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/D_all_results"
fig_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/imgs/"

fig_colors <- c("#1B5E20", "#00A8C6", "#E65100", "#FBC02D", "#8E44AD", "#4B6FAD", "#D81B60")

################################################################################
# Clean raw simulation data
sim_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/250_cam_results/"
design_names <- c(
  "Random", "Slow_80_bias", "Medium_80_bias", "Fast_80_bias", "Slow_bias", "Medium_bias", "Fast_bias"
)
ncams <- 250
file_names <- c(
  "random",
  "slow",
  "med",
  "fast",
  "all_slow",
  "all_med",
  "all_fast"
)

# Collect density estimates for from all results
# (This makes it easier to analyze)
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

skipped_files <- character(0)
D_all <- tibble::tibble()
for (ii in 1:length(file_names)) {
  print(ii)

  file_ii <- paste0(sim_dir, design_names[ii], "_", ncams, "_cam_10min.RData")
  # file_ii <- paste0(sim_dir, design_names[ii], "_", ncams, "_cam.RData")
  if (file.exists(file_ii)) {
    all_results <- loadRData(file_ii)

    D_all <- D_all |>
      dplyr::bind_rows(
        all_results[[4]] |>
          dplyr::bind_rows() |>
          dplyr::mutate(
            SampDesign = paste0(file_names[ii], "_cam")
          )
      )

  } else {
    skipped_files <- c(skipped_files, file_ii)
  }

  file_ii_REST <- paste0(sim_dir, design_names[ii], "_", ncams, "_cam_REST.RData")
  if (file.exists(file_ii_REST)) {
    all_results_REST <- loadRData(file_ii_REST)

    D_all <- D_all |>
      dplyr::bind_rows(
        all_results_REST[[4]] %>%
          dplyr::bind_rows() |>
          dplyr::mutate(
            SampDesign = paste0(file_names[ii], "_cam")
          )
      )

  } else {
    skipped_files <- c(skipped_files, file_ii_REST)
  }
}
# print(skipped_files)
# D_all |>
#   dplyr::count(Model, Covariate, SampDesign)
#
# save(D_all, file = paste0(
#   save_dir,
#   "D_all.RData")
# )

#--------------------------------------------------

load(paste0(
  save_dir,
  "D_all.RData")
)

D_all <- D_all |>
  dplyr::filter(Covariate != "Covariate") %>%
  # dplyr::filter(Est < 250) %>%
  dplyr::mutate(
    SampDesign = dplyr::case_when(
      SampDesign == "slow_cam" ~ "80% High",
      SampDesign == "med_cam" ~ "80% Moderate",
      SampDesign == "fast_cam" ~ "80% Low",
      SampDesign == "all_slow_cam" ~ "100% High",
      SampDesign == "all_med_cam" ~ "100% Moderate",
      SampDesign == "all_fast_cam" ~ "100% Low",
      .default = "Random"

    )
  ) %>%
  dplyr::mutate(
    MAE = abs(tot_animals - Est),
    Var = SD ^ 2,
    log_mu = log(Est ^2 / sqrt(Var + Est^2)),
    log_sigma = sqrt(log(1 + (Var / Est^2))),
    LCL_95 = ifelse(
      Model == "IS",
      Est / exp(1.96 * sqrt(log(1 + (SD / Est) ^ 2))),
      qlnorm(0.025, meanlog = log_mu, sdlog = log_sigma)
    ),
    UCL_95 = ifelse(
      Model == "IS",
      Est * exp(1.96 * sqrt(log(1 + (SD / Est) ^ 2))),
      qlnorm(0.975, meanlog = log_mu, sdlog = log_sigma)
    ),
    CrI = UCL_95 - LCL_95,
    coverage_95 = ifelse(tot_animals >= LCL_95 & tot_animals <= UCL_95, 1, 0)
  )

D_all$SampDesign <- factor(
  D_all$SampDesign,
  levels = c("Random", "80% Low", "100% Low", "80% Moderate", "100% Moderate",
             "80% High", "100% High")
)

IS_random <- D_all |>
  dplyr::filter(Model == "IS" & SampDesign == "Random") %>%
  dplyr::summarise(
    Mean_MAE = mean(MAE, na.rm = T),
    Mean_Est = mean(Est, na.rm = T),
    SD_Est = sd(Est, na.rm = T),
    Mean_SD = mean(SD, na.rm = T),
    Mean_CrI = mean(CrI, na.rm = T),
    coverage_prob = mean(coverage_95, na.rm = T),
    Relative_Bias = mean((Est - tot_animals) / tot_animals, na.rm = T),
    RMSE = sqrt(mean((Est - tot_animals) ^ 2, na.rm = T)),
    .groups = 'drop'
  )

REST_random <- D_all |>
  dplyr::filter(Model == "REST" & Covariate == "Non-Covariate" & SampDesign == "Random") %>%
  dplyr::summarise(
    Mean_MAE = mean(MAE, na.rm = T),
    Mean_Est = mean(Est, na.rm = T),
    SD_Est = sd(Est, na.rm = T),
    Mean_SD = mean(SD, na.rm = T),
    Mean_CrI = mean(CrI, na.rm = T),
    coverage_prob = mean(coverage_95, na.rm = T),
    Relative_Bias = mean((Est - tot_animals) / tot_animals, na.rm = T),
    RMSE = sqrt(mean((Est - tot_animals) ^ 2, na.rm = T)),
    .groups = 'drop'
  )

# Random Mean estimates
D_all %>%
  dplyr::filter(SampDesign == "Random") %>%
  ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Model",
                y = "Posterior Mean") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::annotate("text", x = -Inf, y = Inf,
                    label = "a", hjust = -1, vjust = 1.5,
                    size = 5) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 legend.title = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 legend.position = "none",
                 legend.background = ggplot2::element_blank(),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

# ggplot2::ggsave(
#   paste0(fig_dir,
#          "random_cam.pdf"),
#   plot = ggplot2::last_plot(),
#   # path = file_path,
#   # scale = 1,
#   width = 5,
#   height = 3,
#   # units = c("in", "cm", "mm", "px"),
#   dpi = 600,
#   limitsize = TRUE,
#   bg = NULL
# )

# Random credible interval width
D_all %>%
  dplyr::filter(SampDesign == "Random") %>%
  ggplot2::ggplot(ggplot2::aes(x = Model, y = CrI, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::labs(x = "Model",
                y = "Credible Interval Width") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  ggplot2::annotate("text", x = -Inf, y = Inf,
                    label = "b", hjust = -1, vjust = 1.5,
                    size = 5) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 legend.title = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                 legend.position = "none",
                 legend.background = ggplot2::element_blank(),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

# ggplot2::ggsave(
#   paste0(fig_dir,
#          "_cam_cri_width.pdf"),
#   plot = ggplot2::last_plot(),
#   # path = file_path,
#   # scale = 1,
#   width = 5,
#   height = 3,
#   # units = c("in", "cm", "mm", "px"),
#   dpi = 600,
#   limitsize = TRUE,
#   bg = NULL
# )


################################################################################
# Plot sampling bias results
D_all_separated <- D_all %>%
  dplyr::filter(SampDesign != "Random") %>%
  tidyr::separate(SampDesign, into = c("Percentage", "BiasLevel"), sep = " ", remove = FALSE) %>%
  dplyr::mutate(
    BiasLevel = factor(BiasLevel, levels = c("Low", "Moderate", "High")),
    Percentage = factor(Percentage, levels = c("80%", "100%"))
  )

annotation_df <- data.frame(
  x = Inf,
  y = Inf,
  label = "a",
  BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
)

# All other sample designs
D_all_separated %>%
  dplyr::filter(SampDesign != "Random") |>
  ggplot2::ggplot(ggplot2::aes(x = Percentage, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", linewidth = 0.7) +
  ggplot2::labs(x = "Sampling Bias",
                y = "Posterior Mean") +
  ggplot2::facet_grid(~ BiasLevel) +
  ggplot2::scale_fill_manual(values = fig_colors[1:5]) +
  ggplot2::geom_text(
    data = annotation_df,
    mapping = ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 2.5,
    vjust = 1.5,
    size = 3
  ) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7),
    legend.key.size = ggplot2::unit(0.4, "cm"),
    legend.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"),
    legend.position = c(0.075, 0.815),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.spacing = ggplot2::unit(0, 'lines'),
    axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
    panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth = 0.5),
    legend.background = ggplot2::element_rect(color = "black", linewidth = 0.5),
    legend.spacing.y = ggplot2::unit(0, "mm"),
    legend.box.background = ggplot2::element_blank()
  )

# ggplot2::ggsave(
#   paste0(fig_dir, "bias_cam.pdf"),
#   plot = ggplot2::last_plot(),
#   width = 5,
#   height = 3,
#   dpi = 600,
#   limitsize = TRUE,
#   bg = "white"
# )

################################################################################
# Plot CrI width for PATH
annotation_df_b <- data.frame(
  x = Inf,
  y = Inf,
  label = "b",
  BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
)
annotation_IS <- data.frame(
  x = 2.05,
  y_CrI = IS_random$Mean_CrI + 0.4,
  label = "IS Random",
  BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
)

annotation_REST <- data.frame(
  x = 2,
  y_CrI = REST_random$Mean_CrI + 0.4,
  label = "REST Random",
  BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
)


D_all_separated %>%
  dplyr::filter(SampDesign != "Random" & !(Model %in% c("IS", "REST"))) |>
  ggplot2::ggplot(ggplot2::aes(x = Percentage, y = CrI, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = NA) +
  ggplot2::scale_fill_manual(values = fig_colors[2]) +
  ggplot2::geom_hline(yintercept=IS_random$Mean_CrI, linetype="dashed", linewidth = 0.7, color = fig_colors[1]) +
  ggplot2::geom_hline(yintercept=REST_random$Mean_CrI, linetype="dashed", linewidth = 0.7, color = fig_colors[3]) +
  ggplot2::labs(x = "Sampling Bias",
                y = "Credible Interval Width") +
  ggplot2::facet_grid(~ BiasLevel) +
  ggplot2::expand_limits(y = max(IS_random$Mean_CrI, REST_random$Mean_CrI) + 1.5) +
  ggplot2::geom_text(
    data = annotation_df_b,
    mapping = ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 2.5,
    vjust = 1.5,
    size = 5
  ) +
  ggplot2::geom_text(
    data = annotation_IS,
    mapping = ggplot2::aes(x = x, y = y_CrI, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  ggplot2::geom_text(
    data = annotation_REST,
    mapping = ggplot2::aes(x = x, y = y_CrI, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7),
    legend.key.size = ggplot2::unit(0.4, "cm"),
    legend.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"),
    legend.position = c(0.07, 0.886),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.spacing = ggplot2::unit(0, 'lines'),
    axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
    panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth = 0.5),
    legend.background = ggplot2::element_rect(color = "black", linewidth = 0.5),
    legend.spacing.y = ggplot2::unit(0, "mm"),
    legend.box.background = ggplot2::element_blank()
  )


# ggplot2::ggsave(
#   paste0(fig_dir,
#          "PATH_bias_cam_cri_width.pdf"),
#   plot = ggplot2::last_plot(),
#   width = 5,
#   height = 3,
#   dpi = 600,
#   limitsize = TRUE,
#   bg = "white"
# )

