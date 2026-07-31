library(dplyr)
library(ggplot2)
source("./R/utils.R")
source("./R/plot_funs.R")

tot_animals <- 100

save_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/D_all_results"
fig_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/imgs/"

fig_colors <- c("#1B5E20", "#00A8C6", "#FBC02D", "#E65100", "#8E44AD", "#4B6FAD", "#D81B60")

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

  # file_ii <- paste0(sim_dir, design_names[ii], "_", ncams, "_cam_10min.RData")
  file_ii <- paste0(sim_dir, design_names[ii], "_", ncams, "_cam.RData")
  if (file.exists(file_ii)) {
    all_results <- loadRData(file_ii)

    D_all <- D_all |>
      dplyr::bind_rows(
        all_results[[5]] |>
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
print(skipped_files)
D_all |>
  dplyr::count(Model, Covariate, SampDesign)

# # save(D_all, file = paste0(
# #   save_dir,
# #   "D_all.RData")
# # )

#--------------------------------------------------

load(paste0(
  save_dir,
  "D_all.RData")
)

D_all <- D_all |>
  dplyr::mutate(
    Model = dplyr::case_when(
      Covariate == "Covariate" & Model == "REST" ~ "REST (Cov)",
      Covariate == "Non-Covariate" & Model == "REST" ~ "REST (Non-Cov)",
      .default = Model
    )
  )

D_all <- D_all %>%
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
  )

D_all$SampDesign <- factor(
  D_all$SampDesign,
  levels = c("Random", "80% Low", "100% Low", "80% Moderate", "100% Moderate",
             "80% High", "100% High")
)
IS_random <- D_all |>
  dplyr::filter(Model == "IS" & SampDesign == "Random") |>
  dplyr::summarise(
    Mean_MAE = mean(abs(tot_animals - Est), na.rm = T),
    Mean_Est = mean(Est, na.rm = T),
    Mean_var = mean(SD, na.rm = T)
  )

# Random
D_all %>%
  dplyr::filter(SampDesign == "Random") %>%
  ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Model",
                y = "Posterior Mean") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  # patchwork::plot_annotation(tag_levels = 'a') +
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
#          file_names[1],
#          "_cam.pdf"),
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

D_all %>%
  dplyr::filter(SampDesign == "Random") %>%
  ggplot2::ggplot(ggplot2::aes(x = Model, y = SD, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::labs(x = "Model",
                y = "Posterior Variance") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
  # patchwork::plot_annotation(tag_levels = 'a') +
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
#          file_names[1],
#          "_cam_SD.pdf"),
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
    # BiasLevel = paste(BiasLevel, "Density"),
    BiasLevel = factor(BiasLevel, levels = c("Low", "Moderate", "High")),
    Percentage = factor(Percentage, levels = c("80%", "100%")),
    Est = ifelse(
      Model == "REST (Cov)" & SampDesign %in% c("100% High", "100% Moderate", "100% Low"),
      NA,
      Est
    )
  )

annotation_df <- data.frame(
  x = Inf,
  y = Inf,
  label = "a",
  BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
  # Percentage = factor("100%", levels =  c("80%", "100%"))
)

# All other sample designs
D_all_separated %>%
  dplyr::filter(SampDesign != "Random") |>
  ggplot2::ggplot(ggplot2::aes(x = Percentage, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Sampling Bias",
                y = "Posterior Mean") +
  ggplot2::facet_grid(~ BiasLevel) + #, switch = "x") +
  ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
  ggplot2::geom_text(
    data = annotation_df,
    mapping = ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 2.5,   # Horizontal adjustment (same as before)
    vjust = 1.5,   # Vertical adjustment (same as before)
    size = 5
  ) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 legend.title=element_text(size=10),
                 legend.text=element_text(size=9),
                 legend.position = c(0.092, 0.793),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 panel.spacing = unit(0,'lines'),
                 axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth = 0.5),
                 legend.background = ggplot2::element_rect(color = "black"),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))

# ggplot2::ggsave(
#   paste0(fig_dir,
#          "bias_cam.pdf"),
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
annotation_df_b <- data.frame(
  x = Inf,
  y = Inf,
  label = "b",
  BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
  # Percentage = factor("100%", levels =  c("80%", "100%"))
)
annotation_IS <- data.frame(
  x = Inf,
  y = Inf,
  label = "IS Variance",
  BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
  # Percentage = factor("100%", levels =  c("80%", "100%"))
)


D_all_separated %>%
  dplyr::filter(SampDesign != "Random" & !(Model %in% c("IS", "REST (Non-Cov)"))) |>
  ggplot2::ggplot(ggplot2::aes(x = Percentage, y = SD, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=IS_random$Mean_var, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Sampling Bias",
                y = "Posterior Variance") +
  ggplot2::facet_grid(~ BiasLevel) + #, switch = "x") +
  # ggplot2::guides(
  #   fill = ggplot2::guide_legend(
  #     title = "Model",
  #     override.aes = list(
  #       linetype = c("blank", "dashed"),
  #       fill = c("#00A8C6", NA)
  #     )
  #   )
  # ) +
  ggplot2::geom_text(
    data = annotation_df_b,
    mapping = ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 2.5,   # Horizontal adjustment (same as before)
    vjust = 1.5,   # Vertical adjustment (same as before)
    size = 5
  ) +
  ggplot2::geom_text(
    data = annotation_IS,
    mapping = ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    x = 2,
    y = IS_random$Mean_var + 0.8,
    size = 3
  ) +
  # annotate("text",
  #          x = 6,
  #          y = IS_random$Mean_var + 0.7,
  #          size = 3,
  #          label = "IS Variance",
  #          color = "black") +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 legend.title=element_text(size=10),
                 legend.text=element_text(size=9),
                 legend.position = c(0.092, 0.793),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 panel.spacing = unit(0,'lines'),
                 axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
                 panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth = 0.5),
                 legend.background = ggplot2::element_rect(color = "black"),
                 legend.spacing.y = ggplot2::unit(0, "mm"),
                 legend.box.background = ggplot2::element_rect(colour = "black"))


# ggplot2::ggsave(
#   paste0(fig_dir,
#          "PATH_bias_cam_SD.pdf"),
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

#--------------------------------------------------
# # Try all in one plot?
# D_all %>%
#   dplyr::filter(Est < 250) %>%
#   dplyr::mutate(
#     SampDesign = dplyr::case_when(
#       SampDesign == "slow_cam" ~ "80% High",
#       SampDesign == "med_cam" ~ "80% Moderate",
#       SampDesign == "fast_cam" ~ "80% Low",
#       SampDesign == "all_slow_cam" ~ "100% High",
#       SampDesign == "all_med_cam" ~ "100% Moderate",
#       SampDesign == "all_fast_cam" ~ "100% Low",
#       .default = "Random"
#
#     )
#   ) %>%
#   ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = SampDesign)) +
#   ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
#   ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
#   ggplot2::labs(x = "Sampling Bias",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors[1:7]) +
#   ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#                  # legend.position = "none",
#                  # legend.background = ggplot2::element_blank(),
#                  # legend.spacing.y = ggplot2::unit(0, "mm"),
#                  # legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_all %>%
#   dplyr::filter(Est < 250) %>%
#   dplyr::mutate(
#     SampDesign = dplyr::case_when(
#       SampDesign == "slow_cam" ~ "80% High",
#       SampDesign == "med_cam" ~ "80% Moderate",
#       SampDesign == "fast_cam" ~ "80% Low",
#       SampDesign == "all_slow_cam" ~ "100% High",
#       SampDesign == "all_med_cam" ~ "100% Moderate",
#       SampDesign == "all_fast_cam" ~ "100% Low",
#       .default = "Random"
#
#     )
#   ) %>%
#   ggplot2::ggplot(ggplot2::aes(x = Model, y = SD / Est, fill = SampDesign)) +
#   ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
#   ggplot2::labs(x = "Sampling Bias",
#                 y = "Posterior CVs") +
#   ggplot2::scale_fill_manual(values= fig_colors[1:7]) +
#   ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#                  # legend.position = "none",
#                  # legend.background = ggplot2::element_blank(),
#                  # legend.spacing.y = ggplot2::unit(0, "mm"),
#                  # legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
#
