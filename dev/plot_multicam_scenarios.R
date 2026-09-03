library(dplyr)
library(ggplot2)
tot_animals <- 100

save_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/D_all_results/"
img_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/imgs"

fig_colors <- c("#1B5E20", "#00A8C6", "#E65100", "#FBC02D", "#8E44AD", "#4B6FAD", "#D81B60")

################################################################################
#  Download and clean raw sim data
sim_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/multi_cam_results/"
design_names <- c(
  "Random", "Slow_80_bias", "Medium_80_bias", "Fast_80_bias", "Slow_bias", "Medium_bias", "Fast_bias"
)

ncams <- c(25, 50, 75, 100, 125)

# Collect density estimates for from all results
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

skipped_files <- character(0)
D_all_multi_cam <- c()
for (dn in design_names) {
  for (nc in ncams) {
    print(dn)

    file_ii <- paste0(sim_dir, dn, "_", nc, "_cam_10min.RData")
    if (file.exists(file_ii)) {
      all_results <- loadRData(file_ii)

      D_all_multi_cam <- D_all_multi_cam |>
        dplyr::bind_rows(
          all_results[[5]]
        )

    } else {
      skipped_files <- c(skipped_files, file_ii)
    }

    file_ii_REST <- paste0(sim_dir, dn, "_", nc, "_cam_REST.RData")
    if (file.exists(file_ii_REST)) {
      all_results_REST <- loadRData(file_ii_REST)

      D_all_multi_cam <- D_all_multi_cam |>
        dplyr::bind_rows(
          all_results_REST[[4]] %>%
            dplyr::bind_rows() |>
            dplyr::mutate(
              cam_design = dn,
              cams = nc
            )
        )
    } else {
      skipped_files <- c(skipped_files, file_ii_REST)
    }

  }
}

print(skipped_files)
D_all_multi_cam |>
  dplyr::count(Model, Covariate, cam_design, cams)
#
# save(D_all_multi_cam, file = paste0(
#   save_dir,
#   "D_all_multi_cam.RData")
# )

################################################################################
# Upload base data for IS random at 250 cams
load(paste0(
  save_dir,
  "D_all.RData")
)

IS_random <- D_all |>
  dplyr::filter(Model == "IS" & SampDesign == "random_cam") %>%
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
  ) |> # & Est < 250) |>
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
  ) %>%
  dplyr::mutate(
    cams = 250,
    Model = "IS Random"
  )

REST_random <- D_all |>
  dplyr::filter(Model == "REST" & Covariate == "Non-Covariate" & SampDesign == "random_cam") %>%
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
  ) |> # & Est < 250) |>
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
  ) |>
  dplyr::mutate(
    cams = 250,
    Model = "REST Random"
  )

# Load multicam data
load(paste0(
  save_dir,
  "D_all_multi_cam.RData")
)

# Format data
D_all_multi_cam <- D_all_multi_cam |>
  dplyr::filter(Covariate != "Covariate") %>%
  dplyr::mutate(
    Design = dplyr::case_when(
      cam_design %in% c("Slow_80_bias") ~ "80% High",
      cam_design %in% c("Medium_80_bias") ~ "80% Moderate",
      cam_design %in% c("Fast_80_bias") ~ "80% Low",
      cam_design %in% c("Slow_bias") ~ "100% High",
      cam_design %in% c("Medium_bias") ~ "100% Moderate",
      cam_design %in% c("Fast_bias") ~ "100% Low",
      .default = cam_design
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

# # Outliers defined as 1.5 times the interquartile range
# # So 250 is about the mean upper bound over all results
# outlier <- D_all_multi_cam  %>%
#   # dplyr::group_by(Model, Covariate, Design) %>%
#   dplyr::summarise(
#     q1 = quantile(Est, 0.25, na.rm = T),
#     q3 = quantile(Est, 0.75, na.rm = T),
#     lower_bound = q1 - (1.5 * (q3 - q1)),
#     upper_bound = q3 + (1.5 * (q3 - q1))
#   )

# Factor variables
D_all_multi_cam$Model <- factor(
  D_all_multi_cam$Model,
  levels = c("IS", "PATH", "REST")
)
D_all_multi_cam$Design <- factor(
  D_all_multi_cam$Design,
  levels = c("Random", "80% Low", "100% Low", "80% Moderate", "100% Moderate", "80% High", "100% High")
)

################################################################################
# Summarize results
D_cam_means <- D_all_multi_cam |>
  dplyr::group_by(cams, Model, Design) %>%
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
  ) %>%
  dplyr::mutate(
    Speed = dplyr::case_when(
      Design %in% c("80% High", "100% High") ~ "High",
      Design %in% c("80% Moderate", "100% Moderate") ~ "Moderate",
      Design %in% c("80% Low", "100% Low") ~ "Low",
      Design %in% c("Random") ~ "Random"
    ),
    Bias = dplyr::case_when(
      Design %in% c("80% High", "80% Moderate", "80% Low") ~ "80%",
      Design %in% c("100% High", "100% Moderate", "100% Low") ~ "100%",
      Design %in% c("Random") ~ "Random"
    )
  )

# Plot all in one plot
D_cam_means_fctr <- D_cam_means %>%
  dplyr::rename(Model_char = Model) |>
  dplyr::mutate(Model = paste(Model_char, Design)) |>
  dplyr::filter(
    Model %in% c("IS Random", "REST Random", "PATH Random",
                 "PATH 80% High", "PATH 100% High",
                 "PATH 80% Moderate", "PATH 100% Moderate",
                 "PATH 80% Low", "PATH 100% Low")
  )

D_cam_means_fctr$Model <- factor(
  D_cam_means_fctr$Model,
  levels = c("IS Random", "REST Random", "PATH Random",
             "PATH 80% High", "PATH 100% High",
             "PATH 80% Moderate", "PATH 100% Moderate",
             "PATH 80% Low", "PATH 100% Low"
             )
)

# Colors: Baselines get distinct colors; Biased PATH models grouped by bias level
model_colors <- c(
  "IS Random" = "#1B5E20",
  "REST Random" = "#E65100",
  "Random" = "#00A8C6",
  "80% Low" = "#FBC02D",
  "100% Low" = "#FBC02D",
  "80% Moderate" = "#8E44AD",
  "100% Moderate" = "#8E44AD",
  "80% High" = "#4B6FAD",
  "100% High" = "#4B6FAD"
)

# Linetypes: Baselines = Solid; 80% = Dashed; 100% = Dotted
model_linetypes <- c(
  "IS Random" = "solid",
  "REST Random" = "solid",
  "Random" = "solid",
  "80% Low" = "dashed",
  "100% Low" = "dotted",
  "80% Moderate" = "dashed",
  "100% Moderate" = "dashed",
  "80% High" = "dashed",
  "100% High" = "dotdash"
)

# Linewidth: Make the baselines slightly thicker so they stand out
model_linewidths <- c(
  "IS Random" = 1.2,
  "REST Random" = 1.2,
  "PATH Random" = 1.2,
  "PATH 80% Low" = 0.8,
  "PATH 100% Low" = 0.8,
  "PATH 80% Moderate" = 0.8,
  "PATH 100% Moderate" = 0.8,
  "PATH 80% High" = 0.8,
  "PATH 100% High" = 0.8
)


# D_cam_means_fctr |>
# dplyr::filter(
#   Bias %in% c("Random", "100%")
# ) |>
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_SD,
#                                color = Model,
#                                linetype = Model,
#                                linewidth = Model)) +
#   ggplot2::geom_line() +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::scale_colour_manual(name = "Model", values = model_colors) +
#   ggplot2::scale_linetype_manual(name = "Model", values = model_linetypes) +
#   ggplot2::scale_linewidth_manual(name = "Model", values = model_linewidths) +
#   # ggplot2::scale_shape_manual(name = "Model", values = c(
#   #   "IS Random" = 16, "REST Random" = 17, "PATH Random" = 15,
#   #   "PATH 80% Low" = 1, "PATH 100% Low" = 1,
#   #   "PATH 80% Moderate" = 1, "PATH 100% Moderate" = 1,
#   #   "PATH 80% High" = 1, "PATH 100% High" = 1
#   # )) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Mean Variance") +
#   ggplot2::scale_x_continuous(breaks = c(10, 25, 50, 75, 100),
#                               labels = c("10", "25", "50", "75", "100")) +
#   ggplot2::annotate("text", x = -Inf, y = Inf,
#                     label = "b", hjust = -1, vjust = 1.5,
#                     size = 4) +
#   ggplot2::theme(
#     text = ggplot2::element_text(size = 10),
#     axis.title = ggplot2::element_text(size = 12),
#     axis.text = ggplot2::element_text(size = 10),
#     legend.position = "right",
#     legend.text = ggplot2::element_text(size = 8),
#     legend.title = ggplot2::element_text(size = 10),
#     legend.key.width = ggplot2::unit(1.5, "cm"),    panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
#     legend.background = ggplot2::element_rect(color = NA)
#   )
#
# ggplot2::ggsave(
#   paste0(img_dir, "/cam_vs_SD.pdf"),
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

# # Plot all MAEs
# D_cam_means_fctr |>
#   dplyr::filter(
#     Bias %in% c("Random", "100%")
#   ) |>
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_MAE, color = Model, linetype = Model)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   # custom_colors +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "solid",
#       "PATH Random" = "solid",
#       "REST (Cov) Random" = "solid",
#       "PATH 80% High" = "dashed",
#       "REST (Cov) 80% High" = "dashed",
#       "PATH 100% High" = "dotted",
#       "PATH 80% Moderate" = "dashed",
#       "REST (Cov) 80% Moderate" = "dashed",
#       "PATH 100% Moderate" = "dotted",
#       "PATH 80% Low" = "dashed",
#       "REST (Cov) 80% Low" = "dashed",
#       "PATH 100% Low" = "dotted"
#     )
#   ) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Mean Errors") +
#   ggplot2::scale_x_continuous(breaks = c(10,25,50, 75, 100),
#                               labels = c("10","25","50", "75", "100")) +
#   ggplot2::annotate("text", x = -Inf, y = Inf,
#                     label = "a", hjust = -1, vjust = 1.5,
#                     size = 5) +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = "none",# c(0.85, 0.72),
#     legend.title = ggplot2::element_text(size=13),
#     legend.text = ggplot2::element_text(size = 13),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_blank(), # ggplot2::element_rect(color = "white"), #
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )

# ggplot2::ggsave(
#   paste0(img_dir, "/cam_vs_MAE.pdf"),
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
# Plot RMSE for PATH with different cameras, IS and REST with 250 cameras
D_cam_means_fctr |>
  dplyr::filter(
    Bias %in% c("Random", "100%"),
    !(Model %in% c("IS Random", "REST Random"))
  ) |>
  dplyr::mutate(Percent_effort = cams / 250) |>
  ggplot2::ggplot(ggplot2::aes(x = Percent_effort, y = RMSE, color = Design, shape = Design)) +
  ggplot2::geom_hline(yintercept = IS_random$RMSE, linetype = "dashed", linewidth = 0.7, color = fig_colors[1]) +
  ggplot2::geom_hline(yintercept = REST_random$RMSE, linetype = "dashed", linewidth = 0.7, color = fig_colors[3]) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    data = IS_random,
    mapping = ggplot2::aes(x = 0.1, y = RMSE + 0.4, label = "IS Random"),
    inherit.aes = FALSE,
    color = fig_colors[1],
    size = 3
  ) +
  ggplot2::geom_text(
    data = REST_random,
    mapping = ggplot2::aes(x = 0.1, y = RMSE + 0.4, label = "REST Random"),
    inherit.aes = FALSE,
    color = fig_colors[3],
    size = 3
  ) +
  ggplot2::scale_colour_manual(name = "PATH Sample Design", values = model_colors) +
  # ggplot2::scale_linetype_manual(name = "PATH Sample Design", values = model_linetypes) +
  ggplot2::scale_x_continuous(
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c("10%", "20%", "30%", "40%", "50%")
  ) +
  ggplot2::scale_shape_manual(name = "PATH Sample Design", values = c(16, 15, 17, 18)) +
  ggplot2::labs(
    x = "Relative Effort",
    y = "Root Mean Squared Error (RMSE)"
  ) +
  ggplot2::annotate("text", x = Inf, y = Inf,
                    label = "a", hjust = 2, vjust = 1.5,
                    size = 5, fontface = "bold") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = c(0.85, 0.75),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    legend.key.width = ggplot2::unit(1.2, "cm"),
    legend.background = ggplot2::element_rect(color = "black", fill = "white", linewidth = 0.5),

    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# ggplot2::ggsave(
#   paste0(img_dir, "/cam_vs_RMSE.pdf"),
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


# Plot empirical SD for PATH with different cameras, IS and REST with 250 cameras
D_cam_means_fctr |>
  dplyr::filter(
    Bias %in% c("Random", "100%"),
    !(Model %in% c("IS Random", "REST Random"))
  ) |>
  dplyr::mutate(Percent_effort = cams / 250) |>
  ggplot2::ggplot(ggplot2::aes(x = Percent_effort, y = SD_Est, color = Design, shape = Design)) +
  ggplot2::geom_hline(yintercept = IS_random$SD_Est, linetype = "dashed", linewidth = 0.7, color = fig_colors[1]) +
  ggplot2::geom_hline(yintercept = REST_random$SD_Est, linetype = "dashed", linewidth = 0.7, color = fig_colors[3]) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    data = IS_random,
    mapping = ggplot2::aes(x = 0.1, y = SD_Est + 0.4, label = "IS Random"),
    inherit.aes = FALSE,
    color = fig_colors[1],
    size = 3
  ) +
  ggplot2::geom_text(
    data = REST_random,
    mapping = ggplot2::aes(x = 0.1, y = SD_Est + 0.4, label = "REST Random"),
    inherit.aes = FALSE,
    color = fig_colors[3],
    size = 3
  ) +
  ggplot2::scale_colour_manual(name = "PATH Sample Design", values = model_colors) +
  # ggplot2::scale_linetype_manual(name = "PATH Sample Design", values = model_linetypes) +
  ggplot2::scale_x_continuous(
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c("10%", "20%", "30%", "40%", "50%")
  ) +
  ggplot2::scale_shape_manual(name = "PATH Sample Design", values = c(16, 15, 17, 18)) +
  ggplot2::labs(
    x = "Relative Effort",
    y = "Empirical SD"
  ) +
  ggplot2::annotate("text", x = Inf, y = Inf,
                    label = "b", hjust = 2, vjust = 1.5,
                    size = 5, fontface = "bold") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = c(0.85, 0.75),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    legend.key.width = ggplot2::unit(1.2, "cm"),
    legend.background = ggplot2::element_rect(color = "black", fill = "white", linewidth = 0.5),

    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# ggplot2::ggsave(
#   paste0(img_dir, "/cam_vs_emp_SD.pdf"),
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

# Plot 95% coverage for PATH with different cameras, IS and REST with 250 cameras
D_cam_means_fctr |>
  dplyr::filter(
    Bias %in% c("Random", "100%"),
    !(Model %in% c("IS Random", "REST Random"))
  ) |>
  dplyr::mutate(Percent_effort = cams / 250) |>
  ggplot2::ggplot(ggplot2::aes(x = Percent_effort, y = coverage_prob, color = Design, shape = Design)) +
  ggplot2::geom_hline(yintercept = IS_random$coverage_prob, linetype = "dashed", linewidth = 0.7, color = fig_colors[1]) +
  ggplot2::geom_hline(yintercept = REST_random$coverage_prob, linetype = "dashed", linewidth = 0.7, color = fig_colors[3]) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    data = IS_random,
    mapping = ggplot2::aes(x = 0.1, y = coverage_prob + 0.05, label = "IS Random"),
    inherit.aes = FALSE,
    color = fig_colors[1],
    size = 3
  ) +
  ggplot2::geom_text(
    data = REST_random,
    mapping = ggplot2::aes(x = 0.1, y = coverage_prob + 0.05, label = "REST Random"),
    inherit.aes = FALSE,
    color = fig_colors[3],
    size = 3
  ) +
  ggplot2::scale_colour_manual(name = "PATH Sample Design", values = model_colors) +
  # ggplot2::scale_linetype_manual(name = "PATH Sample Design", values = model_linetypes) +
  ggplot2::scale_x_continuous(
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c("10%", "20%", "30%", "40%", "50%")
  ) +
  ggplot2::scale_shape_manual(name = "PATH Sample Design", values = c(16, 15, 17, 18)) +
  ggplot2::labs(
    x = "Relative Effort",
    y = "Coverage Probability"
  ) +
  ggplot2::annotate("text", x = Inf, y = Inf,
                    label = "b", hjust = 2, vjust = 1.5,
                    size = 5, fontface = "bold") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = c(0.85, 0.75),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    legend.key.width = ggplot2::unit(1.2, "cm"),
    legend.background = ggplot2::element_rect(color = "black", fill = "white", linewidth = 0.5),

    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# ggplot2::ggsave(
#   paste0(img_dir, "/cam_vs_emp_SD.pdf"),
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

# Plot credible intervals for PATH with different cameras, IS and REST with 250 cameras
D_cam_means_fctr |>
  dplyr::filter(
    Bias %in% c("Random", "100%"),
    !(Model %in% c("IS Random", "REST Random"))
  ) |>
  dplyr::mutate(Percent_effort = cams / 250) |>
  ggplot2::ggplot(ggplot2::aes(x = Percent_effort, y = Mean_CrI, color = Design, shape = Design)) +
  ggplot2::geom_hline(yintercept = IS_random$Mean_CrI, linetype = "dashed", linewidth = 0.7, color = fig_colors[1]) +
  ggplot2::geom_hline(yintercept = REST_random$Mean_CrI, linetype = "dashed", linewidth = 0.7, color = fig_colors[3]) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    data = IS_random,
    mapping = ggplot2::aes(x = 0.1, y = Mean_CrI + 2, label = "IS Random"),
    inherit.aes = FALSE,
    color = fig_colors[1],
    size = 3
  ) +
  ggplot2::geom_text(
    data = REST_random,
    mapping = ggplot2::aes(x = 0.1, y = Mean_CrI + 2, label = "REST Random"),
    inherit.aes = FALSE,
    color = fig_colors[3],
    size = 3
  ) +
  ggplot2::scale_colour_manual(name = "PATH Sample Design", values = model_colors) +
  # ggplot2::scale_linetype_manual(name = "PATH Sample Design", values = model_linetypes) +
  ggplot2::scale_x_continuous(
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c("10%", "20%", "30%", "40%", "50%")
  ) +
  ggplot2::scale_shape_manual(name = "PATH Sample Design", values = c(16, 15, 17, 18)) +
  ggplot2::labs(
    x = "Relative Effort",
    y = "Credible Interval Width"
  ) +
  ggplot2::annotate("text", x = Inf, y = Inf,
                    label = "c", hjust = 2, vjust = 1.5,
                    size = 5, fontface = "bold") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = c(0.85, 0.75),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    legend.key.width = ggplot2::unit(1.2, "cm"),
    legend.background = ggplot2::element_rect(color = "black", fill = "white", linewidth = 0.5),

    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# ggplot2::ggsave(
#   paste0(img_dir, "/cam_vs_cri.pdf"),
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
################################################################################
################################################################################
# # Plot effort vs MAE for IS method at 250 cameras
# D_cam_means_fctr |>
#   dplyr::filter(
#     Bias %in% c("Random", "100%")
#   ) |>
#   # dplyr::filter(Model != "IS Random") |>
#   dplyr::mutate(Percent_effort = cams / 250) |>
#   ggplot2::ggplot(ggplot2::aes(x = Percent_effort, y = Mean_MAE, linetype = Model, color = Model)) +
#   ggplot2::geom_hline(yintercept=IS_random$Mean_MAE, linetype="dashed", linewidth = 0.7, color = fig_colors[1]) +
#   ggplot2::geom_hline(yintercept=REST_random$Mean_MAE, linetype="dashed", linewidth = 0.7, color = fig_colors[3]) +
#   ggplot2::geom_text(
#     data = annotation_IS,
#     mapping = ggplot2::aes(x = x, y = y, label = label),
#     inherit.aes = FALSE,
#     size = 2
#   ) +
#   ggplot2::geom_text(
#     data = annotation_REST,
#     mapping = ggplot2::aes(x = x, y = y, label = label),
#     inherit.aes = FALSE,
#     size = 2
#   ) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   # custom_colors +
#   ggplot2::scale_colour_manual(name = "Model", values = model_colors) +
#   ggplot2::scale_linetype_manual(name = "Model", values = model_linetypes) +
#   ggplot2::scale_linewidth_manual(name = "Model", values = model_linewidths) +
#   ggplot2::guides(linetype=ggplot2::guide_legend(title="Model and Sample Design"),
#                   color = ggplot2::guide_legend(title="Model and Sample Design")) +
#   ggplot2::labs(x = "Relative Effort",
#                 y = "Mean Error") +
#   ggplot2::scale_x_continuous(breaks = c(.1,.2, .3, .4, .5),
#                               labels = c("10%","20%", "30%", "40%", "50%")) +
#   ggplot2::annotate("text", x = Inf, y = Inf,
#                     label = "a", hjust = 2, vjust = 1.5,
#                     size = 5) +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = "none",# c(0.85, 0.72),
#     legend.title = ggplot2::element_text(size=13),
#     legend.text = ggplot2::element_text(size = 13),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
#     # legend.background = ggplot2::element_blank(), # ggplot2::element_rect(color = "white"), #
#     # legend.spacing.y = ggplot2::unit(0, "mm"),
#     # legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# # ggplot2::ggsave(
# #   paste0(img_dir, "/effort_vs_MAE.pdf"),
# #   plot = ggplot2::last_plot(),
# #   # path = file_path,
# #   # scale = 1,
# #   width = 5,
# #   height = 3,
# #   # units = c("in", "cm", "mm", "px"),
# #   dpi = 600,
# #   limitsize = TRUE,
# #   bg = NULL
# # )
#
# # Plot effort vs precision for IS method at 250 cameras
# D_cam_means_fctr |>
#   dplyr::filter(
#     Bias %in% c("Random", "100%")
#   ) |>
#   # dplyr::filter(Model != "IS Random") |>
#   dplyr::mutate(Percent_effort = cams / 250) |>
#   ggplot2::ggplot(ggplot2::aes(x = Percent_effort, y = Mean_SD, color = Model, linetype = Model)) +
#   ggplot2::geom_hline(yintercept=IS_random$Mean_SD, linetype="dashed", size = 0.7) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   # custom_colors +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "solid",
#       "PATH Random" = "solid",
#       "REST (Cov) Random" = "solid",
#       "PATH 80% High" = "dashed",
#       "REST (Cov) 80% High" = "dashed",
#       "PATH 100% High" = "dotted",
#       "PATH 80% Moderate" = "dashed",
#       "REST (Cov) 80% Moderate" = "dashed",
#       "PATH 100% Moderate" = "dotted",
#       "PATH 80% Low" = "dashed",
#       "REST (Cov) 80% Low" = "dashed",
#       "PATH 100% Low" = "dotted"
#     )
#   ) +
#   ggplot2::labs(x = "Relative Effort",
#                 y = "Mean Variance") +
#   ggplot2::scale_x_continuous(breaks = c(.1,.2, .3, .4, .5),
#                               labels = c("10%","20%", "30%", "40%", "50%")) +
#   ggplot2::annotate("text", x = Inf, y = Inf,
#                     label = "b", hjust = 2, vjust = 1.5,
#                     size = 5) +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = "none",# c(0.85, 0.72),
#     legend.title = ggplot2::element_text(size=13),
#     legend.text = ggplot2::element_text(size = 13),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_blank(), # ggplot2::element_rect(color = "white"), #
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   ) +
#   facet_wrap(~ Design)
#
# # ggplot2::ggsave(
# #   paste0(img_dir, "/effort_vs_Var.pdf"),
# #   plot = ggplot2::last_plot(),
# #   # path = file_path,
# #   # scale = 1,
# #   width = 5,
# #   height = 3,
# #   # units = c("in", "cm", "mm", "px"),
# #   dpi = 600,
# #   limitsize = TRUE,
# #   bg = NULL
# # )
#
# ################################################################################
# # baseline_val <- IS_random$Mean_SD
#
# # Add the Strata column to your main dataframe
# plot_data <- D_cam_means |>
#   dplyr::filter(Model == "PATH" | (Model %in% c("IS", "REST") & Design == "Random")) |>
#   dplyr::mutate(
#     Percent_effort = cams / 250
#   )
#
# # Add the same Strata column to your highlight points
# # crossing_points <- plot_data |>
# #   dplyr::filter(Mean_SD < baseline_val) |>
# #   dplyr::group_by(Model_char) |>
# #   dplyr::slice_min(order_by = Percent_effort, n = 1) |>
# #   dplyr::ungroup()
#
#
# plot_data |>
#   ggplot2::ggplot(ggplot2::aes(x = Percent_effort, y = Mean_MAE, color = Model, linetype = Bias)) +
#   ggplot2::geom_hline(yintercept = IS_random$Mean_MAE, linetype = "dashed", size = 0.7) +
#   ggplot2::geom_hline(yintercept = REST_random$Mean_MAE, linetype = "dashed", color = "red", size = 0.7) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::facet_wrap(~ Speed) +
#   ggplot2::scale_color_manual(values = fig_colors[1:5]) +
#   ggplot2::labs(x = "Relative Effort",
#                 y = "Mean Error") +
#   ggplot2::scale_x_continuous(breaks = c(.1, .2, .3, .4, .5),
#                               labels = c("10%", "20%", "30%", "40%", "50%")) +
#   ggplot2::geom_text(
#     data = annotation_IS,
#     mapping = ggplot2::aes(x = x, y = y, label = label),
#     inherit.aes = FALSE,
#     x = 2,
#     y = IS_random$Mean_SD + 0.8,
#     size = 3
#   ) +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  legend.title=element_text(size=10),
#                  legend.text=element_text(size=9),
#                  legend.position = c(0.092, 0.793),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  panel.spacing = unit(0,'lines'),
#                  axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth = 0.5),
#                  legend.background = ggplot2::element_rect(color = "black"),
#                  legend.spacing.y = ggplot2::unit(0, "mm"),
#                  legend.box.background = ggplot2::element_rect(colour = "black"))
################################################################################
# tmp_pal <- c("black", "#1B5E20", "#00A8C6", "#4B6FAD", "#8E44AD")
# # fixed_tag_palette <- setNames(tmp_pal[seq_along(c("IS", "PATH"))], c("IS", "PATH"))
#
# # Plot sds
# D_cam_means %>%
#   dplyr::filter(Design == "Random") %>%
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c("IS Random" = "dashed", "PATH Random" = "solid")
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "a",
#                     size = 5) +
#   ggplot2::ggtitle("Random Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_cam_means %>%
#   dplyr::filter(Design %in% c("80% High", "100% High") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "dashed",
#       "PATH 80% High" = "solid",
#       "PATH 100% High" = "solid"
#     )
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "d",
#                     size = 5) +
#   ggplot2::ggtitle("High Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_cam_means %>%
#   dplyr::filter(Design %in% c("80% Moderate", "100% Moderate") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "dashed",
#       "PATH 80% Moderate" = "solid",
#       "PATH 100% Moderate" = "solid"
#     )
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "c",
#                     size = 5) +
#   ggplot2::ggtitle("Moderate Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_cam_means %>%
#   dplyr::filter(Design %in% c("80% Low", "100% Low") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(
#     ggplot2::aes(x = cams, y = Mean_SD, color = Model, linetype = Model)
#   ) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::scale_linetype_manual(
#     values = c(
#       "IS Random" = "dashed",
#       "PATH 80% Low" = "solid",
#       "PATH 100% Low" = "solid"
#     )
#   ) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   ggplot2::scale_fill_manual(values= tmp_pal) +
#   ggplot2::scale_color_manual(values = tmp_pal) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Low Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )

################################################################################
################################################################################
# # Try boxplots
# D_all_alt <- D_all %>%
#   dplyr::filter(Model == "PATH") |>
#   dplyr::bind_rows(
#     D_all |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(MAE = abs(tot_animals - Est)) |>
#   dplyr::mutate(Model = paste(Model, Design))
#
# D_all_alt$Model <- factor(
#   D_all_alt$Model,
#   levels = c("IS Random", "PATH Random",
#              "PATH 80% High", "PATH 100% High",
#              "PATH 80% Moderate", "PATH 100% Moderate",
#              "PATH 80% Low", "PATH 100% Low")
# )
#
# D_all_alt |>
#   ggplot2::ggplot(ggplot2::aes(x = as.factor(cams), y = SD, fill = Model)) +
#   ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior SDs") +
#   # ggplot2::scale_fill_manual(values= fig_colors[1:5]) +
#   # ggplot2::scale_color_manual(values = c('grey0', 'grey40', 'grey60')) +
#   # ggplot2::geom_hline(
#   #   yintercept = D_cam_means |>
#   #     dplyr::filter(Model == "IS" & Design == "Random") |>
#   #     dplyr::pull(Mean_SD),
#   #   linetype="dashed",
#   #   size = 0.7
#   # ) |>
#   ggplot2::annotate("text", x = 0.5, y = 32,
#                     label = "b",
#                     size = 5) +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  legend.title=element_text(size=12),
#                  legend.text=element_text(size=10),
#                  legend.position = c(0.9, 0.77),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#                  legend.background = ggplot2::element_rect(color = "white"),
#                  legend.spacing.y = ggplot2::unit(0, "mm"),
#                  legend.box.background = ggplot2::element_rect(colour = "black"))
#
# ################################################################################
# # Plot means
# D_cam_means %>%
#   dplyr::filter(Design == "Random") %>%
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Random Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_cam_means %>%
#   dplyr::filter(Design %in% c("80% High", "100% High") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("High Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_cam_means %>%
#   dplyr::filter(Design %in% c("80% Moderate", "100% Moderate") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Moderate Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_cam_means %>%
#   dplyr::filter(Design %in% c("80% Low", "100% Low") & Model %in% "PATH") %>%
#   dplyr::bind_rows(
#     D_cam_means |>
#       dplyr::filter(Model == "IS" & Design == "Random")
#   ) |>
#   dplyr::mutate(Model = paste(Model, Design)) |>
#   ggplot2::ggplot(ggplot2::aes(x = cams, y = Mean_Est, color = Model)) + #, shape = Speed)) +
#   ggplot2::geom_line(size = 1) +
#   ggplot2::geom_point(size = 2) +
#   ggplot2::labs(x = "Number of Cameras",
#                 y = "Posterior Means") +
#   ggplot2::scale_fill_manual(values= fig_colors) +
#   ggplot2::scale_color_manual(values = fig_colors) +
#   ggplot2::annotate("text", x = 8, y = 80,# max(D_cam_means$Mean_SD, na.rm = T),
#                     label = "b",
#                     size = 5) +
#   ggplot2::ggtitle("Low Habitat Design") +
#   ggplot2::theme(
#     axis.title=element_text(size = 16),
#     axis.text = ggplot2::element_text(size = 16),
#     legend.position = c(0.85, 0.72),
#     legend.text = ggplot2::element_text(size = 8),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
#     legend.background = ggplot2::element_rect(color = "white"),
#     legend.spacing.y = ggplot2::unit(0, "mm"),
#     legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
#
################################################################################
# # Maybe put in a table - PR requires more data than PD and HAM
# D_omit <- D_all %>%
#   dplyr::group_by(Model, cams, Design) %>%
#   dplyr::summarise(
#     count = sum(Est > 250, na.rm = T) / 1000,
#     .groups = 'drop'
#   ) %>%
#   dplyr::mutate(cams = as.character(cams))
#
# D_omit$cams <- factor(
#   D_omit$cams,
#   levels = c("10", "25", "50", "75", "100")
# )
#
# D_omit %>%
#   ggplot2::ggplot(ggplot2::aes(x = Design, y = count, fill = cams)) +
#   ggplot2::geom_col(position = "dodge") +
#   ggplot2::facet_wrap(~ Model, ncol = 1) +
#   ggplot2::labs(
#     x = "Sampling Design",
#     y = "Proportion of Simulations Omitted"
#   ) +
#   ggplot2::scale_fill_manual(values= fig_colors, name = "Number of Cameras") +
#   ggplot2::theme(text = ggplot2::element_text(size = 16),
#                  legend.title = ggplot2::element_blank(),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
#                  # legend.position = "none",
#                  # legend.background = ggplot2::element_blank(),
#                  # legend.spacing.y = ggplot2::unit(0, "mm"),
#                  # legend.box.background = ggplot2::element_rect(colour = "black")
#   )
#
# D_omit %>%
#   # dplyr::filter(!(Design %in% c("100% Low", "100% Moderate", "100% High"))) %>%
#   ggplot2::ggplot(ggplot2::aes(x = Design, y = count, fill = Model)) +
#   ggplot2::geom_col(position = "dodge", color = 'black') +
#   ggplot2::facet_wrap(~ cams, ncol = 1) +
#   ggplot2::labs(
#     x = "Sampling Design",
#     y = "Proportion of Simulations Omitted"
#   ) +
#   ggplot2::scale_fill_manual(values= fig_colors, name = "Number of Cameras") +
#   ggplot2::scale_color_manual(values= fig_colors, name = "Number of Cameras") +
#   # ggplot2::ylim(0,0.2) +
#   ggplot2::scale_y_continuous(limits = c(0,0.2), expand = c(0,0)) +
#   ggplot2::theme(axis.title=element_text(size = 16),
#                  axis.text = ggplot2::element_text(size = 8),
#                  legend.title = ggplot2::element_blank(),
#                  panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank(),
#                  panel.background = ggplot2::element_blank(),
#                  axis.line = ggplot2::element_line(colour = "black"),
#                  panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
#                  # legend.position = "none",
#                  # legend.background = ggplot2::element_blank(),
#                  # legend.spacing.y = ggplot2::unit(0, "mm"),
#                  # legend.box.background = ggplot2::element_rect(colour = "black")
#   )

# ggplot2::ggsave(
#   paste0(img_dir, "Omitted_counts.png"),
#   plot = ggplot2::last_plot(),
#   # path = file_path,
#   # scale = 1,
#   width = 7,
#   height = 5,
#   # units = c("in", "cm", "mm", "px"),
#   dpi = 600,
#   limitsize = TRUE,
#   bg = NULL
# )

