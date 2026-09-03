library(dplyr)
library(ggplot2)
tot_animals <- 25

sim_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_cam_results/"
fig_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/imgs/"
save_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/D_all_results/"

fig_colors <- c("#1B5E20", "#00A8C6", "#FBC02D", "#E65100", "#8E44AD", "#4B6FAD", "#D81B60")
################################################################################
ncams <- 100

cam_designs <- c("Random", "Ag_all")
ABM_designs <- c("Home Range", "Correlated", "Hybrid")
# Collect density estimates for from all results
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

source("./R/utils.R")
source("./R/plot_funs.R")


D_all_NLCD <- tibble::tibble()
NLCD_dat <- tibble::tibble()
for (ii in 1:length(cam_designs)) {
  for (jj in 1:length(ABM_designs)) {
    print(ii)
    all_results <- loadRData(paste0(sim_dir,
                                    cam_designs[ii],
                                    "_",
                                    ncams,
                                    "_cam_NLCD",
                                    ABM_designs[jj],
                                    ".RData")
    )

    D_all_NLCD <- dplyr::bind_rows(
      D_all_NLCD,
      all_results[[5]] %>%
        dplyr::bind_rows() |>
        dplyr::mutate(
          SampDesign = paste0(cam_designs[ii], "_cam"),
          ABM_Design = ABM_designs[jj]
        )
    )

    NLCD_dat <- dplyr::bind_rows(
      NLCD_dat,
      all_results[[4]] %>%
        dplyr::bind_rows() |>
        dplyr::mutate(
          SampDesign = paste0(cam_designs[ii], "_cam"),
          ABM_Design = ABM_designs[jj]
        )
    )

  }
}

D_all_NLCD |>
  dplyr::count(Model, Covariate, SampDesign, ABM_Design)

# save(D_all_NLCD, file = paste0(
#   save_dir,
#   "D_all_NLCD.RData")
# )

#--------------------------------------------------

load(paste0(
  save_dir,
  "D_all_NLCD.RData")
)

D_all <- D_all_NLCD %>%
  dplyr::filter(Est > 0) |> #Removes IS method with no counts
  dplyr::mutate(
    Model = ifelse(
      SampDesign == "Ag_all_cam",
      paste(Model, "Ag Bias"),
      Model
    ),
    ABM_Design = ifelse(
      ABM_Design == "Correlated",
      "CRW",
      ABM_Design
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

D_cam_means <- D_all |>
  dplyr::group_by(cams, Model, design, ABM_Design) %>%
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

################################################################################
# Plot Mean estimates
D_all %>%
  dplyr::filter(Model %in% c("PATH", "IS", "REST", "PATH Ag Bias")) |>
  dplyr::filter(Est < 250) |>
  dplyr::mutate(Model = paste(Model, Covariate)) |>
  ggplot2::ggplot(ggplot2::aes(x = Model, y = Est, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::geom_hline(yintercept=tot_animals, linetype="dashed", size = 0.7) +
  ggplot2::labs(x = "Model",
                y = "Posterior Mean") +
  ggplot2::facet_grid(~ ABM_Design) +
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
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
#          "NLCD_mean_est.pdf"),
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
# Plot Mean estimates
D_all %>%
  dplyr::filter(Model %in% c("PATH", "IS", "REST", "PATH Ag Bias")) |>
  ggplot2::ggplot(ggplot2::aes(x = Model, y = CrI, fill = Model)) +
  ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
  ggplot2::labs(x = "Model",
                y = "Credible Interval Width") +
  ggplot2::facet_grid(~ ABM_Design) +
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
#          "NLCD_cri.pdf"),
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
# annotation_df_b <- data.frame(
#   x = Inf,
#   y = Inf,
#   label = "b",
#   BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
#   # Percentage = factor("100%", levels =  c("80%", "100%"))
# )
# annotation_IS <- data.frame(
#   x = Inf,
#   y = Inf,
#   label = "IS Variance",
#   BiasLevel = factor("High", levels = c("Low", "Moderate", "High"))
#   # Percentage = factor("100%", levels =  c("80%", "100%"))
# )
#
#
# D_all_separated %>%
#   # dplyr::filter(SampDesign != "Random" & Model == "PATH") |>
#   ggplot2::ggplot(ggplot2::aes(x = Percentage, y = SD, fill = Model)) +
#   ggplot2::geom_boxplot(lwd = 0.5, fatten = .5, outlier.shape = NA) +
#   ggplot2::geom_hline(yintercept=IS_random$Mean_var, linetype="dashed", size = 0.7) +
#   ggplot2::labs(x = "Sampling Bias",
#                 y = "Posterior Variance") +
#   ggplot2::scale_y_continuous(limits = c(0, 400)) +
#   ggplot2::facet_grid(~ SampDesign) #+ #, switch = "x") +
#   # ggplot2::guides(
#   #   fill = ggplot2::guide_legend(
#   #     title = "Model",
#   #     override.aes = list(
#   #       linetype = c("blank", "dashed"),
#   #       fill = c("#00A8C6", NA)
#   #     )
#   #   )
#   # ) +
#   # ggplot2::geom_text(
#   #   data = annotation_df_b,
#   #   mapping = ggplot2::aes(x = x, y = y, label = label),
#   #   inherit.aes = FALSE,
#   #   hjust = 2.5,   # Horizontal adjustment (same as before)
#   #   vjust = 1.5,   # Vertical adjustment (same as before)
#   #   size = 5
#   # ) +
#   # ggplot2::geom_text(
#   #   data = annotation_IS,
#   #   mapping = ggplot2::aes(x = x, y = y, label = label),
#   #   inherit.aes = FALSE,
#   #   x = 2,
#   #   y = IS_random$Mean_var + 0.8,
#   #   size = 3
#   # ) +
#   # # annotate("text",
#   # #          x = 6,
#   # #          y = IS_random$Mean_var + 0.7,
#   # #          size = 3,
#   # #          label = "IS Variance",
#   # #          color = "black") +
#   # ggplot2::theme(text = ggplot2::element_text(size = 16),
#   #                legend.title=element_text(size=10),
#   #                legend.text=element_text(size=9),
#   #                legend.position = c(0.092, 0.793),
#   #                panel.grid.major = ggplot2::element_blank(),
#   #                panel.grid.minor = ggplot2::element_blank(),
#   #                panel.background = ggplot2::element_blank(),
#   #                panel.spacing = unit(0,'lines'),
#   #                axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
#   #                panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth = 0.5),
#   #                legend.background = ggplot2::element_rect(color = "black"),
#   #                legend.spacing.y = ggplot2::unit(0, "mm"),
#   #                legend.box.background = ggplot2::element_rect(colour = "black"))
#
#
# # ggplot2::ggsave(
# #   paste0(fig_dir,
# #          "PATH_bias_cam_SD.pdf"),
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
# #--------------------------------------------------
# # Try all in one plot?
# D_all %>%
#   # dplyr::filter(Est < 250) %>%
#   dplyr::mutate(
#     SampDesign = dplyr::case_when(
#       SampDesign == "Ag_bias_cam" ~ "80% Agriculture",
#       SampDesign == "Ag_all_cam" ~ "100% Agriculture",
#       .default = "Random"
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
#       SampDesign == "Ag_bias_cam" ~ "80% Agriculture",
#       SampDesign == "Ag_all_cam" ~ "100% Agriculture",
#       .default = "Random"
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
