# library(dplyr)
library(ggplot2)
# source("./R/utils.R")
# source("./R/plot_funs.R")

save_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/convergence_results/"

fig_colors <- c("#1B5E20", "#00A8C6", "#FBC02D", "#E65100", "#8E44AD", "#4B6FAD", "#D81B60")

file_list <- list.files(path = save_dir,
                        pattern = "\\.RData$",
                        full.names = TRUE,
                        ignore.case = TRUE)

all_Rhat_ESS <- tibble::tibble()
if (length(file_list) == 0) {
  message("No .RData files found in the specified folder.")
} else {
  for (file in file_list) {
    message("Loading: ", basename(file))
    load(file)

    all_Rhat_ESS <- all_Rhat_ESS |>
      dplyr::bind_rows(save_conv_test)
  }
  message("All files loaded successfully!")
}

all_Rhat_ESS_format <- all_Rhat_ESS |>
  dplyr::filter(Model != "IS") |>
  dplyr::mutate(
    Model = dplyr::case_when(
      Covariate == "Covariate" & Model == "REST" ~ "REST (Cov)",
      Covariate == "Non-Covariate" & Model == "REST" ~ "REST (Non-Cov)",
      .default = Model
    )
  ) |>
  dplyr::filter(!(Model == "REST (Cov)" & cam_design %in% c("Medium_bias", "Slow_bias", "Fast_bias"))) |>
  dplyr::filter(
    Model == "PATH" |
      (Model == "REST (Non-Cov)" & cam_design == "Random") |
      Model == "REST (Cov)"
  ) |>
  # dplyr::filter(SampDesign != "Random" & !(Model %in% c("IS", "REST (Non-Cov)"))) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    min_ESS = min(unlist(ess)),
    is_converged = max_Rhat < 1.1 & min_ESS > 400,
  ) |>
  dplyr::ungroup()

Rhat_ESS_summary <- all_Rhat_ESS_format |>
  dplyr::group_by(cam_design, cams, Model, Covariate) |>
  dplyr::summarise(
    percent_converged = sum(is_converged, na.rm = T) / dplyr::n(),
    .groups = 'drop'
  )

Rhat_ESS_summary |>
  ggplot2::ggplot(ggplot2::aes(x = cams, y = percent_converged, linetype = Model, color = cam_design)) +
  ggplot2::geom_line()

ggplot(Rhat_ESS_summary, aes(x = cams, y = percent_converged, color = Model, group = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, aes(shape = Model)) +
  facet_wrap(~ cam_design, ncol = 3) + # Adjust ncol based on how it looks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  scale_x_continuous(breaks = c(25, 50, 75, 100, 125, 250)) +
  theme_bw() +
  labs(
    title = "MCMC Convergence Rates Across Sampling Designs",
    x = "Number of Cameras",
    y = "Convergence Rate (%)",
    color = "Model",
    shape = "Model"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ggplot(Rhat_ESS_summary, aes(x = as.factor(cams), y = paste(cam_design, Model, sep = " - "))) +
#   geom_tile(aes(fill = percent_converged), color = "white") +
#   scale_fill_viridis_c(option = "magma", limits = c(0, 1), name = "Convergence\nRate (%)") +
#   theme_minimal() +
#   labs(
#     title = "Convergence Heatmap",
#     x = "Number of Cameras",
#     y = "Design and Model Combination"
#   ) +
#   theme(
#     axis.text.y = element_text(size = 9),
#     panel.grid = element_blank()
#   )
