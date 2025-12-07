#' Plot Simulation Trajectories and Camera Locations
#'
#' @description
#' Generates a spatial plot showing the landscape grid, animal trajectories, and
#' camera trap polygons.
#'
#' @param study_design A list of simulation parameters including `bounds` and `dx`.
#' @param cam_design A list containing camera design details, specifically `ncam`.
#' @param cam_locs A data frame of camera locations/vertices.
#' @param animalxy.all A data frame of animal trajectories (output from `ABM_sim`).
#'
#' @return A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot geom_tile geom_path geom_rect geom_polygon labs guides theme element_blank coord_fixed aes scale_fill_manual
#' @export
#'
plot_ABM <- function(study_design,
                     cam_design,
                     cam_locs,
                     animalxy.all
                     ) {
  bounds <- unlist(study_design$bounds)
  ncam <- cam_design$ncam

  # Sample area border
  b.df <- data.frame(c(bounds, study_design$dx))

  # Camera locations
  dt.cam <- data.frame(
    group = rep(1:ncam, each = 3),
    polygon.x = unlist(cam_locs$x),
    polygon.y = unlist(cam_locs$y)
  )


  g <- ggplot() +
    # geom_tile(aes(X-.5*dx, Y-.5*dy, fill = Speed)) +
    scale_fill_manual(values = c("grey20", "grey50", "grey80")) +
    geom_path(
      data = animalxy.all,
      aes(X, Y, col = as.factor(Animal_ID))
    ) +
    geom_rect(
      data = b.df,
      aes(xmin = bounds[1],
          xmax = bounds[2],
          ymin = bounds[1],
          ymax = bounds[2]),
      color = "black",
      fill = NA
    ) +
    geom_polygon(
      data = dt.cam,
      aes(x = polygon.x,
          y = polygon.y,
          group = group),
      color = "black",
      fill = NA
    ) +
    labs(x = "X", y = "Y", fill = "Landscape Type") +
    guides(
      colour = "none",
      alpha = "none"
    ) +
    theme(panel.background = element_blank())
  g + coord_fixed()
}

#' Visualize Animal Space Use Intensity
#'
#' @description
#' rasterizes the animal trajectory points to create a heatmap (levelplot) of
#' space-use density across the landscape.
#'
#' @param study_design A list containing `q` (total cells) and `bounds`.
#' @param animalxy.all A data frame of animal trajectories.
#'
#' @return Displays a plot via `grid.arrange` (returns `NULL` or the grid object invisibly).
#'
#' @export
#'
plot_space_use <- function(study_design,
                           animalxy.all
                           ) {
  q <- study_design$q
  bounds <- unlist(study_design$bounds)

  u.abm.all <- matrix(0, nrow = q^0.5, ncol = q^0.5)

  for (xx in 1:nrow(animalxy.all)) {
    x.round <- ceiling(animalxy.all$X[xx] * (q^0.5 / max(bounds)))
    y.round <- ceiling(animalxy.all$Y[xx] * (q^0.5 / max(bounds)))

    # Transpose for converting matrix to raster
    u.abm.all[x.round, y.round] <- u.abm.all[x.round, y.round] + 1
  }

  u_plot <- levelplot(u.abm.all,
                      layers = 1,
                      main = list(expression("Animal space-use"), cex = 1),
                      cuts = 254,
                      margin = FALSE,
                      scales = list(draw = FALSE),
                      xlab = "x",
                      ylab = "y",
                      col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))

  grid.arrange(u_plot, ncol = 1, nrow = 1)
}

#' Plot Simulation Diagnostic Metrics
#'
#' @description
#' A suite of functions to visualize the performance of the model across multiple
#' runs. These plots compare estimates against truth (`tot_animals`) or compare
#' consistency across models.
#'
#' @param D.all A data frame containing simulation results. Must include columns:
#'   `Model`, `Est` (estimate), `SD`, and `Covariate`.
#' @param study_design (Optional) Used in `plot_multirun_means` to draw the
#'   true abundance dashed line (`tot_animals`).
#'
#' @return A `ggplot` object.
#'
#' @describeIn plot_simulation_diagnostics Boxplot of mean abundance estimates with a reference line for truth.
#' @export
plot_multirun_means <- function(study_design,
                                D.all) {

  ggplot(D.all, aes(x = Model, y = Est, fill = Covariate)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 outlier.shape = NA) +
    # ggplot2::geom_violin() +
    geom_hline(yintercept = study_design$tot_animals, linetype = "dashed", size = 1) +
    labs(
      x = "Model",
      y = "Mean Abundance"
    ) +
    coord_cartesian(ylim = quantile(D.all$Est, c(0.01, 0.99), na.rm = T)) +
    # scale_y_continuous(limits = c(lower, upper)) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.84),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

#' @describeIn plot_simulation_diagnostics Boxplot of standard deviations of the estimates.
#' @export
#'
plot_multirun_sds <- function(D.all) {
  ggplot(D.all, aes(x = Model, y = SD, fill = Covariate)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 outlier.shape = NA) +
    labs(
      x = "Model",
      y = "SD"
    ) +
    coord_cartesian(ylim = quantile(D.all$SD, c(0.01, 0.99), na.rm = T)) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.84),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

#' @describeIn plot_simulation_diagnostics Boxplot of Mean Absolute Percentage Error (MAPE).
#' @export
#'
plot_multirun_mape <- function(D.all, tot_N) {
  D.all <- D.all %>%
    dplyr::mutate(MAPE = abs(tot_N - Est) / Est)

  D.all %>%
    ggplot(aes(x = Model, y = MAPE, fill = Covariate)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 outlier.shape = NA) +
    labs(
      x = "Model",
      y = "MAPE"
    ) +
    coord_cartesian(ylim = quantile(D.all$MAPE, c(0.01, 0.99), na.rm = T)) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.84),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

#' @describeIn plot_simulation_diagnostics Boxplot of Coefficient of Variation (CV).
#' @export
#'
plot_multirun_CV <- function(D.all) {
  ggplot(D.all, aes(x = Model, y = SD / Est, fill = Covariate)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 outlier.shape = NA) +
    labs(
      x = "Model",
      y = "CV"
    ) +
    coord_cartesian(ylim = quantile(D.all$SD / D.all$Est, c(0.01, 0.99), na.rm = T)) +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = c(0.17, 0.84),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}


#' @describeIn plot_simulation_diagnostics Histogram of estimates with mean and SD lines.
#' @export
#'
plot_multirun_hist <- function(D.all) {
  mean_data <- D.all |>
    group_by(Model) |>
    summarise(
      Means = median(Est),
      SDs = sd(Est)
    )

  ggplot(D.all, aes(x = Est, fill = Model)) +
    # geom_density(position = "identity", alpha = 0.3, adjust = 1) +
    geom_histogram(position = "identity", alpha = 0.3, adjust = 1) +
    geom_vline(data = mean_data, aes(xintercept = Means, color = Model)) +
    labs(x = "Estimate", y = "Density", fill = "Model") +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "mm"),
      legend.box.background = element_rect(colour = "black")
    )
}

