################################################################################
# Custom abm simulation parameters
tot_N <- 2

# Hybrid correlated walk / home range (60/40 split)
sim_name <- "Hybrid"
home_range_strength <- list(stats::runif(tot_N, 0.0005, 0.001))
init_placement <- list(c(1, 0, 0))
corr_strength <- 0

################################################################################

library(tidyverse)
library(RColorBrewer)
library(lattice)
library(gridExtra)
library(doParallel)
devtools::load_all()

sim_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_cam_results/"
fig_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/imgs/"

# Initializations
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")
options(ggplot2.discrete.colour = fig_colors)
options(ggplot2.discrete.fill = fig_colors)

################################################################################
# Load NLCD data set
# tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5000NLCDclip.tif"
tif_filename <- "/home/guengrosklos/Desktop/NLCD_data/NLCD_data/LowTag5010NLCDclip.tif"

# Time step for ABM
t_step_size <- 0.25

mu_base <- tibble::tibble(
  LandCover = c("Water", "Development", "Forest", "Agriculture"),
  speed = c(0, 360, 150, 90) # Manually set km/hr and convert to cell/hr
)

# custom_tiles <- tibble::tibble(
#   x = list(111:140),
#   y = list(301:330)
# )
custom_tiles <- tibble::tibble(
  x = list(141:230),
  y = list(299:388)
)

bg_info <- buildBackground(mu_base$Mu, tifFile = tif_filename, custom_tiles)

# Plot
# Convert the matrix to a long-format data frame
df <- reshape2::melt(bg_info$Landscape)
colnames(df) <- c("Row", "Column", "LandCover")

# # Reassign ag to water in isolated section
# df$LandCover[18] <- 1
# Convert the numeric values (1-4) into categorical factors with labels
df$LandCover <- factor(df$LandCover,
                       levels = c(1, 2, 3, 4),
                       labels = c("Water", "Development", "Forest", "Agriculture"))

# Join mu to each cover type
df <- df |>
  dplyr::left_join(
    mu_base,
    by = "LandCover"
  ) |>
  tibble::as_tibble()

################################################################################
# Study design
study_design <- tibble::tibble(
  q = nrow(df), # Number grid cells
  dx = 30,  # Grid cell lengths (m)
  dy = 30,
  t_steps = 500 / t_step_size, # Number of time steps
  dt = t_step_size, # Time step size (hr)
  t_censor = 100,
  bounds = list(c(0, dx * q ^ 0.5)), # Sampling area boundaries
  tot_A = (bounds[[1]][2] - bounds[[1]][1])^2,
  num_groups = tot_N,
  group_sizes = list(rep(1, num_groups)),
  group_spread = 0, # Tightness of grouping behavior (relative to grid size)
  h_range_strength = home_range_strength,
  tot_animals = sum(unlist(group_sizes)),
  Initial_placement = init_placement,
  corr_strength = corr_strength,
  # MCMC parms
  num_runs = 1,
  n_iter = 40000,
  burn_in = 30000,
  covariate_labels = list(c("Agriculture", "Development", "Forest")) # don't include restricted habitats
)

# Landscape design
# Motility rates for the hab types are 4, 2, 0.2, 0.5 for water, development, forest, and agriculture
lscape_design <- tibble::tibble(
  lscape_tag = "Custom", #"Random", #
  Speed_ID = c("Water", "Development", "Forest", "Agriculture"),
  Speed_mean = mu_base$speed,
  Speed_sd = Speed_mean / 10,
  gamma_rate = (Speed_mean / Speed_sd) ^ 2,
  gamma_shape = Speed_mean / Speed_sd ^ 2,
  Speed_mins = mu_base$speed * 0.5,
  Speed_maxes = mu_base$speed * 1.5
) |>
  dplyr::arrange(Speed_ID)

# # Plot gamma distributions for speeds in each habitat type
# lscape_design |>
#   dplyr::filter(Speed_ID != "Water") |>
#   dplyr::group_by(Speed_ID) |>
#   dplyr::mutate(
#     speed_plot = list(rgamma(10000, gamma_rate, gamma_shape))
#   ) |>
#   tidyr::unnest(speed_plot) |>
#   ggplot2::ggplot(ggplot2::aes(x = speed_plot, fill = Speed_ID)) +
#   ggplot2::geom_density()

# Define lscape for camera simulation
lscape_defs <- df |>
  tibble::as_tibble() |>
  dplyr::select(X = Column, Y = Row, Speed = LandCover) |>
  dplyr::mutate(
    Y = sqrt(study_design$q) + 1 - Y
  ) |>
  dplyr::arrange(X, Y) |>
  dplyr::mutate(
    Index = 1:study_design$q
  ) |>
  dplyr::left_join(
    lscape_design |>
      dplyr::select(Speed = Speed_ID, Min = Speed_mins, Max = Speed_maxes, gamma_rate, gamma_shape),
    by = "Speed"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    Value = runif(
      1,
      Min,
      Max
    )
  )

# Create covariate matrix with 0, 1 values
# Recalculate total area to exclude Water cells
study_design <- study_design %>%
  dplyr::mutate(
    tot_A = dx * dy * (sum(lscape_defs$Speed != "Water")),
    num_covariates = length(unlist(covariate_labels)),
    Z = list(create_covariate_mat(
      lscape_defs,
      study_design,
      unlist(covariate_labels)))
  )

all_designs <- tibble::tibble(
  Design_name = c("Random", "Forest_all", "Ag_all"),
  Design =c("Random", "Bias", "Bias"),
  Props = c(
    list(c(1, 1, 1)),
    list(c(0, 0, 1)),
    list(c(1, 0, 0))
  )
)

# Define study design for correlated walk animals
study_design_2 <- study_design |>
  dplyr::mutate(
    num_groups = 1,
    group_sizes = list(rep(1, num_groups)),
    h_range_strength = NULL,
    tot_animals = sum(unlist(group_sizes)),
    Initial_placement = NULL,
    corr_strength = 1
  )

# Adjust original study design
study_design <- study_design |>
  dplyr::mutate(
    num_groups = 1,
    group_sizes = list(rep(1, num_groups)),
    tot_animals = sum(unlist(group_sizes))
  )


# Run agent-based model
animalxy.all <- ABM_sim(study_design,
                        lscape_defs) |>
  dplyr::mutate(
    Strategy = "Home Range"
  )

animalxy.all.2 <- ABM_sim(study_design_2,
                          lscape_defs) |>
  dplyr::mutate(
    Animal_ID = Animal_ID + max(animalxy.all$Animal_ID),
    group_ID = group_ID + max(animalxy.all$group_ID),
    ii = ii + max(animalxy.all$ii),
    Strategy = "CRW"
  )

animalxy.all <- animalxy.all |>
  dplyr::bind_rows(animalxy.all.2)

# save(animalxy.all, file = paste0(sim_dir,"sample_scen_2_walks.RData"))

# load(file = paste0(sim_dir,"sample_scen_2_walks.RData"))

plot_ABM_2(study_design, lscape_defs, animalxy.all)

# ggplot2::ggsave(
#   paste0(fig_dir,
#          "NLCD_move_sample.pdf"),
#   plot = ggplot2::last_plot(),
#   width = 5,
#   height = 3,
#   dpi = 600,
#   limitsize = TRUE,
#   bg = "white"
# )

# Calculate the home range radius and area for each animal
# Shooting for home ranges around 0.5 area_km2 for home range animals
hr_metrics <- animalxy.all |>
  group_by(Animal_ID) |>
  mutate(
    # Find the spatial center (centroid) for the animal
    center_X = mean(X),
    center_Y = mean(Y),
    # Calculate the distance of every point to that center
    dist_to_center = sqrt((X - center_X)^2 + (Y - center_Y)^2)
  ) |>
  summarize(
    # Calculate the 95% core radius (removes outlier extreme steps)
    radius_95_m = quantile(dist_to_center, probs = 0.95),
    # Calculate absolute max radius for reference
    radius_max_m = max(dist_to_center),
    # Calculate the 95% home range area in square kilometers
    area_km2 = (pi * radius_95_m^2) / 1000000,

    .groups = "drop"
  )
