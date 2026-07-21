# TODO:
#  - Run PATH, IS, and REST random, then PATH and REST with full and 80% concentration in high-dens habitat

library(tidyverse)
library(RColorBrewer)
library(lattice)
library(gridExtra)
library(doParallel)
devtools::load_all()

sim_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/sim_results_NLCD/"

# Initializations
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")
options(ggplot2.discrete.colour = fig_colors)
options(ggplot2.discrete.fill = fig_colors)

################################################################################
# Load NLCD data set
# tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5000NLCDclip.tif"
tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5010NLCDclip.tif"

mu_base <- tibble::tibble(
  LandCover = c("Water", "Development", "Forest", "Agriculture"),
  Mu = c(4, 2, 0.02, 0.5),
  # speed = 4 * Mu * 900 / 30 / 30 * 8,
  # speed = c(0, 2, 0.08, 0.5) * 30, # Manually set km/hr to cell/hr
  speed = c(0, 2, 0.5, 0.3) * 30, # Manually set km/hr and convert to cell/hr
  speed_km_hr = 4 * Mu * 900 / 30 / 1000
)

custom_tiles <- tibble::tibble(
  x = list(111:140),
  y = list(301:330)
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
# Run with different number of cameras
# cam_tests <- c(25, 50, 75, 100, 125)
# cam_tests <- c(50, 75, 100)
cam_tests <- c(100)

# tele_sample <- NULL
tele_sample <- tibble::tibble(
  t_sample_freq = 6,
  ID_sample_size = 5
)

# Study design
study_design <- tibble::tibble(
  q = 30^2, # Number grid cells
  dx = 30,  # Grid cell lengths (m)
  dy = 30,
  t_steps = 500, # Number of time steps
  dt = 1, # Time step size (hr)
  t_censor = 2,
  bounds = list(c(0, dx * q ^ 0.5)), # Sampling area boundaries
  tot_A = (bounds[[1]][2] - bounds[[1]][1])^2,
  num_groups = 10,
  group_sizes = list(rep(1, num_groups)),
  group_spread = 0, # Tightness of grouping behavior (relative to grid size)
  h_range_strength = list(stats::runif(num_groups, 0.0005, 0.01)),
  tot_animals = sum(unlist(group_sizes)),
  # Initial_placement = list(c(0.8, 0, 0.2)),
  Initial_placement = list(c(1, 0, 0)),
  # MCMC parms
  num_runs = 10,
  n_iter = 40000,
  burn_in = 30000,
  covariate_labels = list(c("Agriculture", "Development", "Forest")) # don't include restricted habitats
)

# Landscape design
# Motility rates for the hab types are 4, 2, 0.2, 0.5 for water, development, forest, and agriculture
lscape_design <- tibble::tibble(
  lscape_tag = "Custom", #"Random", #
  Speed_ID = c("Water", "Development", "Forest", "Agriculture"),
  # Speed_mins = c(0, 1.27, 0.12, 0.64), # Trying the sqrt of motility to start
  # Speed_maxes = c(0, 1.55, 0.155, .77)
  Speed_mins = mu_base$speed * 0.5,
  Speed_maxes = mu_base$speed * 1.5
  # Probs = c(0.1, 0.1, 0.8)
) |>
  dplyr::arrange(Speed_ID)

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
      dplyr::select(Speed = Speed_ID, Min = Speed_mins, Max = Speed_maxes),
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
  Design_name = c("Random", "Ag_bias", "Ag_all"),
  Design =c("Random", "Bias", "Bias"),
  Props = c(
    list(c(1, 1, 1)),
    list(c(0.8, 0, 0)),
    list(c(1, 0, 0))
  )
)

for (cam_des in 1:nrow(all_designs)) {
  for (cam in 1:length(cam_tests)) {

    # Cam designs
    # aim for cam_A ~ 50 - 100 m^2
    cam_design <- tibble::tibble(
      ncam = cam_tests[cam],
      snap_rate = 1 / 6, # snapshot rate (hours)
      Design_name = all_designs$Design_name[cam_des],
      Design = all_designs$Design[cam_des],
      Props = all_designs$Props[cam_des],
      cam_length = 0.02 * 30, #  study_design$dx  * 0.5, # length of all viewshed sides
      cam_A = cam_length ^ 2 / 2,
      tot_snaps = ncam * study_design$t_steps / snap_rate
    )

    # Initialize summary matrices
    D_all <- vector(mode = "list", length = study_design$num_runs * 2)
    D_all_REST <- vector(mode = "list", length = study_design$num_runs * 2)

    all_data <- tibble::tibble(
      iteration = 1:study_design$num_runs,
      cam_captures = NA,
      count_data = NA,
      encounter_data = NA,
      stay_time_all = NA,
      stay_time_data = NA
    )

    # Multi-run simulations
    for (run in 1:study_design$num_runs) {
      print(paste("Design", cam_des, "Cam", cam, "Run", run, "of",
                  study_design$num_runs))

      # Run agent-based model
      animalxy.all <- ABM_sim(study_design,
                              lscape_defs)

      tele_summary <- Collect_tele_data(animalxy.all, study_design, tele_sample = tele_sample) |>
        dplyr::arrange(Speed)

      # Use largest stay time as reference category
      ref_cat_idx <- which(tele_summary$stay_prop == min(tele_summary$stay_prop))

      # Set reference category for intercept
      study_design$Z[[1]][, ref_cat_idx] <- 1

      # Subtract reference category from stay time proportion
      # prop_adjust <- tele_summary$stay_prop /
      #   tele_summary$stay_prop[ref_cat_idx]
      # prop_adjust[ref_cat_idx] <- tele_summary$stay_prop[ref_cat_idx]
      # kappa.prior.mu.adj <- log(prop_adjust)
      kappa.prior.mu <- log(tele_summary$stay_prop)
      kappa.prior.var <- tele_summary$stay_sd^2 # stay_time_summary$cell_sd ^ 2

      # Place cameras on study area
      cam_locs <- create_cam_samp_design(study_design,
                                         lscape_defs,
                                         cam_design)

      ################################
      # Collect data
      ################################
      "%notin%" <- Negate("%in%")
      all_data$cam_captures[run] <- list(get_cam_captures(
        animalxy.all %>%
          dplyr::filter(t != 0),
        cam_locs,
        study_design
      ))

      habitat_summary <- lscape_defs %>%
        dplyr::filter(Speed != "Water") |>
        dplyr::group_by(Speed) %>%
        dplyr::summarise(
          n_lscape = dplyr::n() * study_design$dx * study_design$dy # Area of each habitat
        ) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(
          cam_locs %>%
            dplyr::group_by(Speed) %>%
            dplyr::summarise(
              ncams = dplyr::n(),
              .groups = 'drop'
            ),
          by = dplyr::join_by(Speed)
        ) %>%
        dplyr::left_join(
          tele_summary,
          by = dplyr::join_by(Speed)
        ) %>%
        dplyr::mutate(
          prop_cams = ncams / sum(ncams, na.rm = T),
          d_coeff = n_lscape * prop_cams / stay_prop / (cam_design$cam_A * study_design$t_steps)
        ) %>%
        replace(is.na(.), 0) %>%
        dplyr::select(Speed, n_lscape, prop_cams, d_coeff) |>
        dplyr::arrange(Speed)

      habitat_summary$Speed <- factor(
        habitat_summary$Speed,
        levels = unlist(study_design$covariate_labels)
      )

      seq_tbl <- tibble::tibble(
        val = seq(1, study_design$t_steps, by = cam_design$snap_rate)
      )

      count_data <- get_count_data(
        cam_locs,
        all_data$cam_captures[[run]],
        animalxy.all %>%
          dplyr::filter(t != 0),
        seq_tbl
      )

      all_data$count_data[[run]] <- list(count_data)

      encounter_data <- get_encounter_data(cam_locs, all_data$cam_captures[[run]])$encounter
      stay_time_data <- get_stay_time_data(cam_locs, all_data$cam_captures[[run]])[[2]] |>
        as.matrix()

      # Run models only if any data points were collected
      if (sum(count_data$count) == 0) {
        D.PR.MCMC.habitat <- NA
        SD.PR.MCMC.habitat <- NA
      } else {
        chain.PATH <- fit.model.mcmc.PATH(
          study_design = study_design,
          cam_design = cam_design,
          cam_locs = cam_locs,
          gamma_start = rep(log(mean(count_data$count)), study_design$num_covariates),
          gamma_prior_var = 10^4,
          gamma_tune = rep(-1, study_design$num_covariates),
          kappa_start = log(exp(kappa.prior.mu) / sum(exp(kappa.prior.mu))),
          kappa_prior_mu = kappa.prior.mu,
          kappa_prior_var = kappa.prior.var,
          kappa_tune = -1, #rep(-1, study_design$num_covariates),
          count_data_in = count_data,
          habitat_summary
        )

        ## Posterior summaries
        # plot(chain.PATH$tot_u[study_design$burn_in:study_design$n_iter])
        D.PATH.MCMC <- mean(chain.PATH$tot_u[study_design$burn_in:study_design$n_iter])
        SD.PATH.MCMC <- sd(chain.PATH$tot_u[study_design$burn_in:study_design$n_iter])

        if (any(colMeans(chain.PATH$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) || any(colMeans(chain.PATH$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
          warning(("Mean Count accept rate OOB"))
          D.PATH.MCMC <- NA
          SD.PATH.MCMC <- NA
        }
      }

################################################################################
      if (sum(encounter_data) == 0) {
        D.REST.MCMC <- NA
        SD.REST.MCMC <- NA
        D.REST.MCMC.cov <- NA
        SD.REST.MCMC.cov <- NA
      } else {
        # REST, no covariates
        chain.REST <- fit.model.mcmc.REST(
          study_design,
          cam_design,
          gamma_start = log(mean(encounter_data)),
          kappa_start = log(mean(stay_time_data,na.rm=T)),
          gamma_prior_var = 10^4,
          kappa_prior_var = 10^4,
          gamma_tune = -1,
          kappa_tune = -1,
          encounter_data_in = encounter_data,
          stay_time_data_in = stay_time_data
        )

        ## Posterior summaries
        # plot(chain.REST$tot_u[study_design$burn_in:study_design$n_iter])
        D.REST.MCMC <- mean(chain.REST$tot_u[study_design$burn_in:study_design$n_iter])
        SD.REST.MCMC <- sd(chain.REST$tot_u[study_design$burn_in:study_design$n_iter])

        if(any(colMeans(chain.REST$accept[study_design$burn_in:study_design$n_iter,])< 0.2) || any(colMeans(chain.REST$accept[study_design$burn_in:study_design$n_iter,])> 0.7)){
          warning(('REST accept rate OOB'))
          D.REST.MCMC <- NA
          SD.REST.MCMC <- NA
        }

        ###################################
        # REST w/ covariates
        ###################################
        chain.REST.cov <- fit.model.mcmc.REST.cov(
          study_design,
          cam_design,
          cam_locs,
          gamma_start = rep(log(mean(encounter_data)), study_design$num_covariates),
          kappa_start = rep(log(mean(stay_time_data,na.rm=T)), study_design$num_covariates),
          gamma_prior_var = 10^4,
          kappa_prior_var = 10^4,
          gamma_tune = rep(-1, study_design$num_covariates),
          kappa_tune = rep(-1, study_design$num_covariates),
          encounter_data_in = encounter_data,
          stay_time_data_in = stay_time_data
        )

        # ## Posterior summaries
        # plot(chain.REST.cov$tot_u[study_design$burn_in:study_design$n_iter])
        D.REST.MCMC.cov <- mean(chain.REST.cov$tot_u[study_design$burn_in:study_design$n_iter])
        SD.REST.MCMC.cov <- sd(chain.REST.cov$tot_u[study_design$burn_in:study_design$n_iter])

        if(any(colMeans(chain.REST.cov$accept[study_design$burn_in:study_design$n_iter,])< 0.2) ||
           any(colMeans(chain.REST.cov$accept[study_design$burn_in:study_design$n_iter,])> 0.7)){
          warning(('REST accept rate OOB'))
          D.REST.MCMC.cov <- NA
          SD.REST.MCMC.cov <- NA
        }

      }

      D_all[[(run - 1) * 2 + 1]] <- tibble::tibble(
        iteration = run,
        design = cam_design$Design_name,
        cams = cam_design$ncam,
        Model = "PATH",
        Covariate = "Non-Covariate",
        Est = D.PATH.MCMC,
        SD = SD.PATH.MCMC
      )

      D_all_REST[[(run - 1) * 2 + 1]] <- tibble::tibble(
        iteration = run,
        design = cam_design$Design_name,
        cams = cam_design$ncam,
        Model = "REST",
        Covariate = "Non-Covariate",
        Est = D.REST.MCMC,
        SD = SD.REST.MCMC
      )

      D_all_REST[[run * 2]] <- tibble::tibble(
        iteration = run,
        design = cam_design$Design_name,
        cams = cam_design$ncam,
        Model = "REST",
        Covariate = "Covariate",
        Est = D.REST.MCMC.cov,
        SD = SD.REST.MCMC.cov
      )

      ################################################################################
      # IS method
      tot_snaps <- study_design$t_steps / cam_design$snap_rate * cam_design$ncam
      IS_mean <- sum(count_data$count) * study_design$tot_A /
        (tot_snaps * cam_design$cam_A)

      M <- cam_design$ncam
      J <- study_design$t_steps
      L <- cam_design$cam_A * M * J
      sum_c <- sum((J * cam_design$cam_A) ^ 2 * (count_data$count /
                                                   (J * cam_design$cam_A) - sum(count_data$count) / L) ^ 2)

      IS_var <- M /(L^2 * (M - 1)) * sum_c

      form <- sprintf("~ %f * x1", study_design$tot_A)
      SE_N = msm::deltamethod(as.formula(form), IS_mean / study_design$tot_A, IS_var)

      D_all[[run * 2]] <- tibble::tibble(
        iteration = run,
        design = cam_design$Design_name,
        cams = cam_design$ncam,
        Model = "IS",
        Covariate = "Non-Covariate",
        Est = IS_mean,
        SD = SE_N
        # all_results = list(chain.PR.habitat)
      )

    }

    D_all <- dplyr::bind_rows(D_all) |>
      dplyr::bind_rows(D_all_REST)

    save_results <- list(
      study_design,
      cam_design,
      lscape_design,
      all_data,
      D_all
    )

    save(save_results, file = paste0(sim_dir,
                                     cam_design$Design_name,
                                     "_",
                                     cam_design$ncam,
                                     "_cam_NLCD.RData")
    )

    rm(save_results, all_data, D_all)
  }
}

# D_all <- dplyr::bind_rows(D_all) |>
#   dplyr::bind_rows(D_all_REST)
# # Omit all_results column (too much data)
# D_all <- D_all %>%
#   dplyr::select(-all_results)
#
# # Omit TTE data
# all_data <- all_data %>%
#   dplyr::select(-c(stay_time_all, TTE_data_all, TTE_data_raw, TTE_data))
#
# D_all$Model <- factor(D_all$Model, levels = c("TDST", "REST", "TTE", "PR", "IS", "STE", "SECR"))
#
# NA_summary <- D_all %>%
#   group_by(Model) %>%
#   summarise(num_NAs = sum(is.na(Est)))

###########################
# Plots
###########################
plot_multirun_means(study_design, D_all %>%
                      dplyr::filter(is.finite(Est)))
plot_multirun_sds(D_all %>%
                    dplyr::filter(is.finite(Est)))
# plot_multirun_mape(D_all %>%
#                       dplyr::filter(is.finite(Est)),
#                    study_design$tot_animals)
# plot_multirun_CV(D_all %>%
#                    dplyr::filter(is.finite(Est)))
# plot_multirun_hist(D_all)

# plot_ABM(study_design,
#          cam_design,
#          cam_locs,
#          animalxy.all)
plot_ABM_2(study_design, lscape_defs, animalxy.all)
# plot_ABM_2(study_design, lscape_defs, animalxy.all |> dplyr::filter(Animal_ID == 1))
plot_space_use(study_design,
               animalxy.all)
#
# plot_count_data(count_data_in = all_data$count_data[[1]],
#                 fill = "Speed")
# plot_encounter_data(encounter_data_in = all_data$encounter_data[[1]],
#                 fill = "Speed")
# plot_staytime_data(stay_time_raw_in = all_data$stay_time_raw[[1]],
#                    fill = "Speed")

# # Plot camera time series
# cam_caps <- all_data$cam_captures[[1]]
# cam_caps <- cam_caps %>%
#   dplyr::filter(in_cam == T) %>%
#     dplyr::mutate(t_hour = ceiling(t)) %>%
#   dplyr::group_by(t_hour) %>%
#     dplyr::arrange(t_hour) %>%
#   dplyr::summarise(
#     group_size = n(),
#     .groups = "drop"
#   )
#
# cam_caps <- cam_caps %>%
#   dplyr::bind_rows(tibble::tibble(
#     t_hour = which(seq(study_design$dt,
#                        study_design$t_steps * study_design$dt,
#                        by = study_design$dt) %notin% cam_caps$t_hour),
#       group_size = 0
#     ) ) %>%
#   arrange(t_hour)
#
# xy_cams <- unnest(cam_locs, cols = c(x, y, vertex)) %>%
#   group_by(lscape_index) %>%
#   dplyr::mutate(
#     x = x[1],
#     y = y[1]
#   ) %>%
#   dplyr::select(cam_ID, lscape_index, x, y) %>%
#   distinct()

# # Plot spatial cam captures
# cd <- all_data$count_data[[1]]
#
# cd %>%
#   ggplot() +
#   geom_point(aes(xy_cams$x, xy_cams$y),
#              color = "black",
#              # shape = 2,
#              size = cd$count*5,
#              stroke = 2) +
#   labs(x = "x", y = "y") +
#   theme_minimal(base_size = 25)
#
# animal_stats <- animalxy.all %>%
#   dplyr::group_by(Animal_ID) %>%
#   dplyr::summarise(
#     tot_dist = sum(trav_dist, na.rm = T)
#   )
#
Data_summary <- all_data %>%
  group_by(iteration) %>%
  dplyr::summarise(
    num_counts = sum(count_data[[1]]$count),
    num_det = cam_captures[[1]] %>%
      dplyr::filter(in_cam == T) %>%
      dplyr::group_by(lscape_index, t) %>%
      dplyr::summarise(
        group = n(),
        .groups = "drop"
      ) %>%
      dplyr::summarise(num_det = sum(group > 0)) %>%
      pull(num_det),
    det_prob = num_det / (study_design$t_steps * cam_design$ncam),
    num_encounters = sum(encounter_data[[1]]$encounter),
    mean_stay_time = mean(as.matrix(stay_time_data[[1]]), na.rm = T),
    max_stay_time = max(stay_time_data[[1]], na.rm = T),
    num_stay_time_outliers_1 = sum(stay_time_data[[1]] > study_design$dt, na.rm = T),
    num_stay_time_outliers_2 = sum(stay_time_data[[1]] > 2 * study_design$dt, na.rm = T),
    num_stay_time_outliers_3 = sum(stay_time_data[[1]] > 3 * study_design$dt, na.rm = T)
  )


# results_fast_cam_alt <- list(
#   # save_animal_data,
#   study_design,
#   cam_design,
#   lscape_design,
#   all_data,
#   D_all
# )
#
# save(results_fast_cam_alt, file = "Sim_results/results_fast_cam_alt.RData")


# save(save_animal_data, file = "Sim_results/save_animal_data.RData")
# save(save_lscape_defs, file = "Sim_results/save_lscape_defs.RData")
