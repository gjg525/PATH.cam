library(tidyverse)
library(RColorBrewer)
library(lattice)
library(gridExtra)
library(doParallel)
# library(PATH.cam)
devtools::load_all()

# Initializations
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")
options(ggplot2.discrete.colour = fig_colors)
options(ggplot2.discrete.fill = fig_colors)

# Study design
study_design <- tibble::tibble(
  q = 30^2, # Number grid cells
  dx = 1,  # Grid cell lengths
  dy = 1,
  t_steps = 500, # Number of time steps
  dt = 1, # Time step size
  t_censor = 2,
  bounds = list(c(0, dx * q ^ 0.5)), # Sampling area boundaries
  tot_A = (bounds[[1]][2] - bounds[[1]][1])^2,
  num_groups = 100,
  group_sizes = list(rep(1, num_groups)),
  group_spread = 0, # Tightness of grouping behavior (relative to grid size)
  tot_animals = sum(unlist(group_sizes)),
  # MCMC parms
  num_runs = 1000,
  n_iter = 40000,
  burn_in = 30000,
  covariate_labels = list(c("Slow", "Medium", "Fast"))
)

# Landscape design
lscape_design <- tibble::tibble(
  lscape_tag = "Random", # "Custom", #
  Speed_ID = c("Slow", "Medium", "Fast"),
  Speed_mins = c(.19, .4, .9),
  Speed_maxes = c(.21, .5, 1.1),
  Probs = c(0.1, 0.1, 0.8)
)

# Cam designs
cam_design <- tibble::tibble(
  ncam = 250,
  Design = "Random",
  Props = list(c(1, 1, 1)), # proportion of cameras placed in each habitat type
  # Design = "Bias",
  # Props = list(c(0.8, 0.1, 0.1)),
  # Props = list(c(0.1, 0.8, 0.1)),
  # Props = list(c(0.1, 0.1, 0.8)),
  # Props = list(c(1, 0, 0)),
  # Props = list(c(0, 1, 0)),
  # Props = list(c(0, 0, 1)),
  # Props = list(c(.3, 0.1, .6)),
  cam_length = study_design$dx * 0.1, # length of all viewshed sides
  cam_A = cam_length ^ 2 / 2,
  tot_snaps = ncam * study_design$t_steps
)

# Initialize summary matrices
D_all <- tibble::tibble(
  iteration = rep(1:study_design$num_runs, each = 2),
  Model = NA,
  Covariate = NA,
  Est = NA,
  SD = NA
)

all_data <- tibble::tibble(
  iteration = 1:study_design$num_runs,
  cam_captures = NA,
  count_data = NA
)

# # camera data summaries
cam_design <- cam_design %>%
  dplyr::mutate(percent_cam_coverage = ncam * cam_A / study_design$tot_A)

# Multi-run simulations
D_all <- vector(mode = "list", length = study_design$num_runs * 2)
D_all_REST <- vector(mode = "list", length = study_design$num_runs * 2)
for (run in 1:study_design$num_runs) {
  print(paste("Run", run, "of", study_design$num_runs))

  # Create custom landscape (tag = "grid", "circ", "squares", "metapop")
  lscape_defs <- lscape_creator(study_design, lscape_design)

  # Run agent-based model
  animalxy.all <- ABM_sim(study_design,
                          lscape_defs)

  # Create covariate matrix with 0, 1 values
  study_design <- study_design %>%
    dplyr::mutate(
      num_covariates = length(unlist(covariate_labels)),
      Z = list(create_covariate_mat(
        lscape_defs,
        study_design,
        unlist(covariate_labels)))
    )

  tele_summary <- Collect_tele_data(animalxy.all, study_design)

  # Use smallest stay time as reference category
  ref_cat_idx <- which(tele_summary$stay_prop == min(tele_summary$stay_prop))

  # Set reference category for intercept
  study_design$Z[[1]][, ref_cat_idx] <- 1

  # Subtract reference category from stay time proportion
  prop_adjust <- tele_summary$stay_prop /
    tele_summary$stay_prop[ref_cat_idx]
  prop_adjust[ref_cat_idx] <- tele_summary$stay_prop[ref_cat_idx]
  kappa.prior.mu.adj <- log(prop_adjust)
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
    dplyr::group_by(Speed) %>%
    dplyr::summarise(
      n_lscape = dplyr::n()
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
    dplyr::select(Speed, n_lscape, prop_cams, d_coeff)

  habitat_summary$Speed <- factor(
    habitat_summary$Speed,
    levels = unlist(study_design$covariate_labels)
  )

  count_data <- get_count_data(
    cam_locs,
    all_data$cam_captures[[run]],
    animalxy.all %>%
      dplyr::filter(t != 0))

  encounter_data <- get_encounter_data(cam_locs, all_data$cam_captures[[run]])$encounter
  stay_time_data <- get_stay_time_data(cam_locs, all_data$cam_captures[[run]])[[2]] |>
    as.matrix()
  ################################################################################
  if (sum(count_data$count) == 0) {
    D.PATH.MCMC <- NA
    SD.PATH.MCMC <- NA
    D.REST.MCMC <- NA
    SD.REST.MCMC <- NA
    D.REST.MCMC.cov <- NA
    SD.REST.MCMC.cov <- NA

  } else {
    chain.PATH <- fit.model.mcmc.PATH(
      study_design = study_design,
      cam_design = cam_design,
      cam_locs = cam_locs,
      gamma_start = runif(3, -10, 10),# rep(log(mean(count_data$count)), study_design$num_covariates),
      gamma_prior_var = 10^4,
      gamma_tune = rep(-1, study_design$num_covariates),
      kappa_start = runif(1, -10, 10), #log(exp(kappa.prior.mu) / sum(exp(kappa.prior.mu))),
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

    ################################################################################
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
      gamma_start = rep(log(mean(encounter_data)), 3),
      kappa_start = rep(log(mean(stay_time_data,na.rm=T)), 3),
      gamma_prior_var = 10^4,
      kappa_prior_var = 10^4,
      gamma_tune = c(-1, -1, -1),
      kappa_tune = c(-1, -1, -1),
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
    Model = "PATH",
    Covariate = "Non-Covariate",
    Est = D.PATH.MCMC,
    SD = SD.PATH.MCMC
    # all_results = list(chain.PATH)
  )

  D_all_REST[[(run - 1) * 2 + 1]] <- tibble::tibble(
    iteration = run,
    Model = "REST",
    Covariate = "Non-Covariate",
    Est = D.REST.MCMC,
    SD = SD.REST.MCMC
  )

  D_all_REST[[run * 2]] <- tibble::tibble(
      iteration = run,
      Model = "REST",
      Covariate = "Covariate",
      Est = D.REST.MCMC.cov,
      SD = SD.REST.MCMC.cov
    )
  ################################################################################
  # IS method
  IS_mean <- sum(count_data$count) / study_design$t_steps / cam_design$ncam /
    cam_design$cam_A * study_design$tot_A

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
    Model = "IS",
    Covariate = "Non-Covariate",
    Est = IS_mean,
    SD = SE_N
    # all_results = list(chain.PATH)
  )

}

D_all <- dplyr::bind_rows(D_all)

# Remove outlier estimates
D_all$Est[D_all$Est > 5 * study_design$tot_animals] <- NA
D_all$SD[D_all$SD > 5 * study_design$tot_animals] <- NA

D_all_REST <- dplyr::bind_rows(D_all_REST)

# Remove outlier estimates
D_all_REST$Est[D_all_REST$Est > 5 * study_design$tot_animals] <- NA
D_all_REST$SD[D_all_REST$SD > 5 * study_design$tot_animals] <- NA

###########################
# Plots
###########################
plot_multirun_means(study_design, D_all_REST %>%
                      dplyr::filter(is.finite(Est)))
plot_multirun_sds(D_all_REST %>%
                    dplyr::filter(is.finite(Est)))

plot_ABM(study_design,
         cam_design,
         cam_locs,
         animalxy.all)
plot_space_use(study_design,
               animalxy.all)

