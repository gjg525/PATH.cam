library(tidyverse)
library(RColorBrewer)
library(lattice)
library(gridExtra)
library(doParallel)
devtools::load_all()

sim_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/sim_results/"
sim_dir_REST <- "G:/My Drive/Missoula_postdoc/PATH_model/sim_results_REST/"

# Initializations
fig_colors <- c("#2ca25f", "#fc8d59", "#67a9cf", "#f768a1", "#bae4b3", "#fed98e")
options(ggplot2.discrete.colour = fig_colors)
options(ggplot2.discrete.fill = fig_colors)

# Run with different number of cameras
cam_tests <- c(25, 50, 75, 100, 125)
# cam_tests <- c(50, 75, 100)
# cam_tests <- c(100)

# Load animal GPS data
# load(file = paste0(sim_dir, "save_animal_data_1.RData"))
load(file = paste0(sim_dir, "save_lscape_defs.RData"))

# save_animal_data <- save_animal_data_1
# rm(save_animal_data_1)

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

all_designs <- tibble::tibble(
  Design_name = c("Random", "Slow_80_bias", "Medium_80_bias", "Fast_80_bias",
                  "Slow_bias", "Medium_bias", "Fast_bias"),
  Design =c("Random", "Bias", "Bias", "Bias", "Bias", "Bias", "Bias"),
  Props = c(
    list(c(1, 1, 1)),
    list(c(.8, .1, .1)),
    list(c(.1, .8, .1)),
    list(c(.1, .1, .8)),
    list(c(1, 0, 0)),
    list(c(0, 1, 0)),
    list(c(0, 0, 1))
  )
)
# # Design for single run with 250 cameras
cam_tests <- 250
all_designs <- tibble::tibble(
  Design_name = c("Random", "Slow_80_bias", "Med_80_bias", "Fast_80_bias",
                  "Slow_bias", "Med_bias", "Fast_bias"),
  Design =c("Random", "Bias", "Bias", "Bias", "Bias", "Bias", "Bias"),
  Props = c(
    list(c(1, 1, 1)),
    list(c(.8, .1, .1)),
    list(c(.1, .8, .1)),
    list(c(.1, .1, .8)),
    list(c(1, 0, 0)),
    list(c(0, 1, 0)),
    list(c(0, 0, 1))
  )
)
for (cam_des in 1:nrow(all_designs)) {
  for (cam in 1:length(cam_tests)) {
    
    # Cam designs
    cam_design <- tibble::tibble(
      ncam = cam_tests[cam],
      Design_name = all_designs$Design_name[cam_des],
      Design = all_designs$Design[cam_des],
      Props = all_designs$Props[cam_des],
      cam_length = study_design$dx * 0.1, # length of all viewshed sides
      cam_A = cam_length ^ 2 / 2,
      tot_snaps = ncam * study_design$t_steps
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
      
      # # # Create custom landscape (tag = "grid", "circ", "squares", "metapop")
      # lscape_defs <- lscape_creator(study_design, lscape_design)
      # 
      # # # Run agent-based model
      # animalxy.all <- ABM_sim(study_design,
      #                         lscape_defs)
      # # If running new ABM simulation on each run, save
      # save_animal_data$data[run] <- list(animalxy.all)
      # 
      # # If running new ABM simulation on each run, save
      # save_lscape_defs$data[run] <- list(lscape_defs)
      
      if (run %in% seq(1, 901, by = 100)) {
        dat_name <- paste0("save_animal_data_", ceiling(run / 100))
        load(file = paste0(sim_dir, dat_name, ".RData"))
        
        save_animal_data <- get(dat_name)
        rm(list = paste0("save_animal_data_", ceiling(run / 100)))
      }
      
      # Load ABM from save file
      animalxy.all <- save_animal_data$data[[(run - 1) %% 100 + 1]] %>% 
        dplyr::rename(Road = road)
      lscape_defs <- save_lscape_defs$data[[run]]
      
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
      
      all_data$count_data[[run]] <- list(count_data)
      
      encounter_data <- get_encounter_data(cam_locs, all_data$cam_captures[[run]])$encounter
      stay_time_data <- get_stay_time_data(cam_locs, all_data$cam_captures[[run]])[[2]] |>
        as.matrix()
      
      # Run models only if any data points were collected
      if (sum(count_data$count) == 0) {
        D.PR.MCMC.habitat <- NA
        SD.PR.MCMC.habitat <- NA
      } else {
        # chain.PATH <- fit.model.mcmc.PATH(
        #   study_design = study_design,
        #   cam_design = cam_design,
        #   cam_locs = cam_locs,
        #   gamma_start = rep(log(mean(count_data$count)), study_design$num_covariates),
        #   gamma_prior_var = 10^4,
        #   gamma_tune = rep(-1, study_design$num_covariates),
        #   kappa_start = log(exp(kappa.prior.mu) / sum(exp(kappa.prior.mu))),
        #   kappa_prior_mu = kappa.prior.mu,
        #   kappa_prior_var = kappa.prior.var,
        #   kappa_tune = -1, #rep(-1, study_design$num_covariates),
        #   count_data_in = count_data,
        #   habitat_summary
        # )
        # 
        # ## Posterior summaries
        # # plot(chain.PATH$tot_u[study_design$burn_in:study_design$n_iter])
        # D.PATH.MCMC <- mean(chain.PATH$tot_u[study_design$burn_in:study_design$n_iter])
        # SD.PATH.MCMC <- sd(chain.PATH$tot_u[study_design$burn_in:study_design$n_iter])
        # 
        # if (any(colMeans(chain.PATH$accept[study_design$burn_in:study_design$n_iter, ]) < 0.2) || any(colMeans(chain.PATH$accept[study_design$burn_in:study_design$n_iter, ]) > 0.7)) {
        #   warning(("Mean Count accept rate OOB"))
        #   D.PATH.MCMC <- NA
        #   SD.PATH.MCMC <- NA
        # }
        
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
      
      ncam_temp <- cam_design$ncam
      
      D_all[[(run - 1) * 2 + 1]] <- tibble::tibble(
        iteration = run,
        Model = "PATH",
        Covariate = "Non-Covariate",
        Est = NA, #D.PATH.MCMC,
        SD = NA, #SD.PATH.MCMC
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
        cam_design = cam_design$Design_name,
        cams = ncam_temp,
        Model = "IS",
        Est = IS_mean,
        SD = SE_N
        # all_results = list(chain.PR.habitat)
      )
      
    }
    
    # save_results <- list(
    #   # save_animal_data,
    #   study_design,
    #   cam_design,
    #   lscape_design,
    #   all_data,
    #   D_all
    # )
    # 
    # save(save_results, file = paste0(sim_dir,
    #                                  cam_design$Design_name,
    #                                  "_",
    #                                  cam_design$ncam,
    #                                  "_cam.RData")
    # )
    
    # rm(save_results, all_data, D_all)
    rm(all_data, D_all)
    
    save_results_REST <- list(
      study_design,
      cam_design,
      lscape_design,
      D_all_REST
    )
    
    save(save_results_REST, file = paste0(sim_dir_REST,
                                     cam_design$Design_name,
                                     "_",
                                     cam_design$ncam,
                                     "_cam_REST.RData")
    )
    
    rm(save_results_REST, D_all_REST)
    
  }
}


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
# plot_space_use(study_design,
#                animalxy.all)
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
