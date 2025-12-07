#' Fit Landscape Density Model using MCMC
#'
#' @description
#' Implements a custom Adaptive Metropolis-Hastings MCMC algorithm to estimate
#' landscape-dependent density parameters (`gamma`) and staying time proportions
#' (`kappa`). The model integrates camera trap counts with landscape covariates
#' to estimate total abundance (`tot_u`).
#'
#' @details
#' The function estimates parameters using a Poisson likelihood for camera counts
#' and Normal priors for coefficients. It employs an adaptive tuning mechanism
#' where proposal variances (`gamma_tune`, `kappa_tune`) are adjusted every 100
#' iterations to target an acceptance rate of approximately 44%.
#'
#' The model assumes:
#' \itemize{
#'   \item \code{d = exp(gamma)}: Density intensity per landscape type.
#'   \item \code{stay_prop = exp(kappa)}: Proportional staying time.
#' }
#'
#' @param study_design A list of simulation/study parameters. Must include:
#' \itemize{
#'   \item \code{n_iter}: Number of MCMC iterations.
#'   \item \code{dx}, \code{dy}: Grid cell dimensions.
#'   \item \code{t_steps}: Number of time steps.
#'   \item \code{covariate_labels}: Names of the covariates (e.g., speed categories).
#'   \item \code{num_covariates}: Integer count of covariates.
#'   \item \code{Z}: Covariate design matrix.
#'   \item \code{q}: Total number of grid cells.
#' }
#' @param cam_design A list containing camera design details, specifically
#'   \code{cam_A} (area of camera view).
#' @param cam_locs A data frame of camera locations (passed but not currently used
#'   in the calculation logic provided).
#' @param gamma_start Vector of initial values for density coefficients.
#' @param gamma_prior_var Numeric, variance for the Normal prior on gamma.
#' @param gamma_tune Vector of initial tuning variances for gamma proposals.
#' @param kappa_start Vector of initial values for staying time coefficients.
#' @param kappa_prior_mu Vector of means for the Normal prior on kappa.
#' @param kappa_prior_var Numeric, variance for the Normal prior on kappa.
#' @param kappa_tune Vector of initial tuning variances for kappa proposals.
#' @param count_data_in A data frame containing observed counts.
#'   \strong{Crucial:} If \code{grouping = "Speed"}, this must contain a column named "Speed".
#'   If \code{grouping} is anything else, it must contain a column named "Road".
#' @param habitat_summary A list/data frame containing landscape summaries:
#' \itemize{
#'   \item \code{n_lscape}: Count of cells per landscape type.
#'   \item \code{prop_cams}: Proportion of cameras in each landscape type.
#' }
#' @param grouping Character string (default "Speed"). determines which column
#'   in \code{count_data_in} is used to subset counts.
#'
#' @return A list containing MCMC chains and diagnostics:
#' \itemize{
#'   \item \code{accept}: Matrix of acceptance (1/0) events for diagnostics.
#'   \item \code{gamma}: Matrix of posterior samples for density coefficients.
#'   \item \code{kappa}: Matrix of posterior samples for staying time coefficients.
#'   \item \code{tot_u}: Vector of estimated total abundance/density across iterations.
#' }
#'
#' @export
#'
fit.model.mcmc.PATH <- function(study_design,
                                      cam_design,
                                      cam_locs,
                                      gamma_start,
                                      gamma_prior_var = 10^4,
                                      gamma_tune,
                                      kappa_start,
                                      kappa_prior_mu,
                                      kappa_prior_var,
                                      kappa_tune,
                                      count_data_in,
                                      habitat_summary,
                                      grouping = "Speed") {

  n_iter <- study_design$n_iter
  cam_A <- cam_design$cam_A
  cell_A <- study_design$dx * study_design$dy
  t_steps <- study_design$t_steps
  covariate_labels <- unlist(study_design$covariate_labels)
  num_covariates <- study_design$num_covariates
  Z <- matrix(unlist(study_design$Z), study_design$q, study_design$num_covariates)

  # Variables that will be saved
  gamma <- matrix(, n_iter + 1, num_covariates)
  kappa <- matrix(, n_iter + 1, num_covariates)
  tot_u <- matrix(, n_iter + 1, 1)
  # accept <- matrix(, n_iter + 1, 2 * num_covariates)
  accept <- matrix(, n_iter + 1, num_covariates + 1)
  gamma[1, ] <- gamma_start
  colnames(gamma) <- paste0("gamma.", covariate_labels)
  kappa[1, ] <- kappa_start
  colnames(kappa) <- paste0("kappa.", covariate_labels)
  colnames(tot_u) <- "Total estimate"
  # colnames(accept) <- c(
  #   paste0("accept.rate.gamma.", covariate_labels),
  #   paste0("accept.rate.kappa.", covariate_labels)
  # )
  colnames(accept) <- c(
    paste0("accept.rate.gamma.", covariate_labels),
    "accept.rate.kappa."
  )

  d <- exp(gamma[1, ])
  stay_prop <- exp(kappa[1, ])

  u <- d * habitat_summary$n_lscape / stay_prop *
    habitat_summary$prop_cams / (cam_design$cam_A * study_design$t_steps)

  tot_u <- sum(u)

  tune_check <- 100
  batch_n <- 0

  # Begin MCMC loop
  # prog_bar <- txtProgressBar(min = 0,max = n_iter,style = 3,width = 50,char = "=")
  for (i in 1:n_iter) {
    # setTxtProgressBar(prog_bar, i)
    # Sample gamma
    for (gg in 1:num_covariates) {
      gamma_star <- gamma[i, ]
      gamma_star[gg] <- rnorm(1, gamma[i, gg], exp(2 * gamma_tune[gg]))
      d_star <- exp(gamma_star)

      if (grouping == "Speed") {
        count_data_in_h <-
          count_data_in$count[count_data_in$Speed == covariate_labels[gg]]
      } else {
        count_data_in_h <-
          count_data_in$count[count_data_in$Road == covariate_labels[gg]]
      }

      if (length(count_data_in_h) == 0) {
        count_data_in_h <- 0
      }

      if (all(d_star > 0)) {
        mh1 <- sum(dpois(count_data_in_h, d_star[gg], log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma_star[gg], 0, gamma_prior_var^0.5, log = TRUE))
        mh2 <- sum(dpois(count_data_in_h, d[gg], log = TRUE), na.rm = TRUE) +
          sum(dnorm(gamma[i, gg], 0, gamma_prior_var^0.5, log = TRUE))
        mh <- exp(mh1 - mh2)

        if (mh > runif(1)) {
          gamma[i, ] <- gamma_star
          accept[i + 1, gg] <- 1
        } else {
          gamma[i, ] <- gamma[i, ]
          accept[i + 1, gg] <- 0
        }
      } else {
        gamma[i, ] <- gamma[i, ]
        accept[i + 1, gg] <- 0
      }

    }
    gamma[i + 1, ] <- gamma[i, ]
    d <- exp(gamma[i + 1, ])

    # Sample kappa
    kappa_star <- kappa[i, ]
    kappa_star <- rnorm(num_covariates, kappa[i, ], exp(2 * kappa_tune))
    # kappa_star <- truncnorm::rtruncnorm(num_covariates, 0, Inf, kappa[i, ], exp(2 * kappa_tune))
    kappa_star <- log(exp(kappa_star) / sum(exp(kappa_star)))

    # Proportional staying time defined by priors
    mh1 <- sum(dnorm(kappa_star,kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    mh2 <- sum(dnorm(kappa[i,],kappa_prior_mu,kappa_prior_var^0.5,log=TRUE))
    mh <- exp(mh1-mh2)

    if (mh > runif(1) & !is.na(mh)) {
      kappa[i, ] <- kappa_star
      accept[i + 1, num_covariates + 1] <- 1
      # u <- u_star
    } else {
      kappa[i, ] <- kappa[i, ]
      accept[i + 1, num_covariates + 1] <- 0
    }

    kappa[i + 1, ] <- kappa[i, ]

    stay_prop <- exp(kappa[i + 1, ])

    u <- d * habitat_summary$n_lscape / stay_prop *
      habitat_summary$prop_cams / (cam_design$cam_A * study_design$t_steps)

    tot_u[i + 1] <- sum(u)

    # Update tuning parms
    if (i %% tune_check == 0) {
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1
      # delta_n <- min(0.01,batch_n^-1)
      accept_gamma_check <- colMeans(accept[(i - tune_check + 1):i, 1:num_covariates], na.rm = T)
      #mean(accept[(i - tune_check + 1):i], na.rm = T)
      gamma_tune[accept_gamma_check > 0.44] <- gamma_tune[accept_gamma_check > 0.44] + delta_n
      gamma_tune[accept_gamma_check <= 0.44] <- gamma_tune[accept_gamma_check <= 0.44] - delta_n

      accept_kappa_check <- mean(accept[(i - tune_check + 1):i, num_covariates + 1], na.rm = T)
      kappa_tune[accept_kappa_check > 0.44] <- kappa_tune[accept_kappa_check > 0.44] + delta_n
      kappa_tune[accept_kappa_check <= 0.44] <- kappa_tune[accept_kappa_check <= 0.44] - delta_n
    }
  }
  # print("MCMC complete")

  list(accept = accept, gamma = gamma, kappa = kappa, tot_u = tot_u)
}
