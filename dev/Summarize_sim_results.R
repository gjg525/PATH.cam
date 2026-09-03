library(dplyr)
library(ggplot2)
devtools::load_all()
tot_animals <- 100

save_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/D_all_results"
fig_dir <- "G:/My Drive/Missoula_postdoc/PATH_model/imgs/"

fig_colors <- c("#1B5E20", "#00A8C6", "#E65100", "#FBC02D", "#8E44AD", "#4B6FAD", "#D81B60")

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
D_dat <- tibble::tibble()
for (ii in 1:length(file_names)) {
  print(ii)

  file_ii <- paste0(sim_dir, design_names[ii], "_", ncams, "_cam_10min.RData")
  # file_ii <- paste0(sim_dir, design_names[ii], "_", ncams, "_cam.RData")
  if (file.exists(file_ii)) {
    all_results <- loadRData(file_ii)

    D_all <- D_all |>
      dplyr::bind_rows(
        all_results[[4]] |>
          dplyr::bind_rows() |>
          dplyr::mutate(
            SampDesign = paste0(file_names[ii], "_cam")
          )
      )

    D_dat <- dplyr::bind_rows(
      D_dat,
      all_results[[4]] %>%
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


################################################################################
all_count_data <- D_dat |>
  dplyr::select(iteration, count_data, SampDesign) |>
  tidyr::unnest(count_data) |>
  tidyr::unnest(count_data)

all_tele_data <- D_dat |>
  dplyr::select(iteration, tele_summary, SampDesign) |>
  tidyr::unnest(tele_summary) |>
  tidyr::unnest(tele_summary) |>
  dplyr::left_join(
    D_dat |>
      dplyr::select(iteration, tele_summary_full, SampDesign) |>
      tidyr::unnest(tele_summary_full) |>
      tidyr::unnest(tele_summary_full) |>
      dplyr::select(res_idx_truth = stay_prop, iteration, Speed, SampDesign),
    by = join_by(iteration, Speed, SampDesign)
  ) |>
  dplyr::mutate(
    res_idx_error = stay_prop - res_idx_truth
  )

all_encounter_data <- D_dat |>
  dplyr::select(iteration, encounter_data, SampDesign) |>
  tidyr::unnest(encounter_data) |>
  tidyr::unnest(encounter_data)

all_stay_data <- D_dat |>
  dplyr::select(iteration, stay_time_data, SampDesign) |>
  tidyr::unnest(stay_time_data) |>
  mutate(
    stay_time_data = purrr::map(stay_time_data, function(mat) {
      n_rows <- nrow(mat)
      if (is.null(n_rows)) {
        return(list(mat))
      }
      purrr::map(seq_len(n_rows), ~ mat[.x, , drop = FALSE])
    })
  ) %>%
  tidyr::unnest_longer(stay_time_data)

all_data <- all_count_data |>
  dplyr::left_join(
    all_tele_data,
    by = dplyr::join_by(iteration, Speed, SampDesign)
  ) |>
  dplyr::bind_cols(
    all_encounter_data |>
      dplyr::select(encounter_data)
  ) |>
  dplyr::bind_cols(
    all_stay_data |>
      dplyr::select(stay_time_data)
  )

all_data |>
  dplyr::filter(count > 0) |>
  ggplot2::ggplot(ggplot2::aes(x = count, fill = Speed)) +
  ggplot2::geom_histogram()

all_data |>
  dplyr::filter(encounter_data > 0) |>
  ggplot2::ggplot(ggplot2::aes(x = encounter_data, fill = Speed)) +
  ggplot2::geom_histogram()
  # ggplot2::facet_grid(~ ABM_Design, scales = "free_x")

all_data %>%
  dplyr::mutate(
    stay_time_data = purrr::map(stay_time_data, as.numeric)
  ) %>%
  tidyr::unnest_longer(stay_time_data) %>%
  dplyr::filter(!is.na(stay_time_data)) |>
  # dplyr::filter(stay_time_data > 0) |>
  ggplot2::ggplot(ggplot2::aes(x = stay_time_data, fill = Speed)) +
  ggplot2::geom_histogram()

# count_summary <- all_count_data |>
#   dplyr::group_by(iteration, Speed, SampDesign) |>
#   dplyr::summarise(
#     count = sum(count),
#     .groups = 'drop'
#   )
#

# tele_summary_data <- tele_data |>
#   dplyr::group_by(Speed) |>
#   dplyr::summarise(
#     mean_res_idx = mean(stay_prop),
#     mean_res_idx_error = mean(res_idx_error),
#     sd_res_idx = sd(stay_prop),
#     .groups = 'drop'
#   )

# habitat_D <- landscape_summary |>
#   dplyr::left_join(
#     tele_summary_data,
#     by = dplyr::join_by(Speed)
#   ) |>
#   dplyr::mutate(
#     N_habitat = tot_animals * mean_res_idx, #average number of animals in habitat
#     D = N_habitat / habitat_area
#   )


all_data_est <- all_habitat_data |>
  dplyr::left_join(
    D_all_NLCD |>
      dplyr::filter(Model == "PATH" & Covariate == "Non-Covariate"),
    by = dplyr::join_by(iteration, SampDesign)
  ) |>
  dplyr::left_join(
    count_summary,
    by = dplyr::join_by(iteration, Speed, SampDesign)
  ) |>
  dplyr::left_join(
    tele_data,
    by = dplyr::join_by(iteration, Speed, SampDesign)
  )


all_data_est_summary <- all_data_est |>
  dplyr::group_by(iteration, SampDesign, Model, Covariate) |>
  dplyr::summarise(
    Est = unique(Est),
    SD = unique(SD),
    count = sum(count, na.rm = T),
    .groups = 'drop'
  )

all_data_est |>
  ggplot2::ggplot(ggplot2::aes(x = log(abs(res_idx_error)), y = Est - tot_animals, color = Speed)) +
  ggplot2::geom_point()

all_data_est |>
  dplyr::mutate(count_ratio = count / d_coeff) |>
  ggplot2::ggplot(ggplot2::aes(
    x = count,
    y = Est,
    color = Speed
  )) +
  ggplot2::geom_point()

all_data_est_summary |>
  ggplot2::ggplot(ggplot2::aes(
    x = count,
    y = Est,
    color = SampDesign
  )) +
  ggplot2::geom_point()

