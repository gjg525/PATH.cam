library(dplyr)
library(ggplot2)
devtools::load_all()
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
# (This makes it easier to analyze)
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

################################################################################
# Download NLCD map
tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5010NLCDclip.tif"

mu_base <- tibble::tibble(
  LandCover = c("Water", "Development", "Forest", "Agriculture"),
  Mu = c(4, 2, 0.02, 0.5),
  speed = c(0, 1.5, 0.7, 0.3) * 30, # Manually set km/hr and convert to cell/hr
  speed_km_hr = 4 * Mu * 900 / 30 / 1000
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

landscape_summary <- df |>
  dplyr::group_by(Speed = LandCover) |>
  dplyr::summarise(
    prop_habitat = dplyr::n() / 8100,
    habitat_area = dplyr::n() * 0.0009, # Assumes 30mx30m grid cell
    .groups = 'drop'
  )

all_habitat_data <- NLCD_dat |>
  dplyr::select(iteration, habitat_summary, SampDesign, ABM_Design) |>
  tidyr::unnest(habitat_summary) |>
  tidyr::unnest(habitat_summary)


################################################################################
all_count_data <- NLCD_dat |>
  dplyr::select(iteration, count_data, SampDesign, ABM_Design) |>
  tidyr::unnest(count_data) |>
  tidyr::unnest(count_data)

all_tele_data <- NLCD_dat |>
  dplyr::select(iteration, tele_summary, SampDesign, ABM_Design) |>
  tidyr::unnest(tele_summary) |>
  tidyr::unnest(tele_summary) |>
  dplyr::left_join(
    NLCD_dat |>
      dplyr::select(iteration, tele_summary_full, SampDesign, ABM_Design) |>
      tidyr::unnest(tele_summary_full) |>
      tidyr::unnest(tele_summary_full) |>
      dplyr::select(res_idx_truth = stay_prop, iteration, Speed, SampDesign, ABM_Design),
    by = join_by(iteration, Speed, SampDesign, ABM_Design)
  ) |>
  dplyr::mutate(
    res_idx_error = stay_prop - res_idx_truth
  )

all_encounter_data <- NLCD_dat |>
  dplyr::select(iteration, encounter_data, SampDesign, ABM_Design) |>
  tidyr::unnest(encounter_data) |>
  tidyr::unnest(encounter_data)

all_stay_data <- NLCD_dat |>
  dplyr::select(iteration, stay_time_data, SampDesign, ABM_Design) |>
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
    by = dplyr::join_by(iteration, Speed, SampDesign, ABM_Design)
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
  ggplot2::geom_histogram() +
  ggplot2::facet_grid(~ ABM_Design, scales = "free_x")

all_data |>
  dplyr::filter(encounter_data > 0) |>
  ggplot2::ggplot(ggplot2::aes(x = encounter_data, fill = Speed)) +
  ggplot2::geom_histogram() +
  ggplot2::facet_grid(~ ABM_Design, scales = "free_x")

all_data %>%
  dplyr::mutate(
    stay_time_data = purrr::map(stay_time_data, as.numeric)
  ) %>%
  tidyr::unnest_longer(stay_time_data) %>%
  dplyr::filter(!is.na(stay_time_data)) |>
  # dplyr::filter(stay_time_data > 0) |>
  ggplot2::ggplot(ggplot2::aes(x = stay_time_data, fill = Speed)) +
  ggplot2::geom_histogram() +
  ggplot2::facet_grid(~ ABM_Design, scales = "free_x")

# count_summary <- all_count_data |>
#   dplyr::group_by(iteration, Speed, SampDesign, ABM_Design) |>
#   dplyr::summarise(
#     count = sum(count),
#     .groups = 'drop'
#   )
#

# tele_summary_data <- tele_data |>
#   dplyr::group_by(Speed, ABM_Design) |>
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
    by = dplyr::join_by(iteration, SampDesign, ABM_Design)
  ) |>
  dplyr::left_join(
    count_summary,
    by = dplyr::join_by(iteration, Speed, SampDesign, ABM_Design)
  ) |>
  dplyr::left_join(
    tele_data,
    by = dplyr::join_by(iteration, Speed, SampDesign, ABM_Design)
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
  ggplot2::geom_point() +
  ggplot2::facet_grid(~ ABM_Design)

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

