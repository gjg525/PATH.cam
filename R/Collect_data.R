#' Extract Camera Trap Captures from Trajectories
#'
#' @description
#' Identifies moments when simulated animals enter the field of view of simulated
#' cameras. It filters trajectories by grid cell, calculates precise entry/exit
#' times, and computes the duration of the stay (`t_stay`).
#'
#' @param animalxy A data frame of animal trajectories (output from \code{\link{ABM_sim}}).
#' @param cam_locs A data frame containing camera locations. Must include columns:
#'   \code{lscape_index}, \code{X}, \code{Y}, and \code{cam_ID}.
#' @param study_design A list of simulation parameters (see \code{\link{ABM_sim}}).
#'
#' @return A tibble containing "capture events" with columns:
#' \itemize{
#'   \item \code{pass_i}: Unique identifier for a single pass event.
#'   \item \code{t_stay}: Duration of time spent in the camera view.
#'   \item \code{in_cam}: Boolean, indicating if the animal was inside the camera view.
#'   \item \code{encounter}: The specific encounter index.
#'   \item \code{t_in}, \code{t_out}: Entry and exit timestamps.
#' }
#' @export
#'
get_cam_captures <- function(animalxy, cam_locs, study_design) {
  t_steps <- study_design$t_steps
  dt <- study_design$dt

  # collect data within each grid cell to reduce subsequent calculations
  cell_captures <- animalxy |>
    dplyr::ungroup() |>
    dplyr::filter(lscape_index %in% cam_locs$lscape_index) |>
    dplyr::mutate(pass_i = cumsum(c(1, (abs(ii[-n()] - ii[-1]) > 1) |
      (abs(lscape_index[-n()] - lscape_index[-1]) > 0)))) |>
    dplyr::group_by(pass_i) |>
    dplyr::mutate(next_i = ifelse(max(t) == t_steps * dt,
      max(ii),
      max(ii) + 1
    )) |>
    dplyr::ungroup()

  # # Use animal's next step to capture all trajectories within cell
  cell_captures <- cell_captures |>
    dplyr::add_row(
      animalxy |>
        dplyr::ungroup() |>
        dplyr::filter(ii %in% unique(cell_captures$next_i)) |>
        dplyr::mutate(pass_i = unique(cell_captures$pass_i))
      ) |>
    dplyr::arrange(ii) |>
    dplyr::select(-next_i)

  # Remove repeat rows (occurs for max time step)
  cell_captures <- dplyr::distinct(cell_captures)

  cell_check <- cell_captures |>
    dplyr::group_by(pass_i) |>
    dplyr::mutate(
      num_animals = length(unique(Animal_ID)),
      num_lscapes = length(unique(lscape_index))
    )

  # Check that only 1 animal is included in each pass
  if (any(cell_check$num_animals > 1)) {
    stop("More than 1 animal in a pass")
  }
  # Check that only 2 cells are included in each pass
  if (any(cell_check$num_lscapes > 2)) {
    stop("More than 2 cells in a pass")
  }

  if (nrow(cell_captures) > 0) {
    # Collect all camera viewshed captures
    # in_cam determines whether an individual starts within the camera viewshed
    cam_captures <- cell_captures |>
      dplyr::group_by(pass_i) |>
      # dplyr::mutate(pass_i = pass_i[t == min(t)]) |>
      dplyr::summarise(
        xy_index = list(cam_locs[cam_locs$lscape_index %in% lscape_index[1], 3:4]),
        cam_intersects = list(
          calc_intersects(
            matrix(unlist(xy_index),
              nrow = length(unlist(xy_index)) / 2,
              ncol = 2
            ),
            cbind(X, Y),
            trav_speeds,
            t
          )
        ),
        t_stay = c(unlist(lapply(cam_intersects, "[[", 2)), NA),
        in_cam = c(unlist(lapply(cam_intersects, "[[", 3)), NA),
        encounter = unlist(lapply(cam_intersects, "[[", 4)),
        t_in = c(unlist(lapply(cam_intersects, "[[", 5)), NA),
        t_out = c(unlist(lapply(cam_intersects, "[[", 6)), NA),
        lscape_index = lscape_index[1],
        Animal_ID = Animal_ID,
        t = t,
        ii = ii,
        X = X,
        Y = Y,
        trav_speeds = trav_speeds,
        Speed = lscape_type[1],
        .groups = "drop"
      ) |>
      dplyr::filter(t_stay > 0) |>
      dplyr::group_by(pass_i) |>
      dplyr::mutate(pass_i = pass_i + 0.001 * cumsum(c(1, abs(encounter[-n()] - encounter[-1]) > 0))) |>
      dplyr::select(-cam_intersects) |>
      ungroup()
  } else {
    cam_captures <- cell_captures
  }
  return(cam_captures)
}

#' Aggregate Animal Counts per Camera
#'
#' @description
#' Summarizes the camera capture data to determine the total number of animals
#' detected at each camera location.
#'
#' @param cam_locs A data frame of camera locations.
#' @param cam_captures A data frame of processed captures (output from \code{\link{get_cam_captures}}).
#' @param animalxy A data frame of raw trajectories (used for verification/context, optional).
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{cam_ID}: Camera identifier.
#'   \item \code{lscape_index}: Landscape cell index.
#'   \item \code{count}: Total number of animals detected.
#'   \item \code{mean_group}: Average group size detected.
#' }
#' @export
#'
get_count_data <- function(cam_locs, cam_captures, animalxy) {
  count_data <- cam_captures |>
    dplyr::filter(in_cam == T) |>
    dplyr::group_by(lscape_index, t) |>
    dplyr::summarise(
      group_size = n(),
      .groups = "drop"
    ) |>
    dplyr::group_by(lscape_index) |>
    dplyr::summarise(
      count = sum(group_size),
      mean_group = mean(group_size),
      .groups = "drop"
    ) |>
    dplyr::full_join(cam_locs |>
                       select(cam_ID, lscape_index, Speed),
                     by = c("lscape_index")) |>
    dplyr::arrange(cam_ID) |>
    dplyr::mutate(count = replace(count, is.na(count), 0))

  return(count_data)
}

#' Aggregate Encounter Frequencies
#'
#' @description
#' Calculates the number of unique "passes" (distinct encounter events) at each
#' camera location, regardless of group size.
#'
#' @param cam_locs A data frame of camera locations.
#' @param cam_captures A data frame of processed captures (output from \code{\link{get_cam_captures}}).
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{lscape_index}: Landscape cell index.
#'   \item \code{encounter}: Number of unique passing events.
#'   \item \code{cam_ID}: Camera identifier.
#'   \item \code{Speed}: Landscape speed category.
#' }
get_encounter_data <- function(cam_locs, cam_captures) {
  encounter_data <- cam_captures |>
    dplyr::group_by(lscape_index) |>
    dplyr::summarise(
      encounter = length(unique(pass_i)),
      .groups = "drop"
    ) |>
    dplyr::full_join(cam_locs |>
                       select(cam_ID, lscape_index, Speed),
                     by = c("lscape_index")) |>
    dplyr::arrange(cam_ID) |>
    dplyr::mutate(encounter = replace(encounter, is.na(encounter), 0))

  return(encounter_data)
}

#' Calculate and Normalize Stay Times
#'
#' @description
#' Computes the duration animals stay within the camera's view. It normalizes
#' these times based on landscape speed categories (Slow, Medium, Fast) and
#' formats the data for downstream analysis.
#'
#' @param cam_locs A data frame of camera locations.
#' @param cam_captures A data frame of processed captures.
#'
#' @return A list containing two objects:
#' \itemize{
#'   \item \code{[[1]]} (stay_time_raw): A data frame of raw stay times per pass.
#'   \item \code{[[2]]} (stay_time_data): A wide-format data frame where columns represent
#'   individual encounters and values represent normalized stay times.
#' }
#'
get_stay_time_data <- function(cam_locs, cam_captures) {

  stay_time_raw <- cam_captures |>
    dplyr::group_by(pass_i) |>
    dplyr::summarise(
      t_stay = sum(t_stay),
      lscape_index = lscape_index[1],
      Speed = Speed[1],
      .groups = "drop"
    )

  stay_time_ref <- tibble::tibble(
    Speed = c("Slow", "Medium", "Fast"),
    mean_stay = 1
  ) %>%
    dplyr::filter(
      !(Speed %in% unique(stay_time_raw$Speed))
    )

  stay_time_normalize <- stay_time_raw %>%
    dplyr::group_by(Speed) %>%
    dplyr::summarise(
      mean_stay = mean(t_stay),
      .groups = 'drop'
    ) %>%
    dplyr::bind_rows(stay_time_ref) %>%
    dplyr::summarise(
      sum_stay = sum(mean_stay)
    ) %>%
    dplyr::pull(sum_stay)

  # Format staying time data
  stay_time_data <- stay_time_raw %>%
    dplyr::mutate(t_stay = t_stay / stay_time_normalize) %>%
    # dplyr::add_row(lscape_index = cam_locs$lscape_index[cam_locs$lscape_index %notin% cam_captures$lscape_index]) |>
    dplyr::full_join(cam_locs |>
                       select(cam_ID, lscape_index),
                     by = "lscape_index") |>
    dplyr::arrange(cam_ID) |>
    dplyr::select(lscape_index, t_stay) |>
    dplyr::group_by(lscape_index) |>
    dplyr::mutate(encounter = 1:n()) |>
    dplyr::ungroup() |>
    pivot_wider(names_from = encounter, values_from = t_stay) |>
    dplyr::select(-lscape_index)

  return(list(stay_time_raw, stay_time_data))
}


