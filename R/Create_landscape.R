#' Generate a Randomized Landscape of Movement Speeds
#'
#' @description
#' Creates a spatial grid representing the landscape and assigns random movement
#' speed categories to each grid cell.
#'
#' @details
#' The function generates a square grid based on the total number of cells (`q`)
#' defined in `study_design`. It then randomly samples (with replacement) from
#' the speed categories defined in `lscape_design` to assign a landscape type
#' (e.g., "Forest", "Open") to each cell. Finally, it draws a specific continuous
#' speed value for each cell from a uniform distribution bounded by that category's
#' min/max limits.
#'
#' @param study_design A list containing simulation parameters. Must include:
#' \itemize{
#'   \item \code{q}: Total number of grid cells (integer). The function assumes a
#'   square grid, so \code{sqrt(q)} should ideally be an integer.
#' }
#' @param lscape_design A list (or data frame) defining the speed categories. Must include vectors for:
#' \itemize{
#'   \item \code{Speed_ID}: Names of the speed categories (e.g., "Fast", "Slow").
#'   \item \code{Speed_mins}: Minimum speed values for each category.
#'   \item \code{Speed_maxes}: Maximum speed values for each category.
#' }
#'
#' @return A tibble representing the landscape grid with columns:
#' \itemize{
#'   \item \code{Index}: Cell ID.
#'   \item \code{X}, \code{Y}: Grid coordinates.
#'   \item \code{Speed}: Assigned speed category.
#'   \item \code{Min}, \code{Max}: Bounds for the assigned category.
#'   \item \code{Value}: The specific randomized speed value for that cell.
#' }
#'
#' @export
lscape_creator <- function(study_design, lscape_design) {
  q <- study_design$q

  speed_bounds <- tibble::tibble(
    Speed = unlist(lscape_design$Speed_ID),
    Min = unlist(lscape_design$Speed_mins),
    Max = unlist(lscape_design$Speed_maxes)
  )

  lscape_defs <- tibble::tibble(
    Index = 1:q,
    X = (Index - 1) %% q^0.5 + 1,
    Y = ceiling(Index / q^0.5)
  )

  # Create dataframe for randomly-distributed speeds
  lscape_defs <- lscape_defs |>
    dplyr::bind_cols(dplyr::sample_n(speed_bounds, q, replace = T)) |>
    dplyr::mutate(
      Value = runif(
        q,
        Min,
        Max
      )
    )

  return(lscape_defs)
}
