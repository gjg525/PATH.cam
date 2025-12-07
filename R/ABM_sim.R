#' Run Agent-Based Model Simulation for Animal Movement
#'
#' @description
#' Simulates animal movement across a heterogeneous landscape using a parallelized
#' agent-based model (ABM). Animals move in groups based on correlated random walks
#' (CRW), adjusting their speed and direction based on landscape features and
#' boundary interactions.
#'
#' @details
#' The simulation utilizes parallel processing (defaulting to 3 cores) to iterate
#' through animal groups. Movement directions are determined using a Von Mises
#' distribution (`Rfast::rvonmises`). If an animal crosses a grid cell boundary
#' into a new landscape type, its speed is adjusted using a truncated normal
#' distribution based on the new cell's friction/speed value.
#'
#' Boundary conditions are reflective: if an animal hits the simulation bounds,
#' it reflects back into the simulation area.
#'
#' @param study_design A list containing the simulation parameters. Must include:
#' \itemize{
#'   \item \code{dx}, \code{dy}: Dimensions of grid cells (numeric).
#'   \item \code{q}: Total number of grid cells (numeric). Assumes a square grid where sqrt(q) is an integer.
#'   \item \code{t_steps}: Number of time steps to simulate (integer).
#'   \item \code{group_sizes}: Vector of integers defining the size of each animal group.
#'   \item \code{group_spread}: Initial spatial spread of groups (numeric).
#'   \item \code{dt}: Time increment per step (numeric).
#'   \item \code{bounds}: Vector defining the spatial extent (e.g., \code{c(min, max)}).
#' }
#' @param lscape_defs A data frame defining the landscape properties. Must include:
#' \itemize{
#'   \item \code{Value}: Base movement speeds associated with landscape cells.
#'   \item \code{Speed}: (Optional) Speed category labels used for the final output.
#' }
#' @param animalxy.0 (Optional) A data frame of initial animal positions. If \code{NULL} (default),
#'   positions are generated using internal helper functions \code{place_animals}
#'   and \code{create_activity_mat}.
#'
#' @return A data frame (tibble) containing the complete simulation trajectory with columns:
#' \itemize{
#'   \item \code{Animal_ID}: Unique identifier for the animal.
#'   \item \code{group_ID}: Identifier for the group.
#'   \item \code{X}, \code{Y}: Spatial coordinates.
#'   \item \code{t}: Time stamp.
#'   \item \code{lscape_index}: Index of the grid cell the animal is currently in.
#'   \item \code{trav_dist}: Distance traveled from the previous step.
#'   \item \code{trav_speeds}: Calculated speed for the step.
#'   \item \code{lscape_type}: The speed category of the current landscape cell.
#' }
#'
#' @export
#'
ABM_sim <- function(
    study_design,
    lscape_defs,
    animalxy.0 = NULL){

  dx <- study_design$dx
  dy <- study_design$dy
  q <- study_design$q
  t_steps <- study_design$t_steps
  speeds <- matrix(lscape_defs$Value, q^0.5, q^0.5)
  group_sizes <- unlist(study_design$group_sizes)
  group_spread <- study_design$group_spread
  dt <- study_design$dt
  bounds <- unlist(study_design$bounds)

  if (is.null(animalxy.0)) {
    # Place animals on landscape
    animalxy.0 <- place_animals(study_design, lscape_defs)

    # Create active/non-active switch for each animal
    animalxy.0 <- create_activity_mat(study_design, animalxy.0)
  }

  # corner locations of each grid cell
  grid.ints <- seq(min(bounds), max(bounds), by = max(bounds) / q^0.5)

  # register cluster to be used by %dopar%
  # n_cores <- parallel::detectCores() # Check number of cores available
  my_cluster <- parallel::makeCluster(3, type = "PSOCK")
  doParallel::registerDoParallel(cl = my_cluster)

  # Loop through groups in parallel
  animalxy.par <- foreach::foreach(nc = 1:length(group_sizes)) %dopar% {
    # Get group size for current group
    group_size <- group_sizes[nc]

    # Temporary location matrix for each group of animals
    group.par <- animalxy.0 |>
      dplyr::filter(group_ID == nc) |>
      dplyr::mutate(
        t = 0
      ) |>
      dplyr::select(Animal_ID, group_ID, X, Y, t)

    # active/non-active matrix
    activity_par <- animalxy.0 |>
      dplyr::filter(group_ID == nc) |>
      dplyr::pull(activity_mat)

    # Define group turning angle
    theta_group <- runif(1, 0, 2 * pi)

    # Initialize temporary locations for all animals in group
    group_X <- matrix(NA, group_size, t_steps + 1)
    group_Y <- matrix(NA, group_size, t_steps + 1)

    group_X[, 1] <- group.par$X
    group_Y[, 1] <- group.par$Y

    # Loop through all time steps
    for (i in 1:t_steps) {
      # Shift turning angle with turning bias
      theta_group <- Rfast::rvonmises(
        1,
        theta_group,
        0
      )

      # Define turning angles for all individuals
      theta.all <- Rfast::rvonmises(group_size, theta_group, 20)

      # Loop through all individuals in group
      for (ci in 1:group_size) {
        # Time step
        t.step <- dt

        active_par <- unlist(activity_par[ci])

        step_count <- 0
        # Individual moves across landscape until time is depleted
        while (t.step > 0) {
          # x,y indices of individual
          X.ind <- ceiling(group_X[ci, i] * (q^0.5 / max(bounds)))
          Y.ind <- ceiling(group_Y[ci, i] * (q^0.5 / max(bounds)))

          # Determine movement speed in current space
          # Note: Come back to define speeds with upper/lower bounds
          # v <- speeds[X.ind, Y.ind]
          v <- truncnorm::rtruncnorm(1,
                                     a = 0,
                                     b = Inf,
                                     speeds[X.ind, Y.ind],
                                     speeds[X.ind, Y.ind] / 10)
          step.size <- v * t.step

          # Step length components
          dX <- step.size * cos(as.numeric(theta.all[ci]))
          dY <- step.size * sin(as.numeric(theta.all[ci]))

          # Determine coordinate `corners' of grid for line segment calculations
          r.corners <- rbind(
            cbind(grid.ints[X.ind], grid.ints[Y.ind]),
            cbind(grid.ints[X.ind], grid.ints[Y.ind]),
            cbind(grid.ints[X.ind], grid.ints[Y.ind + 1]),
            cbind(grid.ints[X.ind + 1], grid.ints[Y.ind])
          )

          # step lengths to add to corners
          r.length <- rbind(cbind(0, dy), cbind(dx, 0), cbind(dx, 0), cbind(0, dy))

          # check for intersections with grid
          # Don't use intersections when start point is on the line (a,b = 0),
          # or if a full step lands on a line (a,b = 1)
          a <- c()
          b <- c()
          int.check <- c()
          for (rr in 1:4) {
            a[rr] <- det(rbind((c(group_X[ci, i], group_Y[ci, i]) - r.corners[rr, ]), r.length[rr, ])) /
              det(rbind(r.length[rr, ], c(dX, dY)))
            b[rr] <- det(rbind((c(group_X[ci, i], group_Y[ci, i]) - r.corners[rr, ]), c(dX, dY))) /
              det(rbind(r.length[rr, ], c(dX, dY)))
          }

          int.check <- a < 1 & a > 0 & b < 1 & b > 0 |
            (a == 0 & b < 1 & b > 0) |
            (a == 1 & b < 1 & b > 0) |
            (b == 0 & a < 1 & a > 0) |
            (b == 1 & a < 1 & a > 0)

          # If individual crosses into new landscape cell, change movement speed
          if (any(int.check == T)) {
            # only use first intersection crossing
            if (sum(int.check == T) > 1) {
              # coordinates for each intersection
              int.coords.all <- t(matrix(c(group_X[ci, i], group_Y[ci, i]), 2, sum(int.check == T))) +
                a[int.check == TRUE] * t(matrix(c(dX, dY), 2, sum(int.check == T)))

              vec.dist.all <- sqrt((t(matrix(c(group_X[ci, i]), 1, sum(int.check == T))) - int.coords.all[, 1])^2 +
                (t(matrix(c(group_Y[ci, i]), 1, sum(int.check == T))) - int.coords.all[, 2])^2)
              int.T <- which(int.check == T)
              int.ind <- int.T[which(vec.dist.all == min(vec.dist.all))]
            } else {
              int.ind <- which(int.check == T)
            }

            # X,Y coordinates at intersection
            temp.X <- r.corners[int.ind, 1] + b[int.ind] * r.length[int.ind, 1]
            temp.Y <- r.corners[int.ind, 2] + b[int.ind] * r.length[int.ind, 2]

            # Set individual near grid line, adjust for new speed
            # Determine whether the lscape type is "blocked" or not
            # Do this better
            tx <- ceiling((temp.X + 0.00001 * sign(cos(theta.all[ci])) *
              (temp.X == r.corners[int.ind, 1])) * (q^0.5 / max(bounds)))
            ty <- ceiling((temp.Y + 0.00001 * sign(sin(theta.all[ci])) *
              (temp.Y == r.corners[int.ind, 2])) * (q^0.5 / max(bounds)))
            tx <- min(c(tx, q^0.5))
            tx <- max(c(tx, 1))
            ty <- min(c(ty, q^0.5))
            ty <- max(c(ty, 1))
            if (temp.X < (max(bounds)) & temp.X > (min(bounds)) &
              temp.Y < (max(bounds)) & temp.Y > (min(bounds))
            ) {
              t.step <- t.step - sqrt((group_X[ci, i] - (temp.X))^2 +
                                        (group_Y[ci, i] - (temp.Y))^2) / v

              # perturb individual so that they are not on intersection
              group_X[ci, i] <- temp.X + 0.00001 * sign(cos(theta.all[ci])) * (temp.X == r.corners[int.ind, 1])
              group_Y[ci, i] <- temp.Y + 0.00001 * sign(sin(theta.all[ci])) * (temp.Y == r.corners[int.ind, 2])

              # Track all intermediate steps
              group.par[nrow(group.par) + 1, ] <- cbind(
                min(group.par$Animal_ID) + ci - 1,
                group.par$group_ID[1],
                group_X[ci, i],
                group_Y[ci, i],
                i * dt - t.step
              )
            } else {
              if (temp.X >= (max(bounds)) || temp.X <= (min(bounds))) {
                t.step <- t.step - sqrt((group_X[ci, i] - (temp.X))^2 +
                  (group_Y[ci, i] - (temp.Y))^2) / v

                theta.x <- cos(theta.all[ci])
                theta.y <- sin(theta.all[ci])

                # Reflective angle
                theta.all[ci] <- (atan(-theta.y / theta.x) + pi * (theta.x > 0)) %% (2 * pi)

                # perturb individual so that they are not on intersection
                group_X[ci, i] <- temp.X + 0.00001 * sign(cos(theta.all[ci])) * (temp.X == r.corners[int.ind, 1])
                group_Y[ci, i] <- temp.Y + 0.00001 * sign(sin(theta.all[ci])) * (temp.Y == r.corners[int.ind, 2])

                # Track all intermediate steps
                group.par[nrow(group.par) + 1, ] <- cbind(
                  min(group.par$Animal_ID) + ci - 1,
                  group.par$group_ID[1],
                  group_X[ci, i],
                  group_Y[ci, i],
                  i * dt - t.step
                )
                # Designate one individual to determine group turning angle
                if (ci == 1) {
                  theta_group <- theta.all[ci]
                }
              }
              if (temp.Y >= (max(bounds)) || temp.Y <= (min(bounds))) {
                t.step <- t.step - sqrt((group_X[ci, i] - (temp.X))^2 +
                  (group_Y[ci, i] - (temp.Y))^2) / v

                theta.x <- cos(theta.all[ci])
                theta.y <- sin(theta.all[ci])
                # Reflective angle
                theta.all[ci] <- (atan(-theta.y / theta.x) + pi * (theta.x < 0)) %% (2 * pi)

                # perturb individual so that they are not on intersection
                group_X[ci, i] <- temp.X + 0.00001 * sign(cos(theta.all[ci])) * (temp.X == r.corners[int.ind, 1])
                group_Y[ci, i] <- temp.Y + 0.00001 * sign(sin(theta.all[ci])) * (temp.Y == r.corners[int.ind, 2])

                group.par[nrow(group.par) + 1, ] <- cbind(
                  min(group.par$Animal_ID) + ci - 1,
                  group.par$group_ID[1],
                  group_X[ci, i],
                  group_Y[ci, i],
                  i * dt - t.step
                )

                # designate one individual to determine group movement
                if (ci == 1) {
                  theta_group <- theta.all[ci]
                }
              }
            }
          } else {
            # Individuals take step within cell
            group_X[ci, i + 1] <- group_X[ci, i] + dX
            group_Y[ci, i + 1] <- group_Y[ci, i] + dY
            t.step <- 0

            # Track all intermediate steps
            group.par[nrow(group.par) + 1, ] <- cbind(
              min(group.par$Animal_ID) + ci - 1,
              group.par$group_ID[1],
              group_X[ci, i + 1],
              group_Y[ci, i + 1],
              i * dt
            )
          }

          step_count <- step_count + 1
          if (step_count > 2000) {
            stop("Animal stuck in simulation. Check for convergent landscape directions.")
          }
        }
      }
      # group direction is determined from 1 representative animal
      theta_group <- theta.all[1]
    }
    animalxy.par <- group.par
  }

  # Stop cluster
  stopCluster(my_cluster)

  # Convert list to data frame
  animalxy.all <- dplyr::bind_rows(animalxy.par)

  # Convert animalxy 2D coordinates to single digit coordinates
  animalxy.all <- animalxy.all |>
    dplyr::ungroup() |>
    dplyr::mutate(
      lscape_index = ceiling(animalxy.all$X / dx) +
        floor(animalxy.all$Y / dx) * q^0.5,
      ii = 1:nrow(animalxy.all)
  ) |>
    dplyr::group_by(Animal_ID) |>
    dplyr::mutate(
      trav_dist = c(NA, sqrt((X[2:n()] - X[1:(n() - 1)])^2 +
        (Y[2:n()] - Y[1:(n() - 1)])^2)),
      trav_speeds = trav_dist / c(NA, t[2:n()] - t[1:(n() - 1)])
    )

  # Escapee
  if (max(animalxy.all$X) > max(bounds) || min(animalxy.all$X) < min(bounds)) {
    stop("Animal out of bounds")
  }

  # Add speed info
  animalxy.all$lscape_type <- lscape_defs$Speed[animalxy.all$lscape_index]

  return(animalxy.all)
  }
