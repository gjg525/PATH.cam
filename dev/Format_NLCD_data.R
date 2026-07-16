library(ggplot2)

tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5000NLCDclip.tif"
tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5004NLCDclip.tif"
tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5006NLCDclip.tif"
tif_filename <- "G:/My Drive/Missoula_postdoc/PATH_model/NLCD_data/LowTag5010NLCDclip.tif"

# buildBackground <- function(parameterSet, tifFile) {
#   # Read raster and convert to matrix matching MATLAB's orientation
#   r <- terra::rast(tifFile)
#   AA <- terra::as.matrix(r, wide = TRUE)
#
#   # flipud
#   BB <- AA[nrow(AA):1, ]
#
#   Background <- BB[331:362, 548:ncol(BB)]
#
#   water <- 11
#   development <- c(21, 22, 23, 24)
#   forest <- c(31, 41, 42, 43, 52)
#   agland <- c(71, 81, 82, 90, 95)
#
#   Landscape <- matrix(0, nrow = nrow(Background), ncol = ncol(Background))
#   Landscape[Background %in% water] <- 1
#   Landscape[Background %in% development] <- 2
#   Landscape[Background %in% forest] <- 3
#   Landscape[Background %in% agland] <- 4
#
#   habTypes <- sort(unique(as.vector(Landscape)))
#   habCount <- numeric(length(habTypes))
#
#   for (ii in seq_along(habTypes)) {
#     habCount[ii] <- sum(Landscape == habTypes[ii])
#   }
#
#   ParSets <- array(0, dim = c(nrow(Landscape), ncol(Landscape), 13))
#   Mu <- matrix(0, nrow = nrow(Landscape), ncol = ncol(Landscape))
#
#   for (mm in seq_along(habTypes)) {
#     mask <- as.numeric(Landscape == habTypes[mm])
#
#     # Mu <- Mu + parameterSet[1, mm] * mask
#     Mu <- Mu + parameterSet[mm] * mask
#     # ParSets[,,1]  <- ParSets[,,1] + parameterSet[1, mm] * mask
#     # ParSets[,,2]  <- ParSets[,,2] + parameterSet[2, mm] * mask
#     # ParSets[,,3]  <- ParSets[,,3] + parameterSet[3, mm] * mask
#     # ParSets[,,4]  <- ParSets[,,4] + parameterSet[4, mm] * mask
#     # ParSets[,,5]  <- ParSets[,,5] + parameterSet[5, mm] * mask
#     # ParSets[,,6]  <- ParSets[,,6] + parameterSet[6, mm] * mask
#     # ParSets[,,7]  <- ParSets[,,7] + parameterSet[7, mm] * mask
#     # ParSets[,,8]  <- ParSets[,,8] + parameterSet[8, mm] * mask
#     # ParSets[,,9]  <- ParSets[,,9] + parameterSet[9, mm] * mask
#     # ParSets[,,10] <- ParSets[,,10] + parameterSet[10, mm] * mask
#     # ParSets[,,11] <- ParSets[,,11] + parameterSet[11, mm] * mask
#     # ParSets[,,12] <- ParSets[,,12] + parameterSet[12, mm] * mask
#     # ParSets[,,13] <- ParSets[,,13] + parameterSet[13, mm] * mask
#   }
#
#   # mu_vec <- parameterSet[1, seq_along(habTypes)]
#   # mubar <- sum(habCount) / sum(habCount / mu_vec)
#   # mubar <- 0.001 * mubar
#
#   # return(list(Landscape = Landscape, ParSets = ParSets, Mu = Mu, mubar = mubar))
#   return(list(Landscape = Landscape, Mu = Mu))
# }

# --------------------------------------------------- #
# 1. Initialize parameters
# --------------------------------------------------- #
params <- list(
  scale_diffusion = 1
  # scale_birth = 1,
  # scale_infection_I = 1,
  # scale_infection_H = 1,
  # scale_natural_mortality = 1,
  # scale_disease_mortality = 1,
  # scale_env = 1,
  # scale_cull_I_rate = 0,
  # scale_cull_S_rate = 0,
  # scale_cull_threshold = 1
)

mu_base     <- c(4, 2, 0.02, 0.5)
# alpha_base  <- c(.01, .01, .01, .01)
# betaI_base  <- c(.7, .5, .1, .2)
# betaH_base  <- c(.01, .01, .01, .01)
# omega_base  <- c(.01, .01, .01, .01)
# lambda_base <- c(.1, .1, .1, .1)
# theta_base  <- c(1, 1, 1, 1)
# gamma_base  <- c(.01, .01, .01, .01)
# eta_base    <- c(.00001, .00001, .00001, .00001)
# nu_base     <- c(.01, .01, .01, .01)
# kappa_base  <- c(0.01, 0.01, 0.01, 0.01)
# phi_base    <- c(1.0, 1.0, 1.0, 1.0)

# Scale parameters
mu     <- params$scale_diffusion * mu_base
# alpha  <- params$scale_birth * alpha_base
# betaI  <- params$scale_infection_I * betaI_base
# betaH  <- params$scale_infection_H * betaH_base
# omega  <- params$scale_natural_mortality * omega_base
# lambda <- params$scale_disease_mortality * lambda_base
# theta  <- params$scale_env * theta_base
# gamma  <- params$scale_env * gamma_base
# eta    <- params$scale_env * eta_base
# nu     <- params$scale_env * nu_base
# kappaI <- params$scale_cull_I_rate * kappa_base
# kappaS <- params$scale_cull_S_rate * kappa_base
# phi    <- params$scale_cull_threshold * phi_base

parameterSet <- mu #rbind(mu, alpha, betaI, betaH, omega, lambda,
                   #    theta, gamma, eta, nu, kappaI, kappaS, phi)

# --------------------------------------------------- #
# 2. Build Background
# --------------------------------------------------- #
bg_info <- buildBackground(parameterSet, tifFile = tif_filename)

# Plot
# Convert the matrix to a long-format data frame
df <- reshape2::melt(bg_info$Landscape)
colnames(df) <- c("Row", "Column", "LandCover")

# Convert the numeric values (1-4) into categorical factors with labels
df$LandCover <- factor(df$LandCover,
                       levels = c(1, 2, 3, 4),
                       labels = c("Water", "Development", "Forest", "Agriculture"))

landscape_plot <- ggplot(df, aes(x = Column, y = Row, fill = LandCover)) +
  geom_raster() +
  # Assign professional, colorblind-friendly or intuitive colors
  scale_fill_manual(
    values = c(
      "Water"       = "#1F78B4",  # A nice deep blue
      "Development" = "#E31A1C",  # Red (common for development) or use "#666666" for Grey
      "Forest"      = "#33A02C",  # Deep green
      "Agriculture" = "#FDBF6F"   # Earthy yellow/tan
    ),
    # values = c(
    #   "Water"       = "#FFFFFF",  # A nice deep blue
    #   "Development" = "#AAAAAA",  # Red (common for development) or use "#666666" for Grey
    #   "Forest"      = "#555555",  # Deep green
    #   "Agriculture" = "#000000"   # Earthy yellow/tan
    # ),
    drop = FALSE, # Ensures all categories appear in the legend even if 0 pixels exist
    name = "Land Cover Type"
  ) +
  # Matrices plot from bottom-up in ggplot, so we reverse the Y-axis to match your matrix view
  scale_y_reverse() +
  # Ensure the grid cells are perfectly square
  coord_fixed() +
  # Add titles (optional for scientific papers, often handled in the caption)
  labs(
    x = "Easting (Grid Columns)",
    y = "Northing (Grid Rows)"
  ) +
  # Apply a clean, minimalist theme suitable for publications
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(), # Remove grid lines
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Add a neat bounding box
    axis.text = element_blank(),  # Remove axis numbers (optional, remove if you want coordinates)
    axis.ticks = element_blank(), # Remove axis ticks
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(color = "black", size = 0.2) # Add thin borders around legend keys
  )

# View the plot
print(landscape_plot)

# # Save as a high-resolution PNG
# ggsave("Landscape_Classification.png",
#        plot = landscape_plot,
#        width = 8,
#        height = 5,
#        dpi = 300,
#        bg = "white")
#
# # Save as a PDF (often preferred for vector graphics in papers)
# ggsave("Landscape_Classification.pdf",
#        plot = landscape_plot,
#        width = 8,
#        height = 5,
#        bg = "white")
