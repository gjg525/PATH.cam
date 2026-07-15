buildBackground <- function(parameterSet, tifFile) {
  # Read raster and convert to matrix matching MATLAB's orientation
  r <- terra::rast(tifFile)
  AA <- terra::as.matrix(r, wide = TRUE)

  # flipud
  BB <- AA[nrow(AA):1, ]

  # Background <- BB[331:362, 548:ncol(BB)]
  Background <- BB[333:362, 548:(ncol(BB) - 2)]

  water <- 11
  development <- c(21, 22, 23, 24)
  forest <- c(31, 41, 42, 43, 52)
  agland <- c(71, 81, 82, 90, 95)

  Landscape <- matrix(0, nrow = nrow(Background), ncol = ncol(Background))
  Landscape[Background %in% water] <- 1
  Landscape[Background %in% development] <- 2
  Landscape[Background %in% forest] <- 3
  Landscape[Background %in% agland] <- 4

  habTypes <- sort(unique(as.vector(Landscape)))
  habCount <- numeric(length(habTypes))

  for (ii in seq_along(habTypes)) {
    habCount[ii] <- sum(Landscape == habTypes[ii])
  }

  ParSets <- array(0, dim = c(nrow(Landscape), ncol(Landscape), 13))
  Mu <- matrix(0, nrow = nrow(Landscape), ncol = ncol(Landscape))

  for (mm in seq_along(habTypes)) {
    mask <- as.numeric(Landscape == habTypes[mm])

    Mu <- Mu + parameterSet[mm] * mask
  }

  # return(list(Landscape = Landscape, ParSets = ParSets, Mu = Mu, mubar = mubar))
  return(list(Landscape = Landscape, Mu = Mu))
}
