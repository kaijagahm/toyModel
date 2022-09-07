# Supporting functions for the modelCode_rewiring.R script.
# AUTHOR: Kaija Gahm
# MODIFIED: 2022-06-15
# Just separating these out to keep the other code cleaner and avoid too much scrolling.

library(checkmate)
library(dplyr)

# Define functions --------------------------------------------------------

# Small utility functions for deduplicating edges
# uniqueEdges -------------------------------------------------------------
uniqueEdges <- function(n, triangle = "upper"){
  checkmate::assertNumeric(n, len = 1)
  df <- expand.grid(from = 1:n, to = 1:n)
  if(triangle == "upper"){
    df <- df[which(df[,1] < df[,2]), , drop = FALSE] 
  }else if(triangle == "lower"){
    df <- df[which(df[,1] < df[,2]), , drop = FALSE] 
  }else{
    stop("Triangle must be 'upper' or 'lower'")
  }
  return(df)
}

# dedup -------------------------------------------------------------------
dedup <- function(mat, triangle = "upper"){
  checkmate::assertMatrix(mat)
  if(triangle == "upper"){
    mat <- mat[which(mat[,1] < mat[,2]), , drop = FALSE]
  }else if(triangle == "lower"){
    mat <- mat[which(mat[,1] > mat[,2]), , drop = FALSE]
  }else{
    stop("Triangle must be 'upper' or 'lower'")
  }
  return(mat)
}

# update.network ----------------------------------------------------------
# Function for baseline network dynamics
# Repeat this number of times specified for desired burn.in, in a for loop. Each time, spitting out the network, and the history of the edges.
update.network <- function(ind, # starting index for the history list.
                           network.history, # list of history, since this function depends on being able to look a few timesteps back.
                           mod00 = -0.2, 
                           mod11 = 0.2, 
                           mod10 = -0.1, 
                           mod01 = 0.1,
                           probMatrix = NULL){ 
  
  # argument checks
  checkmate::assertNumeric(mod00, len = 1)
  checkmate::assertNumeric(mod11, len = 1)
  checkmate::assertNumeric(mod10, len = 1)
  checkmate::assertNumeric(mod01, len = 1)
  checkmate::assertList(network.history)
  checkmate::assertNumeric(ind, len = 1)
  N <- nrow(network.history[[ind-1]])
  checkmate::assertMatrix(probMatrix, null.ok = FALSE, nrows = N, ncols = N)
  
  # Establish a network history, two steps back
  if(ind == 1){
    stop("network must have at least some history--ind cannot be = 1")
  }else if(ind == 2){
    prev <- network.history[[ind-1]]
    # create a matrix of zeroes if the history doesn't exist that far back
    prevprev <- matrix(data = 0, 
                       nrow = nrow(prev), ncol = ncol(prev))
  }else if(ind > 2){
    prev <- network.history[[ind-1]]
    prevprev <- network.history[[ind-2]]
  }
  new <- probMatrix
  
  # sort edges by history, two back
  h00 <- dedup(which(prev == prevprev & prev == 0, arr.ind = T), "upper")
  h11 <- dedup(which(prev == prevprev & prev == 1, arr.ind = T), "upper")
  h01 <- dedup(which(prevprev < prev, arr.ind = T), "upper")
  h10 <- dedup(which(prevprev > prev, arr.ind = T), "upper")

  # Modify the new probability matrix (upper triangle only)
  new[h00] <- new[h00]+(new[h00]*mod00)
  new[h11] <- new[h11]+(new[h11]*mod11)
  new[h10] <- new[h10]+(new[h10]*mod10)
  new[h01] <- new[h01]+(new[h01]*mod01)
  
  # Fix any probabilities that are too small or too big
  new[new > 1] <- 1
  new[new < 0] <- 0
  
  # Create adjacency matrix using rbinom
  newAdj <- matrix(rbinom(N*N, 1, new), N, N)
  
  # Symmetrize the matrix
  newAdj <- as.matrix(sna::symmetrize(newAdj, rule = "upper")) # copy the upper triangle over the lower triangle
  
  return(newAdj)
}
