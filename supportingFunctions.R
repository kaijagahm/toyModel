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
                           add00 = c(0.5, 7), 
                           lose01 = 0.3, 
                           add10 = 0.3, 
                           lose11 = c(0.3, 0.3)){ 
  
  # argument checks
  checkmate::assertNumeric(add00, len = 2)
  checkmate::assertNumeric(lose11, len = 2)
  checkmate::assertNumeric(lose01, len = 1)
  checkmate::assertNumeric(add10, len = 1)
  checkmate::assertList(network.history)
  checkmate::assertNumeric(ind, len = 1)
  
  # Establish a network history, two steps back
  if(ind == 1){
    stop("network must have at least some history--ind cannot be = 1")
  }else if(ind == 2){
    prev <- network.history[[ind-1]]
    # create a matrix of zeroes if the history doesn't exist that far back
    prevprev <- matrix(data = 0, 
                       nrow = nrow(prev), ncol = ncol(prev))
    new <- prev
  }else if(ind > 2){
    prev <- network.history[[ind-1]]
    prevprev <- network.history[[ind-2]]
    new <- prev
  }
  
  # sort edges by history, two back
  h00 <- dedup(which(prev == prevprev & prev == 0, arr.ind = T), "upper")
  h11 <- dedup(which(prev == prevprev & prev == 1, arr.ind = T), "upper")
  h01 <- dedup(which(prevprev < prev, arr.ind = T), "upper")
  h10 <- dedup(which(prevprev > prev, arr.ind = T), "upper")
  N <- nrow(prev)

  # Modify the new adjacency matrix (upper triangle only)
  new[h00] <- rbinom(n = nrow(h00), size = 1, 
                     prob = rbeta(n = nrow(h00), 
                                  shape1 = add00[1], shape2 = add00[2]))
  new[h11] <- rbinom(n = nrow(h11), size = 1,
                     prob = rbeta(n = nrow(h11),
                                  shape1 = lose11[1], shape2 = lose11[2]))
  new[h01] <- rbinom(n = nrow(h01), size = 1,
                     prob = lose01)
  new[h10] <- rbinom(n = nrow(h10), size = 1,
                     prob = add10)
  
  # Symmetrize the matrix
  new <- as.matrix(sna::symmetrize(new, rule = "upper")) # copy the upper triangle over the lower triangle
  
  return(new)
}
