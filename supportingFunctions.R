# Supporting functions for the modelCode_rewiring.R script.
# AUTHOR: Kaija Gahm
# MODIFIED: 2022-06-15
# Just separating these out to keep the other code cleaner and avoid too much scrolling.

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(data.table)
library(checkmate)
library(sna)

# Define functions --------------------------------------------------------

# Small utility functions for deduplicating edges
# uniqueEdges -------------------------------------------------------------
uniqueEdges <- function(n, triangle = "upper"){
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
dedup <- function(df, triangle = "upper"){
  if(triangle == "upper"){
    df <- df[which(df[,1] < df[,2]), , drop = FALSE]
  }else if(triangle == "lower"){
    df <- df[which(df[,1] > df[,2]), , drop = FALSE]
  }else{
    stop("Triangle must be 'upper' or 'lower'")
  }
  return(df)
}

# update.network ----------------------------------------------------------
# Function for baseline network dynamics
# Repeat this number of times specified for desired burn.in, in a for loop. Each time, spitting out the network, and the history of the edges.
update.network <- function(ind, # starting index for the history list.
                           network.history, # list of history, since this function depends on being able to look a few timesteps back.
                           add00 = 0.1, 
                           lose01 = 0.3, 
                           add10 = 0.3, 
                           lose11 = 0.1){ 
  
  # establish history two steps back
  prev <- network.history[[ind-1]]
  prevprev <- network.history[[ind-2]]
  new <- prev
  
  # sort edges by history, two back
  h00 <- dedup(which(prev == prevprev & prev == 0, arr.ind = T), "upper")
  h11 <- dedup(which(prev == prevprev & prev == 1, arr.ind = T), "upper")
  h01 <- dedup(which(prevprev < prev, arr.ind = T), "upper")
  h10 <- dedup(which(prevprev > prev, arr.ind = T), "upper")
  rands <- matrix(runif(N*N, 0, 1), nrow = N) # select random numbers from here
  
  # Modify the new adjacency matrix (upper triangle only)
  new[h00] <- ifelse(rands[h00] < add00, 1, 0)
  new[h11] <- ifelse(rands[h11] < lose11, 0, 1)
  new[h01] <- ifelse(rands[h01] < lose01, 0, 1)
  new[h10] <- ifelse(rands[h10] < add10, 1, 0)
  
  # Symmetrize the matrix
  new <- as.matrix(symmetrize(new, rule = "upper")) # copy the upper triangle over the lower triangle
  
  return(new)
}
