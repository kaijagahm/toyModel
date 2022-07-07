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
  rands <- matrix(runif(N*N, 0, 1), nrow = N) # select random numbers from here
  
  # Modify the new adjacency matrix (upper triangle only)
  new[h00] <- ifelse(rands[h00] < rbeta(1, shape1 = add00[1], shape2 = add00[2]), 1, 0)
  new[h11] <- ifelse(rands[h11] < rbeta(1, shape1 = lose11[1], shape2 = lose11[2]), 0, 1)
  new[h01] <- ifelse(rands[h01] < runif(1), 0, 1) # using uniform probabilities for 
  new[h10] <- ifelse(rands[h10] < runif(1), 1, 0)
  
  # Symmetrize the matrix
  new <- as.matrix(sna::symmetrize(new, rule = "upper")) # copy the upper triangle over the lower triangle
  
  return(new)
}


# remove.network.node -----------------------------------------------------
# Function to remove a node from the network.
remove.network.node <- function(network, previous, n.removed = 1, id = NULL, 
                                pm, # both bereaved 
                                ps, # one bereaved, one not
                                pa, # neither bereaved
                                histMultiplier) {
  # Calculate population size in the current network
  N <- nrow(network)
  
  # Assuming no node is specified, remove a random node. `del` is the node to remove.
  if(is.null(id)){
    del <- sample(1:N, n.removed, replace = FALSE)
  }else{
    del <- id
  }
  
  # First, capture the edges involving the removed individual
  edges <- network[del, -del] # row ([del,]), excluding self col ([,-del])
  bereaved <- which(edges == 1) # nodes that were connected to the removed individual
  non.bereaved <- which(edges == 0)

  # Remove the node
  network <- network[-del,] # rows 
  network <- network[,-del] # columns
  N <- N-1 # update the population size.
  
  # Now update the network.
  # First, allocate edges between mutually bereaved nodes (2nd-degree rewiring)
  # get potential edges, as edge list
  potentials <- which(network[bereaved, bereaved, # edges between second-degree connections
                              drop = FALSE] == 0, # keep format. Only edges that didn't already exist.
                      arr.ind = T)
  
  
  # then allocate a new edge vs not
  if(length(potentials) > 0){
    potentials <- dedup(potentials, triangle = "upper") # only the upper triangle

    # for each edge, decide whether it forms or not (0 or 1)
    probs <- abs(rnorm(nrow(potentials), mean = pm, sd = 0.1))
    multiply <- previous[bereaved, bereaved][potentials]*histMultiplier
    probs.adjusted <- probs*multiply
    probs.adjusted[probs.adjusted > 1] <- 1
    probs.adjusted[probs.adjusted < 0] <- 0
    new.edge <- rbinom(n = nrow(potentials), size = 1, prob = probs.adjusted)

    # update the network
    network[bereaved, bereaved][potentials] <- new.edge # update edges between bereaved with either a 0 or a 1
    network <- sna::symmetrize(network, rule = "upper") # copy upper triangle
  }
  
  # Second, allocate edges between bereaved and non-bereaved nodes
  
  potentials <- which(network[bereaved, non.bereaved, drop = FALSE] == 0, 
                      arr.ind = T)
  if(length(potentials) > 0){
    potentials <- dedup(potentials, triangle = "upper")
    # for each edge, decide whether it forms or not (0 or 1)
    probs <- abs(rnorm(nrow(potentials), mean = ps, sd = 0.1))
    multiply <- previous[bereaved, non.bereaved][potentials]*histMultiplier
    probs.adjusted <- probs*multiply
    probs.adjusted[probs.adjusted > 1] <- 1
    probs.adjusted[probs.adjusted < 0] <- 0
    new.edge <- rbinom(n = nrow(potentials), size = 1, prob = probs.adjusted)
    
    network[bereaved, non.bereaved][potentials] <- new.edge
    network <- sna::symmetrize(network, rule = "upper")
  }
  
  # Finally, allocate edges between mutually non-bereaved nodes
  potentials <- which(network[non.bereaved, non.bereaved, drop = FALSE] == 0, 
                      arr.ind = T)
  if(length(potentials) > 0){
    potentials <- dedup(potentials, triangle = "upper")
    
    # for each edge, decide whether it forms or not (0 or 1)
    probs <- abs(rnorm(nrow(potentials), mean = pa, sd = 0.1))
    multiply <- previous[non.bereaved, non.bereaved][potentials]*histMultiplier
    probs.adjusted <- probs*multiply
    probs.adjusted[probs.adjusted > 1] <- 1
    probs.adjusted[probs.adjusted < 0] <- 0
    new.edge <- rbinom(n = nrow(potentials), size = 1, prob = probs.adjusted)
    
    network[non.bereaved, non.bereaved][potentials] <- new.edge
    network <- sna::symmetrize(network, rule = "upper")
  }
  
  return(list(network = network, del = del))
}