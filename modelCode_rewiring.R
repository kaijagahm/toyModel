# Model code, version 2
# Architecture of this version is based more directly on Farine's 2nd-degree rewiring model (2021)
# AUTHOR: Kaija Gahm (kgahm@ucla.edu)
# REVISION DATE: 2022-05-14

# Code overview -----------------------------------------------------------
# 1. random starting network
# 2. burn-in period. This is where the baseline dynamics come in. For Farine, the baseline dynamics were the Ilany & Akcay model of social inheritance. For me, the equivalent burn-in will be the random network fluctuations, without adding or removing any nodes. Edges get removed and added based on probabilities.

# Load packages -----------------------------------------------------------
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(data.table)
library(checkmate)
library(sna) # for symmetrizing matrices and other social network analysis

# Section 1: Define functions ---------------------------------------------
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
update.network <- function(ind, history, 
                           # p gain an edge given not connected in either of the 
                           # previous two time steps
                           add00 = 0.1, 
                           # p lose an edge given connected in prev time step but 
                           # not prev prev
                           lose01 = 0.3, 
                           # p gain an edge given connected in prev prev time step 
                           # but not prev
                           add10 = 0.3, 
                           # p lose edge given connected in previous 2 time steps
                           lose11 = 0.1){ 
  # get history two steps back
  prev <- history[[ind-1]]
  prevprev <- history[[ind-2]]
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
  new <- symmetrize(new, rule = "upper") # copy the upper triangle over the lower triangle
  
  
  
  # determine fate of each edge: operate on all the edges, and then disregard one half of the triangle
  
  # update network
  network <- network # and do some stuff involving edgeInfo.
  
  return(list(network = network, edgeInfo = edgeInfo))
  
}

# remove.network.node -----------------------------------------------------
# Function to remove a node from the network.
## n.removed tells how many nodes to remove. Default is 1. Later, expand this to multiple nodes and define them according to a probability density function.
## id is the id of the node to remove. Default is NULL --> remove a random node.
## network is the starting network, to be perturbed.
remove.network.node <- function(network, n.removed = 1, id = NULL) {
  # Check that `network` is an adjacency matrix.
  
  # Calculate pop size
  N <- length(V(network))
  
  # Assuming no node is specified, remove a random node. `del` is the node to remove.
  if(is.null(id)){
    del <- sample(1:N, n.removed, replace = FALSE)
  }else{
    del <- id
  }
  
  # First, capture the edges involving the removed individual
  edges <- network[del, -del]
  
  # Remove the node
  network <- network[-del,] # rows
  network <- network[,-del] # columns
  N <- N-1 # update population size.
  
  # Now update the network:
  # First get potential edges
  potentials <- which(network[which(edges==1), which(edges==1), # edges between second-degree connections
                              drop = FALSE] == 0, # keep format. Only edges that didn't exist.
                      arr.ind = T) # get row and column indices
  
  # then allocate a new edge vs not
  if(length(potentials) > 0){
    potentials <- potentials[which(potentials[,1] < potentials[,2]),,drop=FALSE]  # to avoid duplicates
    new.edge <- sample(c(0,1),nrow(potentials),prob=c(1-pm,pm),replace=T)
    network[which(edges==1),which(edges==1)][cbind(potentials[,1],potentials[,2])] <- new.edge
    network[which(edges==1),which(edges==1)][cbind(potentials[,2],potentials[,1])] <- new.edge
  }
  
  # then randomly allocate edges between newly disconnected nodes and other nodes
  potentials <- which(network[which(edges==1),which(edges==0),drop=FALSE]==0,arr.ind=T)
  if (length(potentials) > 0) {
    new.edge <- sample(c(0,1),nrow(potentials),prob=c(1-ps,ps),replace=T)
    network[which(edges==1),which(edges==0)][cbind(potentials[,1],potentials[,2])] <- new.edge
    network[which(edges==0),which(edges==1)][cbind(potentials[,2],potentials[,1])] <- new.edge
  }
  
  # return
  return(network)
}

# Section 2: Generate networks and simulate node loss --------------------
# Network parameters
N <- 50 # Nodes in the network
nodes.removed <- 1 # Nodes to remove
n.removed <- 1 # How many to remove at a time
edge.prob <- 0.1 # initial probability of edges in the random starting network

# Simulation parameters
n.rep <- 100 # number of repetitions
burn.in <- 50 # number of days to burn in

# Results storage
den.orig <- matrix(NA, n.rep)
den.rewired <- matrix(NA, n.rep)
mean.deg.orig <- matrix(NA, n.rep)
mean.deg.rewired <- matrix(NA, n.rep)
# assort
# clust

for(zz in 1:n.rep){
  # Generate a random starting network
  network.orig <- rgraph(N, tprob = edge.prob, 
                         mode = "graph") # gives undirected graph
  # XXX connect unconnected nodes
  
  # Run the baseline model
  history <- vector(mode = "list", length = burn.in)
  history[[1]] <- matrix(0, N, N) 
  history[[2]] <- matrix(0, N, N)
  
  for(i in 3:burn.in){
    output <- update.network(ind = i, history = history)
    network.orig <- output
  }
  
  # Save original params
  assort.orig[zz,] <- assortment.continuous(network.orig, traits.orig, weighted=FALSE)$r
  den.orig[zz,] <- gden(network.orig, mode="graph")
  mean.deg.orig[zz,] <- mean(degree(network.orig, gmode="graph",ignore.eval=TRUE))
  clust.orig[zz,] <- gtrans(network.orig,mode="graph")
  # Do a removal and rewire
  # Save final params
  # End
}







