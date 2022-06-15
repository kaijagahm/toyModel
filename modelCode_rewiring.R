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
source("supportingFunctions.R") # all the functions that will be used in the main model

# remove.network.node -----------------------------------------------------
# Function to remove a node from the network.
remove.network.node <- function(network, n.removed = 1, id = NULL, pm) {
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
  
  # Remove the node
  network <- network[-del,] # rows 
  network <- network[,-del] # columns
  N <- N-1 # update the population size.
  
  # Now update the network:
  # First get potential edges, as edge list
  potentials <- which(network[bereaved, bereaved, # edges between second-degree connections
                              drop = FALSE] == 0, # keep format. Only edges that didn't already exist.
                      arr.ind = T)
  
  # then allocate a new edge vs not
  if(nrow(potentials) > 0){
    potentials <- dedup(potentials, triangle = "upper") # only the upper triangle
    
    # for each edge, decide whether it forms or not (0 or 1)
    new.edge <- sample(c(0,1), nrow(potentials), 
                       prob = c(1-pm, pm), replace = T)
    
    # update the network
    network[bereaved, bereaved][potentials] <- new.edge # update edges between bereaved with either a 0 or a 1
    network <- symmetrize(network, rule = "upper") # copy upper triangle
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
                         mode = "graph") # gives undirected graph, already symmetrized.

  # Run the baseline model
  ## set up the history
  network.history <- vector(mode = "list", length = burn.in)
  network.history[[1]] <- matrix(0, N, N) # blank network so we can look 2 timesteps back
  network.history[[2]] <- network.orig # baseline network
  
  ## run the loop, starting at index 3
  for(i in 3:burn.in){
    output <- update.network(ind = i, network.history)
    network.history[[i]] <- output # update history
  }
  
  # Save parameters as they stand for the end of the burn-in, before perturbation
  # assort.orig[zz,] <- assortment.continuous(network.orig, traits.orig, weighted=FALSE)$r
  # den.orig[zz,] <- gden(network.orig, mode="graph")
  # mean.deg.orig[zz,] <- mean(degree(network.orig, gmode="graph",ignore.eval=TRUE))
  # clust.orig[zz,] <- gtrans(network.orig,mode="graph")
  
  # Removal/perturbation and rewiring
  # Save final params
  # End
}







