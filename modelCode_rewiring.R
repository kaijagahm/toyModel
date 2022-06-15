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







