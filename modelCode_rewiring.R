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
burn.out <- 50 # number of days to continue after the perturbation
pm <- 0.3
ps <- 0.1
pa <- 0.2

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
  network.history <- vector(mode = "list", length = burn.in+1+burn.out)
  network.history[[1]] <- matrix(0, N, N) # blank network so we can look 2 timesteps back
  network.history[[2]] <- network.orig # baseline network
  
  ## run the loop, starting at index 3
  for(i in 3:burn.in){
    output <- update.network(ind = i, network.history)
    network.history[[i]] <- output # update history
  }
  
  # Removal/perturbation and rewiring
  rewired <- remove.network.node(network = network.history[[burn.in]], 
                                 n.removed = 1, pm = pm, ps = ps, pa = pa)
  network.history[[burn.in+1]] <- rewired
  
  # XXX can't do this because the networks are now different sizes.
  # # Continue baseline dynamics following removal
  # ## run the loop, starting at index burn.in + 2
  # for(i in burn.in+2:length(network.history)){
  #   output <- update.network(ind = i, network.history)
  #   network.history[[i]] <- output # update history
  # }
  # # End
}







