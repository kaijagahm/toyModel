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

# Section 1: Define functions ---------------------------------------------
# Function for baseline network dynamics
# Repeat this number of times specified for desired burn.in, in a for loop. Each time, spitting out the network, and the history of the edges.
update.network <- function(network, ) { 
  # assign edge histories
  edgeInfo <- NULL
  
  # determine fate of each edge
  edgeInfo <- NULL
  
  # update network
  network <- network # and do some stuff involving edgeInfo.
  
  return(list(network = network, edgeInfo = edgeInfo))
  
}

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





