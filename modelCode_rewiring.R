# Model code, version 2
# Architecture of this version is based more directly on Farine's 2nd-degree rewiring model (2021)
# AUTHOR: Kaija Gahm (kgahm@ucla.edu)
# REVISION DATE: 2022-07-06 (added parameters derived in parameterizingTheModel.Rmd)

# Load packages -----------------------------------------------------------
library(tidyverse)
library(igraph)
library(ggraph)
library(gganimate)
library(tidygraph)
library(data.table)
library(checkmate)
library(graphlayouts)
library(patchwork)
library(sna) # for symmetrizing matrices and other social network analysis
source("supportingFunctions.R") # all the functions that will be used in the main model

# Section 2: Generate networks and simulate node loss --------------------
# Network parameters
N <- 50 # Nodes in the network
nodes.removed <- 1 # Nodes to remove
n.removed <- 1 # How many to remove at a time
edge.prob <- 0.04 # This is derived from the average network density, when taking a 5-day increment, from parameterizingTheModel.Rmd.

# Simulation parameters
n.rep <- 100 # number of repetitions
burn.in <- 50 # number of days to burn in
burn.out <- 50 # number of days to continue after the perturbation
pm <- 0.3
ps <- 0.1
pa <- 0.2
add00 <- c(0.4721719, 7.3144796) # beta distribution parameters derived from parameterizingTheModel.Rmd.
lose01 <- 0.3 
add10 <- 0.2
lose11 <- c(0.3283134, 0.3062181) # beta distribution parameters derived from parameterizingTheModel.Rmd.
histMultiplier <- 1.2

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
    output <- update.network(ind = i, network.history, add00 = add00, 
                             add10 = add10, lose01 = lose01, lose11 = lose11)
    network.history[[i]] <- output # update history
  }
  
  # Removal/perturbation and rewiring
  rewired.list <- remove.network.node(network = network.history[[burn.in]], 
                                      previous = network.history[[burn.in-1]],
                                      histMultiplier = histMultiplier,
                                      n.removed = 1, pm = rnorm(1, pm), ps = rnorm(1, ps), pa = rnorm(1, pa))
  rewired.network <- rewired.list$network # the actual rewired network
  rewired.del <- rewired.list$del # removed individuals. Used for amending history matrices in order to continue baseline dynamics moving forward
  
  # Continue baseline dynamics following removal
  revised.history <- lapply(network.history, function(x){
    x <- x[-rewired.del,]
    x <- x[,-rewired.del]
    return(x)
  })
  
  revised.history <- append(revised.history, 
                            list(rewired.network))
  revised.history <- append(revised.history, rep(NA, burn.out))
  
  for(i in (burn.in+2):length(revised.history)){
    output <- update.network(ind = i, revised.history, add00 = add00, 
                             add10 = add10, lose01 = lose01, lose11 = lose11)
    revised.history[[i]] <- output # update history
  }
}

# Create a composite list of all the adjacency matrices, using the full networks before loss
full.history <- c(network.history, 
                  revised.history[(burn.in+1):length(revised.history)])

# Make these into networks
full.history.graphs <- lapply(full.history, function(x){
  graph_from_adjacency_matrix(x, mode = "undirected", diag = FALSE)
})


# Animate -----------------------------------------------------------------
# Code adapted from https://www.r-bloggers.com/2021/09/animating-network-evolutions-with-gganimate/.

# Create single layout
xy <- layout_nicely(full.history.graphs[[1]]) # all nodes before loss

pList <- vector("list", length(full.history.graphs))
for (i in 1:length(full.history.graphs)) {
  g <- full.history.graphs[[i]]
  pList[[i]] <- ggraph(g, layout = "manual", x = xy[1:length(g), 1], 
                       y = xy[1:length(g), 2]) +
    geom_edge_link0(edge_width = 0.4, alpha = 0.2) +
    geom_node_point(shape = 19, size = 6) +
    theme_graph() +
    theme(legend.position = "bottom")
}

library(animation)
saveGIF({
  for (i in 1:length(pList)) plot(pList[[i]])
}, movie.name = "networkSim.gif", interval = 0.1)





