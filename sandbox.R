# Parameters

# n = number of individuals
# m = number of individuals who die (mortality)
# remove = proportion of dyads removed per day
# create = proportion of dyads created per day
# Removal and connection should depend just on the previous day (maybe later we can allow memory beyond that)

library(igraph)
library(tidyverse)
n <- 30
density <- 0.15
seed <- 3
remove <- 0.1
create <- 0.05
days <- 10

set.seed(seed)
initialConnections <- sample(0:1, size = n*n, replace = T, prob = c(1-density, density))

# Create the initial adjacency matrix
am_init <- matrix(data = initialConnections, nrow = n, ncol = n)
init_graph <- graph_from_adjacency_matrix(am_init, weighted = NULL, mode = "undirected", diag = FALSE)

# Name nodes
V(init_graph)$name <- paste0("n", 1:n)

# Loss and gain of edges --------------------------------------------------
generateNetworks <- function(n, density, seed, remove, create, days){
  # Set the seed
  set.seed(seed)
  initialConnections <- sample(0:1, size = n*n, replace = T, prob = c(1-density, density))
  
  # Create the initial adjacency matrix
  am_init <- matrix(data = initialConnections, nrow = n, ncol = n)

  # For each day, remove `remove` connections and create `create` connections
  ams <- vector(mode = "list", length = days) # create a list to store the adjacency matrices
  ams[[1]] <- am_init
  for(i in 2:days){
    mat_prev <- ams[[i-1]]
    mat_new <- mat_prev
    for(col in 1:ncol(mat_prev)){
      for(row in 1:nrow(mat_prev)){
        if(mat_prev[row, col] == 0){
          mat_new[row, col] <- sample(0:1, size = 1, replace = TRUE, 
                                      prob = c(1-create, create))
        }else{
          mat_new[row, col] <- sample(0:1, size = 1, replace = TRUE, 
                                      prob = c(1-remove, remove))
        }
      }
    }
    ams[[i]] <- mat_new # assign new matrix to the corresponding place in the list.
  }
  
  # Make the adjacency matrices into graphs
  graphs <- lapply(ams, graph_from_adjacency_matrix, weighted = NULL, mode = "undirected", diag = FALSE)
  
  # Return both the adjacency matrices and the graphs
  return(list("ams" = ams, "graphs" = graphs))
}

# Create random coords for the nodes, which will be consistent when we plot
coords <- layout_with_fr(networks$graphs[[1]])

networks <- generateNetworks(n = 20, density = 0.1, seed = 3, remove = 0.1, create = 0.05, days = 10)

# Plot all the networks sequentially, individually
networkPlots <- lapply(networks$graphs, function(x){
  plot(x, layout = coords[,1:2], rescale=T)
})

# Plot all the networks as a facetted plot

cowplot::plot_grid()
