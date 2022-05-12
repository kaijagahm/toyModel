# Agent-based model

# Time steps; each time step has a certain number of interactions (i.e. 10). Probability distribution for how frequently they interact.
library(tidyverse)
library(igraph)
library(ggraph)
library(dils) # for filling in the edge list
library(tidygraph)
library(data.table)

# XXX START HERE: next step is to include some sort of memory for which nodes were connected to which other ones in the previous time step, and have the rewiring incorporate that.
# I'm just struggling to work with the different data structures for networks.

connectIsolatedUpper <- function(mat){
  for(r in 1:(nrow(mat)-1)){
    for(c in 2:ncol(mat)){
      if(all(mat[r, c:ncol(mat)] == 0)){
        rand <- sample(c:ncol(mat), size = 1)
        mat[r, rand] <- 1
      }
    }
  }
  return(mat)
}

# MODEL:
# Returns a list of igraph network objects (graphs), one for each day in 1:tmax. 
# Graphs are BINARY and UNDIRECTED, based on the UPPER TRIANGLE of the adjacency matrix. 
runSim <- function(tmax = 10, # length of time over which to run the simulation
                   n = 50, # number of individuals
                   pNew = 0.1, # probability of gaining an edge given that it doesn't exist
                   pLose = 0.1, # probability of losing an edge given that it exists
                   allowIsolated = TRUE){ # allow nodes to be unconnected? If F, randomly connects each unconnected node with one other node.
  
  # STORAGE
  # Initialize lists to store graphs and ajacency matrices for each time step
  gs <- vector(mode = "list", length = tmax) # storage for graphs for each day--this will be the output.
  ams <- vector(mode = "list", length = tmax) # storage for adjacency matrices for each day
  
  # SETUP
  # Initialize the adjacency matrix for the first day
  amDay1 <- matrix(sample(0:1, n*n, replace = TRUE, prob = c(1-pNew, pNew)), n, n)
  if(allowIsolated == FALSE){
    amDay1 <- connectIsolatedUpper(amDay1)
  }
  ams[[1]] <- amDay1
  
  # RUN SIMULATION
  # Loop through the days
  for(day in 2:tmax){
    cat(paste("day", day, "\n")) # minimal feedback while running the model, for sanity check on speed
    previousAM <- ams[[day-1]] # previous day's matrix to operate on
    
    # Create this day's adjacency matrix
    newAM <- previousAM
    
    # For each edge, determine whether it is created, destroyed, or left alone.
    newAM <- apply(newAM, c(1,2), function(x, new = pNew, lose = pLose){
      switch <- runif(1)
      if(x == 0 & switch < new){
        return(1)
      }else if(x == 0 & switch >= new){
        return(0)
      }else if(x == 1 & switch < lose){
        return(0)
      }else if(x == 1 & switch >= lose){
        return(1)
      }
    })
    
    # If we're not allowing isolated nodes, randomly connect each node one time.
    if(allowIsolated == FALSE){
      # Add a random edge for any isolated nodes, only paying attention to the upper triangle
      newAM <- connectIsolatedUpper(newAM)
    }
    
    # Save this day's adjacency matrix
    ams[[day]] <- newAM
  } # close day
  
  # Make each of the adjacency matrices into a graph, using only the upper triangle
  gs <- lapply(ams, function(x){
    graph_from_adjacency_matrix(x, mode = "upper", diag = FALSE)
  })
  
  return(gs)
}

sim60 <- runSim(tmax = 10, n = 60, pNew = 0.1, pLose = 0.9, allowIsolated = T)
sim50 <- runSim(tmax = 10, n = 50, pNew = 0.1, pLose = 0.9, allowIsolated = T)
sim40 <- runSim(tmax = 10, n = 40, pNew = 0.1, pLose = 0.9, allowIsolated = T)
sim30 <- runSim(tmax = 10, n = 30, pNew = 0.1, pLose = 0.9, allowIsolated = T)
sim10 <- runSim(tmax = 10, n = 10, pNew = 0.1, pLose = 0.9, allowIsolated = T)

# Create coordinates to use for plotting based on the optimal layout on the first day.
layoutCoords <- layout_with_fr(sim10[[1]])

# Make a bunch of plots with the same layout
lapply(sim10, function(x){
  x %>% ggraph(layout = layoutCoords)+
    geom_edge_link(edge_width = 0.2)+
    geom_node_point(size = 5)})

getNodeStats <- function(gs, type = "df"){
  # Check to make sure the "type" argument is valid
  if(!(type %in% c("df", "graphs", "list"))){
    stop("Argument 'type' must be 'df', 'graphs', or 'list'.")
  }
  
  # Calculate stats
  stats <- lapply(gs, function(x){
    x %>% 
      as_tbl_graph() %>%
      activate(nodes) %>%
      mutate(centr = centrality_eigen(),
             deg = degree(.))
  })
  
  # Extract just the node data, making a list of data frames
  statsDFList <- lapply(stats, function(x){
    x %>% activate(nodes) %>% 
      as.data.frame()
  })
  
  # Compress the list of data frames into a single df
  statsDF <- data.table::rbindlist(statsDFList, idcol = "Day") %>% 
    as.data.frame()
  
  # Return different things based on what the user wants
  if(type == "df"){
    return(statsDF)
  }else if(type == "list"){
    return(statsDFList)
  }else if(type == "graphs"){
    return(stats)
  }
}

# Run the simulation with various values of n
nodeDataDF60 <- getNodeStats(sim60, type = "df") %>% mutate(n = 60)
nodeDataDF50 <- getNodeStats(sim50, type = "df") %>% mutate(n = 50)
nodeDataDF40 <- getNodeStats(sim40, type = "df") %>% mutate(n = 40)
nodeDataDF30 <- getNodeStats(sim30, type = "df") %>% mutate(n = 30)
nd <- rbind(nodeDataDF60, nodeDataDF50, nodeDataDF40, nodeDataDF30)

# How does the average degree change over time, with different n?
nd %>%
  group_by(n, Day) %>%
  summarize(mnDeg = mean(deg)) %>%
  ggplot(aes(x = Day, y = mnDeg, col = n, group = n))+
  geom_point()+
  geom_line() # saturates a bit above 1/3 n
