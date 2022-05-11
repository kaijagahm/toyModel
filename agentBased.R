# Agent-based model

# Time steps; each time step has a certain number of interactions (i.e. 10). Probability distribution for how frequently they interact.
library(tidyverse)
library(igraph)
library(ggraph)
library(dils) # for filling in the edge list
library(tidygraph)
library(data.table)

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

# Function to run the model
runSim <- function(tmax = 10, n = 50, pNew = 0.1, pLose = 0.1, allowIsolated = TRUE){
  
  # Initialize adjacency matrix
  
  # Initialize lists to store the graphs, interactions
  # 1. Graphs
  gs <- vector(mode = "list", length = tmax) # graphs
  
  # Initialize the graph for the first day
  amDay1 <- matrix(sample(0:1, n*n, replace = TRUE, prob = c(1-pNew, pNew)), n, n)
  if(allowIsolated == FALSE){
    amDay1 <- connectIsolatedUpper(amDay1)
  }
  
  # Okay, this is our initial adjacency matrix.
  # Run the simulation for the number of days
  # idea: if an edge is lost, resulting in the degree of that node dropping to zero, then some high chance of creating a new edge for that node, selecting from any other nodes that its previous partner is connected to.
  ams <- vector(mode = "list", length = tmax) # storage for edgelists for each day
  ams[[1]] <- amDay1
  
  for(day in 2:tmax){
    cat(paste("day", day, "\n"))
    previousAM <- ams[[day-1]] # here's what we're working with
    
    # Create this day's edge list
    newAM <- previousAM
    
    switchfun <- function(x){
      switch <- runif(1)
      if(x == 0 & switch < pNew){
        return(1)
      }else if(x == 0 & switch >= pNew){
        return(0)
      }else if(x == 1 & switch < pLose){
        return(0)
      }else if(x == 1 & switch >= pLose){
        return(1)
      }
    }
    newAM <- apply(newAM, c(1,2), switchfun)
    
    if(allowIsolated == FALSE){
      # Add a random edge for any isolated nodes, only paying attention to the upper triangle
      newAM <- connectIsolatedUpper(newAM)
    }
    
    # Save this day's adjacency matrix
    ams[[day]] <- newAM
  } # close day
  
  # Make each of the edge lists into a graph, using only the upper triangle
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
