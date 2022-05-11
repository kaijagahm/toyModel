# Agent-based model

# Time steps; each time step has a certain number of interactions (i.e. 10). Probability distribution for how frequently they interact.
library(tidyverse)
library(igraph)
library(ggraph)
library(dils) # for filling in the edge list
library(tidygraph)
library(data.table)

# Function to run the model
runSim <- function(tmax = 10, n = 50, pint = 0.4, pNew = 0.1, pLose = 0.1){
  
  # ***There are going to be three outputs: nodes, ams, gs. ***
  # Initialize nodes to zeros state (done once at the beginning of the simulation)
  # 1. NODES ***
  nodes <- data.frame(nodeID = 1:n, # assign each node an ID
                      intProb = pint, # interaction probability
                      timeFromLastInt = 0, # number of time steps since node's last encounter with another node
                      nInt = 0) # number of interactions this node had
  
  # Initialize lists to store the graphs, interactions
  # 1. Graphs
  gs <- vector(mode = "list", length = tmax) # graphs
  
  # 2. Interactions
  ints <- vector(mode = "list", length = tmax) # interaction edge lists
  
  # Initialize the graph for the first day
  # Randomize which nodes interact on this day
  interactions <- matrix(nrow = 0, ncol = 2) # empty adjacency matrix to fill with interactions
  # Interact each node with each other node
  for(nodeA in 1:n){ 
    for(nodeB in (nodeA+1):n){ # don't need to do the dyads twice; hence nodeA+1
      # Ignore impossible values of nodeB
      if(nodeB > n) next
      if(nodeB == nodeA) next # no self edges!
      # Calculate the probability of the two nodes meeting
      probMeeting <- nodes[nodeA, "intProb"]*nodes[nodeB, "intProb"]
      if(runif(1) < probMeeting){ # random draw between 0 and 1
        # make an edge list
        interactions <- rbind(interactions, c(nodeA, nodeB)) # the nodes and their weight
      }
    }
  }
  
  # We don't want to allow any isolated nodes. 
  # If there's a node that isn't connected to anyone, add one random edge.
  for(node in 1:n){
    if(!(node %in% interactions[,1]) & !(node %in% interactions[,2])){
      allNodes <- 1:n # all nodes
      others <- allNodes[allNodes != node] # remove self from vector of possibilities
      toAdd <- sample(others, 1) # pick one node to join the unconnected node to
      interactions <- rbind(interactions, c(node, toAdd)) # add to edge list
    }
  }
  
  # Convert interactions to a graph
  g <- graph_from_data_frame(interactions, vertices = nodes[,"nodeID"], directed = F) %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(name = 1:n) # name the nodes
  
  # Get initial adjacency matrix
  amDay1 <- get.adjacency(g) %>% as.matrix()
  
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
    
    # Save this day's adjacency matrix
    ams[[day]] <- newAM
  } # close day
  
  # Make each of the edge lists into a graph, using only the upper triangle
  gs <- lapply(ams, function(x){
    graph_from_adjacency_matrix(x, mode = "upper", diag = FALSE)
  })
  
  return(gs)
}

sim60 <- runSim(tmax = 10, n = 60, pint = 0.4, pNew = 0.3, pLose = 0.5)
sim50 <- runSim(tmax = 10, n = 50, pint = 0.4, pNew = 0.3, pLose = 0.5)
sim40 <- runSim(tmax = 10, n = 40, pint = 0.4, pNew = 0.3, pLose = 0.5)
sim30 <- runSim(tmax = 10, n = 30, pint = 0.4, pNew = 0.3, pLose = 0.5)

# Create coordinates to use for plotting based on the optimal layout on the first day.
layoutCoords <- layout_with_fr(sim[[1]])

# Make a bunch of plots with the same layout
lapply(sim, function(x){
  x %>% ggraph(layout = layoutCoords)+
    geom_edge_link(edge_width = 0.2)+
    geom_node_point(col = "steelblue", size = 5)+
    geom_node_text(aes(label = name, vjust = 0.5), col = "black")
})

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
