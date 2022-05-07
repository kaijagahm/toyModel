# Agent-based model

# Time steps; each time step has a certain number of interactions (i.e. 10). Probability distribution for how frequently they interact.
library(tidyverse)
library(igraph)
library(ggraph)
library(dils) # for filling in the edge list
library(tidygraph)
library(data.table)

# Define parametesr
tmax <- 10 # simulation time
n <- 60 # number of nodes
pint <- 0.4 # mean of the initial interaction probability distribution (base probability). How likely to interact.
iter <- 1

# Start interaction loop
runSim <- function(tmax = 10, n = 50, pint = 0.4){
  
  # ***There are going to be three outputs: nodes, ams, gs. ***
  # Initialize nodes to zeros state (done once at the beginning of the simulation)
  # 1. NODES ***
  nodes <- data.frame(nodeID = 1:n, # assign each node an ID
                      intProb = pint, # interaction probability
                      timeFromLastInt = 0, # number of time steps since node's last encounter with another node
                      nInt = 0) # number of interactions this node had
  
  # Initialize lists to store the graphs, interactions, and node-level data
  # 1. Graphs
  gs <- vector(mode = "list", length = tmax) # graphs
  
  # 2. Interactions
  ints <- vector(mode = "list", length = tmax) # interaction edge lists
  
  # 3. Node-level data
  nodeData <- vector(mode = "list", length = tmax) # node data
  
  # Run the simulation for the number of days
  for(day in 1:tmax){
    cat(paste("day", day, "\n"))
    
    # Randomize which nodes interact on this day
    # empty adjacency matrix to fill with interactions
    interactions <- matrix(nrow = 0, ncol = 2) 
    
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
    
    # We don't want to allow any isolated nodes. If there's a node that isn't connected to anyone, add one random edge.
    for(node in 1:n){
      if(!(node %in% interactions[,1]) & !(node %in% interactions[,2])){
        allNodes <- 1:n # all nodes
        others <- allNodes[allNodes != node] # remove self from vector of possibilities
        toAdd <- sample(others, 1) # pick one node to join the unconnected node to
        interactions <- rbind(interactions, c(node, toAdd)) # add to edge list
      }
    }
    
    # Convert interactions to a full adjacency matrix, including any nodes that didn't interact
    g <- graph_from_data_frame(interactions, vertices = nodes[,"nodeID"], directed = F) %>%
      as_tbl_graph() %>%
      activate(nodes) %>%
      mutate(name = 1:n) # name the nodes
    
    # Calculate per-node stats for this day of interactions
    interactingNodes <- as.data.frame(interactions) %>%
      pivot_longer(cols = 1:2) %>%
      pull() %>% as.integer() %>% sort() # get non-unique vector of all nodes that interacted
    
    for(r in 1:nrow(nodes)){
      if(nodes[r, "nodeID"] %in% interactingNodes){
        nodes[r, "nInt"] <- nodes[r, "nInt"] + sum(interactingNodes == r) # add number of interactions from today
      }
      if(nodes[r, "nodeID"] %in% interactingNodes){
        nodes[r, "timeFromLastInt"] <- 0
      }else{
        nodes[r, "timeFromLastInt"] <- nodes[r, "timeFromLastInt"] + 1
      }
    }
    
    # Save to list
    gs[[day]] <- g
    ints[[day]] <- interactions
    nodeData[[day]] <- nodes
    
  } # close day
  
  outputList <- list("gs" = gs, "interactions" = ints, "nodeData" = nodeData)
  return(outputList)
}

sim <- runSim(tmax = 10, n = 60, pint = 0.4)

# Create coordinates to use for plotting based on the optimal layout on the first day.
layoutCoords <- layout_with_fr(sim$gs[[1]])

# Make a bunch of plots with the same layout
lapply(sim$gs, function(x){
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

stats <- getNodeStats(sim$gs, type = "graphs")
nodeData <- getNodeStats(sim$gs, type = "list")
nodeDataDF <- getNodeStats(sim$gs, type = "df")

# Plot degree distributions
# Now use this to plot
nodeDataDF %>%
  ggplot(aes(x = deg, group = Day, color = Day))+
  geom_density()+
  theme_minimal()+
  theme(axis.title.y = element_blank())+
  xlab("Degree")
# The probability of association doesn't change day to day, so it makes sense that all our distributions are basically the same.


# How does mean degree change depending on interaction prob? --------------
probsToTest <- seq(from = 0.01, to = 0.99, by = 0.01)
simOutputs <- vector(mode = "list", length = length(probsToTest))

for(i in 1:length(probsToTest)){
  sim <- runSim(tmax = 10, n = 60, pint = probsToTest[i])
  simOutputs[[i]] <- sim
}

# Get all the graphs
gsList <- lapply(simOutputs, function(x){
  x[["gs"]]
})

# Compute the node stats
nodeDataDFs <- lapply(gsList, function(x){
  getNodeStats(x, type = "df")
})

nodeDataAllSims <- map2(.x = nodeDataDFs, .y = probsToTest, 
                        .f = function(.x, .y){
                          .x %>% mutate(intProb = .y)
                        }) %>%
  data.table::rbindlist() %>% as.data.frame()

# Now we can finally plot
# Degree
nodeDataAllSims %>%
  filter(Day == 1) %>% # we only need the first day's data because all days are the same
  ggplot(aes(x = intProb, y = deg))+
  geom_point(size = 0.6, alpha = 0.2)+
  geom_smooth()+
  theme_minimal()+
  ylab("Degree")+ # high degree means the node is connected to many other nodes
  xlab("Interaction probability")+
  theme(text = element_text(size = 20))

# Centrality
nodeDataAllSims %>%
  filter(Day == 1) %>% # we only need the first day's data because all days are the same
  ggplot(aes(x = intProb, y = centr))+
  geom_point(size = 0.6, alpha = 0.2)+
  geom_smooth()+
  theme_minimal()+
  ylab("Eigenvector centrality")+ # high EC means the node is connected to many well-connected nodes
  xlab("Interaction probability")+
  theme(text = element_text(size = 20))

# Network edge density
densityList <- lapply(gsList, function(x){
  lapply(x, edge_density) %>% 
    unlist() %>% 
    as.data.frame() %>%
    setNames("density")
})

densities <- densityList %>% 
  map2(., probsToTest, function(.x, .y){
    .x %>%
      mutate(intProb = .y)
  }) %>%
  data.table::rbindlist() %>%
  as.data.frame()

densities %>%
  ggplot(aes(x = intProb, y = density))+
  geom_point(size = 0.6, alpha = 0.2)+
  geom_smooth()+
  theme_minimal()+
  ylab("Network Density")+ 
  xlab("Interaction probability")+
  theme(text = element_text(size = 20))

