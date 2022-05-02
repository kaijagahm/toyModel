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
m <- 100 # number of interactions in each time step

# Start interaction loop
runSim <- function(iter, tmax, n, pint, m){
  outputs <- vector(mode = "list", length = iter) # store stuff
  
  for(i in 1:iter){ # for each simulation
    # ***There are going to be three outputs: nodes, ams, gs. ***
    # Initialize nodes to zeros state (done once at the beginning of the simulation)
    # 1. NODES ***
    nodes <- data.frame(nodeID = 1:n, # assign each node an ID
                        intProb = pint, # interaction probability
                        timeFromLastInt = 0, # number of time steps since node's last encounter with another node
                        nInt = 0) # number of interactions this node had
    
    # Initialize lists to store the adjacency matrices and graphs
    # 2. AMS ***
    ams <- vector(mode = "list", length = tmax) # adjacency matrices
    # 3. GS *** 
    gs <- vector(mode = "list", length = tmax) # graphs
    
    # Run the simulation for the number of days
    for(day in 1:tmax){
      cat(paste("day", day, "\n"))
      
      # Randomize which nodes interact on this day
      # empty adjacency matrix to fill with interactions
      am <- matrix(data = 0, nrow = n, ncol = n) 
      interactions <- matrix(nrow = 0, ncol = 3) 
      
      # Interact each node with each other node
      for(nodeA in 1:n){ 
        for(nodeB in (nodeA+1):n){ # don't need to do the dyads twice; hence nodeA+1
          # Ignore impossible values of nodeB
          if(nodeB > n) next
          # Calculate the probability of the two nodes meeting
          probMeeting <- nodes[nodeA, "intProb"]*nodes[nodeB, "intProb"]
          if(runif(1) < probMeeting){ # random draw between 0 and 1
            am[nodeA, nodeB] <- am[nodeA, nodeB] + 1 # increment the adjacency matrix
            # also make an edge list
            interactions <- rbind(interactions, c(nodeA, nodeB, 1)) # the nodes and their weight
          }
        }
      }
      
      # Check that we haven't exceeded the number of allowed interactions per day
      # If we have, remove some interactions at random.
      if(nrow(interactions) > m){
        # how many interactions do we need to remove?
        howManyToRemove <- nrow(interactions) - m
        # select random rows to remove
        whichToRemove <- sample(1:nrow(interactions),
                                size = howManyToRemove)
        
        # remove from adjacency matrix first
        for(k in 1:length(whichToRemove)){
          if(am[interactions[k,1], interactions[k,2]] > 0){
            am[interactions[k,1], interactions[k,2]] <- am[interactions[k,1], interactions[k,2]] - 1
          }
        }
        # ...and now remove the rows from the edge list
        interactions <- interactions[-whichToRemove,]
      }
      
      # We don't want to allow any isolated nodes. If there's a node that isn't connected to anyone, add one random edge.
      for(node in 1:n){
        if(rowSums(am)[node] == 0){
          allNodes <- 1:n # all nodes
          others <- allNodes[allNodes != node] # remove self from vector of possibilities
          toAdd <- sample(others, 1) # pick one node to join the unconnected node to
          am[node, toAdd] <- 1 # add to adjacency matrix
          interactions <- rbind(interactions, c(node, toAdd, 1)) # add to edge list
        }
      }
      
      # Save adjacency matrix
      ams[[day]] <- am
      # Name nodes
      gs <- graph_from_adjacency_matrix(am, mode = "undirected", diag = FALSE) %>%
        as_tbl_graph() %>%
        activate(nodes) %>%
        mutate(name = 1:n)
      # Save to list
      gs[[day]] <- gs
      
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
    } # close day
    
    # Save the outputs for this iteration of the simulation
    outputs[[i]] <- list(nodes, ams, gs)
  } # close simulation iteration
  
  return(outputs) # a massive list of lists
}

sim <- runSim(iter = 1, tmax = 10, n = 60, pint = 0.4, m = 100)


# Create coordinates to use for plotting based on the optimal layout on the first day.
layoutCoords <- layout_with_fr(gs[[1]])

# Make a bunch of plots with the same layout
lapply(gs, function(x){
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
      mutate(centr = centrality_degree(),
             deg = degree(.))
  })
  
  # Extract just the node data, making a list of data frames
  statsDFList <- lapply(stats, function(x){
    x %>% activate(nodes) %>% as.data.frame()
  })
  
  # Compress the list of data frames into a single df
  statsDF <- data.table::rbindlist(statsDFList, idcol = "Day")
  
  # Return different things based on what the user wants
  if(type == "df"){
    return(statsDF)
  }else if(type == "list"){
    return(statsDFList)
  }else if(type == "graphs"){
    return(stats)
  }
}

collapseAM <- function(ams){
  weighted <- Reduce("+", ams) %>% 
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(name = 1:n) %>%
    mutate(deg = degree(.)) %>%
    mutate(str = strength(.))
  
  return(weighted)
}

stats <- getNodeStats(gs, type = "graphs")
nodeData <- getNodeStats(gs, type = "list")
nodeDataDF <- getNodeStats(gs, type = "df")
weighted <- collapseAM(ams)

# Plot the weighted adjacency matrix
weighted %>%
  ggraph(layout = "auto")+
  geom_edge_link(aes(width = weight))+
  geom_node_point(aes(col = deg), size = 5)+
  scale_edge_width(range = c(0.1, 1))+
  geom_node_text(aes(label = name), col = "white", size = 2)+
  theme_graph()

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
probsToTest <- c(0.1, 0.3, 0.5, 0.7)


