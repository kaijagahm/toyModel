# Agent-based model

# Time steps; each time step has a certain number of interactions (i.e. 10). Probability distribution for how frequently they interact.
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(dils) # for filling in the edge list

# Define parametesr
tmax <- 10 # simulation time
n <- 40 # number of nodes
pint <- 0.4 # mean of the initial interaction probability distribution (base probability). How likely to interact.
iter <- 1
m <- 100 # number of interactions in each time step

# Start interaction loop
for(i in 1:iter){ # for each simulation
  
  # Initialize nodes to zeros state (done once at the beginning of the simulation)
  nodes <- data.frame(nodeID = 1:n, # assign each node an ID
              intProb = pint, # interaction probability
              timeFromLastInt = 0, # number of time steps since node's last encounter with another node
              nInt = 0) # number of interactions this node had
  
  # Initialize lists to store the adjacency matrices and graphs
  ams <- vector(mode = "list", length = day)
  gs <- vector(mode = "list", length = day)
  
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
    gs[[day]] <- graph_from_adjacency_matrix(am, mode = "undirected", diag = FALSE)
    
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
} # close simulation iteration

gs <- lapply(gs, function(x){
  x %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(name = 1:n)
})

# Create coordinates to use for plotting based on the optimal layout on the first day.
layoutCoords <- layout_with_fr(gs[[1]])

# Make a bunch of plots with the same layout
lapply(gs, function(x){
  x %>% ggraph(layout = layoutCoords)+
    geom_edge_link(edge_width = 0.2)+
    geom_node_point(col = "steelblue", size = 5)+
    geom_node_text(aes(label = name, vjust = 0.5), col = "black")
  })

# Calculate node-level stats 
stats <- lapply(gs, function(x){
  x %>% 
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(centr = centrality_degree(),
           deg = degree(.))
})
