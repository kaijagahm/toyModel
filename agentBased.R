# Agent-based model

# Time steps; each time step has a certain number of interactions (i.e. 10). Probability distribution for how frequently they interact.
library(tidyverse)
library(igraph)

# Define parametesr
tmax <- 10 # simulation time
n <- 100 # number of nodes
pint <- 0.03 # mean of the initial interaction probability distribution (base probability). How likely to interact.
iter <- 1
m <- 10 # number of interactions in each time step

# Initialize bookkeeping
intMat <- rep(list(matrix(data = 0, nrow = n, ncol = n)),
              iter) # an interaction matrix for each iteration of the loop
# Each element of the intMat is a single iteration. For each element, rows and columns are individuals, and this becomes an adjacency matrix.

# Start interaction loop
for(i in 1:iter){ # for each simulation
  
  # Initialize nodes to zeros state (done once at the beginning of the simulation)
  nodes <- data.frame(nodeID = 1:n, # assign each node an ID
              intProb = pint, # interaction probability
              timeFromLastInt = 0, # number of time steps since node's last encounter with another node
              nInt = 0) # number of interactions this node had
  
  # Main loop:
  for(day in 1:tmax){
    cat(paste("day", day, "\n"))
    
    # Randomize which nodes interact on this day
    # empty edge list to fill with interactions
    interactions <- matrix(nrow = 0, ncol = 2) 
    
    # Interact each node with each other node
    for(nodeA in 1:n){ 
      for(nodeB in (nodeA+1):n){ # don't need to do the dyads twice; hence nodeA+1
        
        # Ignore impossible values of nodeB
        if(nodeB > n) next
        # Calculate the probability of the two nodes meeting
        probMeeting <- nodes[nodeA, "intProb"]*nodes[nodeB, "intProb"]
        if(runif(1) < probMeeting){ # random draw between 0 and 1
          # add the interaction between A and B to the interactions matrix.
          interactions <- rbind(interactions, 
                                c(nodeA, nodeB))
        }
      }
    }
    
    # Get only unique interactions, using igraph (direction doesn't matter)
    uniqInteractions <- interactions %>%
      graph_from_data_frame() %>% 
      simplify() %>%
      as_data_frame()
    
    # Check that we haven't exceeded the number of allowed interactions per day
    # If we have, remove some interactions at random.
    if(nrow(uniqInteractions) > m){
      # how many interactions do we need to remove?
      howManyToRemove <- m - nrow(uniqInteractions)
      # select random rows to remove
      whichToRemove <- sample(1:nrow(uniqInteractions), 
                              size = howManyToRemove)
      # ...and remove the rows.
      uniqInteractions <- uniqInteractions[-whichToRemove,]
    }
    
    # Calculate per-node stats for this day of interactions
    interactingNodes <- uniqInteractions %>%
      pivot_longer(cols = everything()) %>%
      pull() %>% as.integer() # get non-unique vector of all nodes that interacted
    
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
    
    # Create 
    
  } # close day
} # close simulation iteration
