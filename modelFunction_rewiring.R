# Function to run the model
# Architecture of this version is based more directly on Farine's 2nd-degree rewiring model (2021)
# AUTHOR: Kaija Gahm (kgahm@ucla.edu)
# REVISION DATE: 2022-07-06 (added parameters derived in parameterizingTheModel.Rmd)

library(here)
source(here("supportingFunctions.R")) # all the functions that will be used in the main model

runModel <- function(N = 50, # Nodes in the network
                     nodes.removed = 1, # Nodes to remove
                     n.removed = 1, # How many to remove at a time
                     edge.prob = 0.04, # This is derived from the average network density, when taking a 5-day increment, from parameterizingTheModel.Rmd.
                     burn.in = 50,
                     burn.out = 50,
                     pm = 0.3,
                     ps = 0.1, 
                     pa = 0.2,
                     add00 = c(0.4721719, 7.3144796), # beta distribution parameters derived from parameterizingTheModel.Rmd.
                     lose01 = 0.3, 
                     add10 = 0.2,
                     lose11 = c(0.3283134, 0.3062181), # beta distribution parameters derived from parameterizingTheModel.Rmd.
                     histMultiplier = 1.2,
                     doRemoval = TRUE){
  # ARGUMENT CHECKS
  # XXX update these.
  
  # RANDOM STARTING NETWORK
  network.orig <- sna::rgraph(N, tprob = edge.prob, 
                              mode = "graph") # gives undirected graph, already symmetrized.
  
  # BURN-IN
  ## Empty list to hold networks
  network.history <- vector(mode = "list", length = burn.in + 1)
  ## Create pre-history
  network.history[[1]] <- matrix(0, N, N) # blank network to enable looking 2 timesteps back
  network.history[[2]] <- network.orig # baseline network
  
  ## Create the rest of the burn-in-period networks, starting at the 3rd element.
  ## The last element (network.history[[burn.in+1]]) is the one that will be modified by removing node(s).
  for(i in 3:burn.in+1){
    output <- update.network(ind = i, network.history, add00 = add00, 
                             add10 = add10, lose01 = lose01, lose11 = lose11)
    network.history[[i]] <- output # update history
  }
  
  # Removal/perturbation and rewiring
  rewired.list <- remove.and.rewire(network = network.history[[burn.in]], 
                                    previous = network.history[[burn.in-1]],
                                    histMultiplier = histMultiplier,
                                    n.removed = 1, pm = rnorm(1, pm), 
                                    ps = rnorm(1, ps), pa = rnorm(1, pa))
  rewired.network <- rewired.list$network # the actual rewired network
  rewired.del <- rewired.list$del # removed individuals. Used for amending history matrices in order to continue baseline dynamics moving forward
  
  # Continue baseline dynamics following removal
  ## Set up a list: the old network history, the rewired network, and some placeholders for the burnout
  history <- append(network.history, list(rewired.network))
  history <- append(history, rep(NA, burn.out)) # add spots for the burn out
  
  for(i in (burn.in+2):length(history)){
    output <- update.network(ind = i, history, add00 = add00, 
                             add10 = add10, lose01 = lose01, lose11 = lose11)
    history[[i]] <- output # update history
  }
  
  # Make sure column names are assigned so we can keep track of the nodes!
  history <- lapply(history, function(x){
    colnames(x) <- paste0("v", 1:ncol(x))
    return(x)
  })
  
  # Now we have to remove the deleted node(s) from all networks following the removal.
  for(i in (burn.in+1):length(history)){
    history[[i]] <- history[[i]][-rewired.del,]
    history[[i]] <- history[[i]][,-rewired.del]
  }
  return(history)
}
