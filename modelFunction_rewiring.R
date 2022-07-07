# Function to run the model
# Architecture of this version is based more directly on Farine's 2nd-degree rewiring model (2021)
# AUTHOR: Kaija Gahm (kgahm@ucla.edu)
# REVISION DATE: 2022-07-06 (added parameters derived in parameterizingTheModel.Rmd)

source("supportingFunctions.R") # all the functions that will be used in the main model

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
    # Generate a random starting network
    network.orig <- sna::rgraph(N, tprob = edge.prob, 
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
    
    if(doRemoval == FALSE){
      return(network.history)
    }
    
    if(doRemoval == TRUE){
      # Removal/perturbation and rewiring
      rewired.list <- remove.network.node(network = network.history[[burn.in]], 
                                          previous = network.history[[burn.in-1]],
                                          histMultiplier = histMultiplier,
                                          n.removed = 1, pm = rnorm(1, pm), 
                                          ps = rnorm(1, ps), pa = rnorm(1, pa))
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
      return(revised.history)
    }
}
