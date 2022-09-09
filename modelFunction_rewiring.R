# Function to run the model
# AUTHOR: Kaija Gahm (kgahm@ucla.edu)
# REVISION DATE: 2022-08-18 (completely revamped! Added the ability to remove multiple individuals, and changed how rewiring is done. Also finally separated removal and rewiring into two separate time slices.)
# Previous iterations of this model were based more directly on Farine (2021)'s second-degree rewiring idea. The current model expands beyond removing a single individual. Now that we remove two individuals, we keep the *idea* of second-degree rewiring by adjusting the probability of each edge according to how many friends each of the involved nodes lost.

# I think it will end up that second-degree edges are more likely to form if you keep the coefficients as 1 and -1. But need to actually test that hypothesis.

# PACKAGES ----------------------------------------------------------------
library(here)
source(here("supportingFunctions.R")) # all the functions that will be used in the main model

# MODEL FUNCTION DEF ------------------------------------------------------
runModel <- function(N = 50, # Number of nodes in the starting network. Must be an integer > 2. Default is 50.
                     n.removed = 1, # Number of nodes to remove. Must be an integer >= 0 and <= N-1. Default is 1.
                     socAlpha = 14,
                     socBeta = 8,
                     burn.in = 20, # How many iterations of baseline dynamics to run before node removals.
                     recovery = 10, # How many iterations of baseline dynamics to run after node removals.
                     mod00 = -0.2, 
                     mod01 = 0.1, 
                     mod10 = -0.1, 
                     mod11 = 0.2,
                     id = NULL,
                     coefBereavement = 1){ # beta distribution parameters derived from parameterizingTheModel.Rmd.)
  
  # ARGUMENT CHECKS ---------------------------------------------------------
  checkmate::assertInteger(as.integer(N), lower = 2, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(socAlpha, len = 1, lower = 0)
  checkmate::assertNumeric(socBeta, len = 1, lower = 0)
  checkmate::assertInteger(as.integer(n.removed), lower = 0, upper = N-1, any.missing = FALSE, len = 1)
  checkmate::assertInteger(as.integer(burn.in), lower = 2, any.missing = FALSE, len = 1)
  checkmate::assertInteger(as.integer(recovery), lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(mod00, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(mod01, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(mod10, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(mod11, len = 1, any.missing = FALSE)
  if(!is.null(id)){
    checkmate::assertInteger(as.integer(id), len = n.removed, lower = 1, upper = N, unique = T)
  }
  checkmate::assertNumeric(coefBereavement, len = 1, any.missing = FALSE)

  # BURN-IN -----------------------------------------------------------------
  ## Setup: assign individual sociabilities
  ## By trial and error (using this app https://homepage.divms.uiowa.edu/~mbognar/applets/beta.html), I found alpha (14) and beta (8) parameters for a beta distribution that yields a mean of approximately 0.63 and seems to have a reasonable spread (this is based on nothing besides eyeballing it, to be totally honest).
  ## Why do I want the mean to be 0.63? Because initially I had the network density at 0.4, and the square root of 0.4 is around 0.63. If I'm right about this, drawing sociabilities from a distribution centered around 0.63 and then randomly interacting those individuals should yield an overall network density of around 0.4. Let's test it out.
  soc <- data.frame(indiv = 1:N,
                    sociability = rbeta(N, shape1 = socAlpha, shape2 = socBeta))
  probMatrix <- soc$sociability %*% t(soc$sociability)
  
  ## Empty list to hold networks
  network.history <- vector(mode = "list", length = burn.in + 1)
  
  ## Create random networks for the first two time steps 
  network.history[[1]] <- sna::symmetrize(matrix(rbinom(N*N, 1, probMatrix), N, N), 
                                          rule = "upper")
  network.history[[2]] <- sna::symmetrize(matrix(rbinom(N*N, 1, probMatrix), N, N), 
                                          rule = "upper")
  
  ## Starting at element 3, update the probMatrix with the `mod` values. Then create new adjacency matrices.
  ## The last element (network.history[[burn.in+1]]) is the one that will be modified by removing node(s).
  for(i in 3:(burn.in+1)){
    output <- update.network(ind = i, network.history, 
                             mod00 = mod00, 
                             mod11 = mod11, 
                             mod10 = mod10, 
                             mod01 = mod01,
                             probMatrix = probMatrix)
    network.history[[i]] <- output # update history
  }
  
  # assign names so we can keep track of where the removal happened
  names(network.history) <- c(paste0("history_", 1:(burn.in-2)), "back2", "back1", "removed")
  # keep in mind that even though we've called the (burn.in+1)th slice "removed", it does not actually have the node(s) removed (aka set to NA) yet. We'll do that next.
  
  # SELECT WHICH NODES TO REMOVE --------------------------------------------
  if(is.null(id)){
    del <- sample(1:N, n.removed, replace = FALSE)
  }else{
    del <- id
  }
  
  # REMOVE NODES ------------------------------------------------------------
  network.history[["removed"]][del,] <- NA # set rows to NA
  network.history[["removed"]][,del] <- NA # set cols to NA
  
  # SEE CONNECTIONS OF REMOVED NODES ----------------------------------------
  ## Note: decided that for now, we're considering "connections" to be "the individuals what were connected to the now-removed node(s) in timestep `back1`", as opposed to "the individuals that would have been connected to the now-removed node(s) in timestep `removed`. I'm not positive that this is correct, but I'm going with it for now.
  conns <- network.history[["back1"]][del,]
  if(n.removed == 1){ # if we only removed one individual, `conns` is a vector. Have to convert it back to a matrix.
    conns <- matrix(conns, nrow = 1, byrow = TRUE)
  } # if we removed more than one individual, `conns` is already a matrix.
  
  # REWIRING ----------------------------------------------------------------
  ## CALCULATE BEREAVEMENT ---------------------------------------------------
  ## What proportion of its connections did each individual lose?
  nConnsLost <- colSums(conns)
  nConnsHad <- colSums(network.history[["back1"]]) 
  propConnsLost <- nConnsLost/nConnsHad
  propConnsLost[is.nan(propConnsLost)] <- 0
  propConnsLost[del] <- NA # set NA's for the removed nodes
  propConnsLostDF <- data.frame(node = 1:length(propConnsLost),
                                  propConnsLost = propConnsLost)
  
  ## Compute alterations to baseline probabilities based on loss, incorporating `coefBereavement`.
  withLoss <- soc$sociability*(1+propConnsLost)*coefBereavement
  
  # Compute new probability matrix
  probMatrixWithLoss <- withLoss %*% t(withLoss)
  
  ## MODIFY PROBS WITH INFORMATION ABOUT EDGE HISTORY
  # Define histories XXX this is redundant with the internal identification of histories in the `update.network` function.
  prev <- network.history[["removed"]]
  prevprev <- network.history[["back1"]]
  h00 <- dedup(which(prev == prevprev & prev == 0, arr.ind = T), "upper")
  h11 <- dedup(which(prev == prevprev & prev == 1, arr.ind = T), "upper")
  h01 <- dedup(which(prevprev < prev, arr.ind = T), "upper")
  h10 <- dedup(which(prevprev > prev, arr.ind = T), "upper")
  
  newProbMatrix <- probMatrixWithLoss
  
  newProbMatrix[h00] <- newProbMatrix[h00]*(1+mod00)
  newProbMatrix[h11] <- newProbMatrix[h11]*(1+mod11)
  newProbMatrix[h01] <- newProbMatrix[h01]*(1+mod01)
  newProbMatrix[h10] <- newProbMatrix[h10]*(1+mod10)
  
  ## CREATE ADJ MATRIX FOR REWIRED NETWORK -----------------------------------
  rewiredAdj <- suppressWarnings(matrix(rbinom(N*N, 1, newProbMatrix), N, N))
  
  ## SYMMETRIZE REWIRED NETWORK ----------------------------------------------
  rewiredAdj <- as.matrix(sna::symmetrize(rewiredAdj, rule = "upper"))

  # ADD `REWIRED` TO NETWORK.HISTORY LIST --------------------------------------
  network.history <- append(network.history, list("rewired" = rewiredAdj))
  
  # CONTINUE BASELINE DYNAMICS UNTIL THE END -----------------------------------
  network.history <- append(network.history, rep(NA, recovery))
  names(network.history)[(which(names(network.history) == "rewired")+1):length(network.history)] <- paste0("recovery_", 1:recovery)
  
  for(i in (burn.in+1):length(network.history)){
    output <- update.network(ind = i, network.history, 
                             mod00 = mod00, 
                             mod11 = mod11, 
                             mod10 = mod10, 
                             mod01 = mod01,
                             probMatrix = probMatrix)
    network.history[[i]] <- output # update history
  }
  
  # NAME VERTICES --------------------------------------------------------------
  network.history <- lapply(network.history, function(x){
    colnames(x) <- paste0("v", 1:ncol(x))
    return(x)
  })
  
  # REMOVE DELETED NODES -------------------------------------------------------
  ## (instead of just setting them to NA)
  ## If we don't do this, then they will just be treated as isolated nodes, which isn't what we need to do.
  network.history.nodesRemoved <- network.history
  for(i in (burn.in+1):length(network.history.nodesRemoved)){
    network.history.nodesRemoved[[i]] <- network.history.nodesRemoved[[i]][-del,] # remove rows
    network.history.nodesRemoved[[i]] <- network.history.nodesRemoved[[i]][,-del] # remove cols
  }
  
  # MAKE GRAPHS -------------------------------------------------------------
  ## Make the graphs (now with the right number of nodes)
  graphs <- lapply(network.history.nodesRemoved, function(x){
    igraph::graph_from_adjacency_matrix(x, mode = "undirected", diag = FALSE, add.colnames = "label")
  })
  
  # RETURN OUTPUTS AS LIST --------------------------------------------------
  # Network history (with NA's); network history (NAs removed); graphs; indices of deleted nodes
  return(list("network.history.nas" = network.history,
              "network.history.nodesRemoved" = network.history.nodesRemoved,
              "graphs" = graphs,
              "whichRemoved" = del))
}