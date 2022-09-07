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
                     edge.prob = 0.04, # Initial edge probability. This is derived from the average network density, when taking a 5-day increment, from parameterizingTheModel.Rmd.
                     burn.in = 20, # How many iterations of baseline dynamics to run before node removals.
                     recovery = 10, # How many iterations of baseline dynamics to run after node removals.
                     add00 = c(0.4721719, 7.3144796), # beta distribution parameters derived from parameterizingTheModel.Rmd.
                     lose01 = 0.3, # Probability of losing an edge with history h01. Derived from real data.
                     add10 = 0.2, # Probability of gaining an edge with history h10. Derived from real data.
                     lose11 = c(0.3283134, 0.3062181),
                     id = NULL,
                     coefAdd = 1,
                     coefLose = -1){ # beta distribution parameters derived from parameterizingTheModel.Rmd.)
  
  # ARGUMENT CHECKS ---------------------------------------------------------
  checkmate::assertInteger(as.integer(N), lower = 2, any.missing = FALSE, len = 1)
  checkmate::assertInteger(as.integer(n.removed), lower = 0, upper = N-1, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(edge.prob, len = 1, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assertInteger(as.integer(burn.in), lower = 2, any.missing = FALSE, len = 1)
  checkmate::assertInteger(as.integer(recovery), lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(add00, len = 2, any.missing = FALSE)
  checkmate::assertNumeric(lose01, len = 1, any.missing = FALSE, upper = 1, lower = 0)
  checkmate::assertNumeric(add10, len = 1, any.missing = FALSE, upper = 1, lower = 0)
  checkmate::assertNumeric(lose11, len = 2, any.missing = FALSE)
  if(!is.null(id)){
    checkmate::assertInteger(as.integer(id), len = n.removed, lower = 1, upper = N, unique = T)
  }
  checkmate::assertNumeric(coefAdd, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(coefLose, len = 1, any.missing = FALSE)
  
  # BURN-IN -----------------------------------------------------------------
  ## Setup: assign individual sociabilities
  ## By trial and error (using this app https://homepage.divms.uiowa.edu/~mbognar/applets/beta.html), I found alpha (14) and beta (8) parameters for a beta distribution that yields a mean of approximately 0.63 and seems to have a reasonable spread (this is based on nothing besides eyeballing it, to be totally honest).
  ## Why do I want the mean to be 0.63? Because initially I had the network density at 0.4, and the square root of 0.4 is around 0.63. If I'm right about this, drawing sociabilities from a distribution centered around 0.63 and then randomly interacting those individuals should yield an overall network density of around 0.4. Let's test it out.
  soc <- data.frame(indiv = 1:N,
                    sociability = rbeta(N, shape1 = 14, shape2 = 8))
  probMatrix <- soc$sociability %*% t(soc$sociability)
  
  ## Empty list to hold networks
  network.history <- vector(mode = "list", length = burn.in + 1)
  
  ## Create random networks for the first two time steps 
  network.history[[1]] <- sna::symmetrize(matrix(rbinom(N*N, 1, probMatrix), N, N), 
                                          rule = "upper")
  network.history[[2]] <- sna::symmetrize(matrix(rbinom(N*N, 1, probMatrix), N, N), 
                                          rule = "upper")
  ## Starting at element 3, 
  
  ## {to remove} Create the rest of the burn-in-period networks, starting at the 3rd element.
  ## The last element (network.history[[burn.in+1]]) is the one that will be modified by removing node(s).
  for(i in 3:(burn.in+1)){
    output <- update.network(ind = i, network.history, add00 = add00, 
                             add10 = add10, lose01 = lose01, lose11 = lose11)
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
  
  # SEE FRIENDS OF REMOVED NODES ----------------------------------------
  ## Note: decided that for now, we're considering "friends" to be "the ones what were connected to the now-removed node(s) in timestep `back1`", as opposed to "the ones that would have been connected to the now-removed node(s) in timestep `removed`. I'm not positive that this is correct, but I'm going with it for now.
  friendships <- network.history[["back1"]][del,]
  if(n.removed == 1){ # if we only removed one individual, `friendships` is a vector. Have to convert it back to a matrix.
    friendships <- matrix(friendships, nrow = 1, byrow = TRUE)
  } # if we removed more than one individual, `friendships` is already a matrix.
  
  # REWIRING ----------------------------------------------------------------
  ## CALCULATE BEREAVEMENT ---------------------------------------------------
  ## What proportion of its friends did each individual lose?
  nFriendsLost <- colSums(friendships)
  nFriendsHad <- colSums(network.history[["back1"]]) # XXX but see, this still isn't making sense to me. Because then where do the connections from network.history[["removed"]] fit?? Come back to this.
  propFriendsLost <- nFriendsLost/nFriendsHad
  propFriendsLost[is.nan(propFriendsLost)] <- 0
  propFriendsLost[del] <- NA # set NA's for the removed nodes
  propFriendsLostDF <- data.frame(node = 1:length(propFriendsLost),
                                  propFriendsLost = propFriendsLost)
  
  ## CALCULATE BASELINE PROBABILITIES ----------------------------------------
  back1 <- network.history[["back1"]]
  removed <- network.history[["removed"]]
  
  allEdges <- as.matrix(uniqueEdges(N))
  rem <- removed[allEdges]
  b1 <- back1[allEdges]
  allEdges <- as.data.frame(allEdges) %>%
    mutate(back1 = b1,
           removed = rem) %>%
    mutate(history = paste0("h", back1, removed),
           baselineProb = case_when(history == "h00" ~ rbeta(1, shape1 = add00[1], shape2 = add00[2]),
                                    history == "h10" ~ add10,
                                    history == "h11" ~ rbeta(1, shape1 = lose11[1], shape2 = lose11[2]),
                                    history == "h01" ~ lose01))
  
  ## MODIFY PROBS WITH INFORMATION ABOUT FRIENDS LOST ------------------------
  # Join information about how many friends each individual lost (need to do two joins, since there are two individuals involved in each edge.)
  allEdges <- allEdges %>%
    left_join(., propFriendsLostDF, by = c("from" = "node")) %>%
    rename("mod1" = propFriendsLost) %>%
    left_join(., propFriendsLostDF, by = c("to" = "node")) %>%
    rename("mod2" = propFriendsLost) %>%
    # XXX TALK THIS THROUGH WITH SOMEONE TO MAKE SURE IT'S CLEAR.
    mutate(newProb = case_when(history %in% c("h00", "h10") ~ 
                                 baselineProb + ((mod1+mod2)*coefAdd*baselineProb),
                               history %in% c("h11", "h01") ~ 
                                 1-(baselineProb + ((mod1+mod2)*coefLose*baselineProb)),
                               TRUE ~ baselineProb)) # XXX change this to NA so it's clearer if something actually has gone wrong.
  # XXX add a check to make sure there are no NA's--there should be none.
  
  ## DRAW NEW EDGES ----------------------------------------------------------
  newEdges <- suppressWarnings(rbinom(1:nrow(allEdges), 1, prob = allEdges$newProb))
  #  note that all these probabilities have been converted to probabilities of "success", aka probability of the edge *existing*.
  allEdges$rewired <- newEdges
  
  ## CREATE ADJ MATRIX FOR REWIRED NETWORK -----------------------------------
  rewired <- removed # initialize network of same size
  rewired[as.matrix(allEdges[,c("from", "to")])] <- allEdges$rewired
  
  ## SYMMETRIZE REWIRED NETWORK ----------------------------------------------
  rewired <- as.matrix(sna::symmetrize(rewired, rule = "upper"))
  
  ## MAKE SURE ALL EDGES THAT SHOULD BE NA ARE NA ----------------------------
  rewired[del,] <- NA # set rows to NA
  rewired[,del] <- NA # set cols to NA
  
  # ADD `REWIRED` TO NETWORK.HISTORY LIST --------------------------------------
  network.history <- append(network.history, list("rewired" = rewired))
  
  # CONTINUE BASELINE DYNAMICS UNTIL THE END -----------------------------------
  network.history <- append(network.history, rep(NA, recovery))
  names(network.history)[(which(names(network.history) == "rewired")+1):length(network.history)] <- paste0("recovery_", 1:recovery)
  
  for(i in (burn.in+3):length(network.history)){
    output <- update.network(ind = i, network.history, add00 = add00, 
                             add10 = add10, lose01 = lose01, lose11 = lose11)
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
    igraph::graph_from_adjacency_matrix(x, mode = "undirected", add.colnames = "label")
  })
  
  # RETURN OUTPUTS AS LIST --------------------------------------------------
  # Network history (with NA's); network history (NAs removed); graphs; indices of deleted nodes
  return(list("network.history.nas" = network.history,
              "network.history.nodesRemoved" = network.history.nodesRemoved,
              "graphs" = graphs,
              "whichRemoved" = del))
}