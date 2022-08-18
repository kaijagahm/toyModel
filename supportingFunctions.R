# Supporting functions for the modelCode_rewiring.R script.
# AUTHOR: Kaija Gahm
# MODIFIED: 2022-06-15
# Just separating these out to keep the other code cleaner and avoid too much scrolling.

library(checkmate)
library(dplyr)

# Define functions --------------------------------------------------------

# Small utility functions for deduplicating edges
# uniqueEdges -------------------------------------------------------------
uniqueEdges <- function(n, triangle = "upper"){
  checkmate::assertNumeric(n, len = 1)
  df <- expand.grid(from = 1:n, to = 1:n)
  if(triangle == "upper"){
    df <- df[which(df[,1] < df[,2]), , drop = FALSE] 
  }else if(triangle == "lower"){
    df <- df[which(df[,1] < df[,2]), , drop = FALSE] 
  }else{
    stop("Triangle must be 'upper' or 'lower'")
  }
  return(df)
}

# dedup -------------------------------------------------------------------
dedup <- function(mat, triangle = "upper"){
  checkmate::assertMatrix(mat)
  if(triangle == "upper"){
    mat <- mat[which(mat[,1] < mat[,2]), , drop = FALSE]
  }else if(triangle == "lower"){
    mat <- mat[which(mat[,1] > mat[,2]), , drop = FALSE]
  }else{
    stop("Triangle must be 'upper' or 'lower'")
  }
  return(mat)
}

# update.network ----------------------------------------------------------
# Function for baseline network dynamics
# Repeat this number of times specified for desired burn.in, in a for loop. Each time, spitting out the network, and the history of the edges.
update.network <- function(ind, # starting index for the history list.
                           network.history, # list of history, since this function depends on being able to look a few timesteps back.
                           add00 = c(0.5, 7), 
                           lose01 = 0.3, 
                           add10 = 0.3, 
                           lose11 = c(0.3, 0.3)){ 
  
  # argument checks
  checkmate::assertNumeric(add00, len = 2)
  checkmate::assertNumeric(lose11, len = 2)
  checkmate::assertNumeric(lose01, len = 1)
  checkmate::assertNumeric(add10, len = 1)
  checkmate::assertList(network.history)
  checkmate::assertNumeric(ind, len = 1)
  
  # Establish a network history, two steps back
  if(ind == 1){
    stop("network must have at least some history--ind cannot be = 1")
  }else if(ind == 2){
    prev <- network.history[[ind-1]]
    # create a matrix of zeroes if the history doesn't exist that far back
    prevprev <- matrix(data = 0, 
                       nrow = nrow(prev), ncol = ncol(prev))
    new <- prev
  }else if(ind > 2){
    prev <- network.history[[ind-1]]
    prevprev <- network.history[[ind-2]]
    new <- prev
  }
  
  # sort edges by history, two back
  h00 <- dedup(which(prev == prevprev & prev == 0, arr.ind = T), "upper")
  h11 <- dedup(which(prev == prevprev & prev == 1, arr.ind = T), "upper")
  h01 <- dedup(which(prevprev < prev, arr.ind = T), "upper")
  h10 <- dedup(which(prevprev > prev, arr.ind = T), "upper")
  N <- nrow(prev)

  # Modify the new adjacency matrix (upper triangle only)
  new[h00] <- rbinom(n = nrow(h00), size = 1, 
                     prob = rbeta(n = nrow(h00), 
                                  shape1 = add00[1], shape2 = add00[2]))
  new[h11] <- rbinom(n = nrow(h11), size = 1,
                     prob = rbeta(n = nrow(h11),
                                  shape1 = lose11[1], shape2 = lose11[2]))
  new[h01] <- rbinom(n = nrow(h01), size = 1,
                     prob = lose01)
  new[h10] <- rbinom(n = nrow(h10), size = 1,
                     prob = add10)
  
  # Symmetrize the matrix
  new <- as.matrix(sna::symmetrize(new, rule = "upper")) # copy the upper triangle over the lower triangle
  
  return(new)
}


# remove.and.rewire -----------------------------------------------------
# Function to remove a node from the network.
remove.and.rewire <- function(network, previous, previousPrevious,
                                n.removed = 1, id = NULL, 
                                pm, # both bereaved 
                                ps, # one bereaved, one not
                                pa, # neither bereaved
                                histMultiplier,
                                coefAdd = 1, # defaults to 1: proportion of friends is the coefficient to increase by.
                                coefLose = 1) {
  # SETUP
  # Calculate population size in the current network
  N <- nrow(network)
  
  # Select nodes to remove
  if(is.null(id)){
    del <- sample(1:N, n.removed, replace = FALSE)
  }else{
    del <- id
  }
  
  # First, capture the edges involving the removed individual(s)
  edges <- network[del,] # row ([del,]), excluding self col ([,-del])
  if(n.removed == 1){ # if we only removed one individual, `edges` is a vector. Have to convert it back to a matrix.
    edges <- matrix(edges, nrow = 1, byrow = TRUE)
  } # if we removed more than one individual, `edges` is already a matrix.
  
  # UPDATE EDGES TO MAINTAIN HISTORY
  h00 <- dedup(which(previous == previousPrevious & previous == 0, arr.ind = T), "upper")
  h11 <- dedup(which(previous == previousPrevious & previous == 1, arr.ind = T), "upper")
  h01 <- dedup(which(previousPrevious < previous, arr.ind = T), "upper")
  h10 <- dedup(which(previousPrevious > previous, arr.ind = T), "upper")
  
  network[h00] <- rbinom(n = nrow(h00), size = 1, 
                     prob = rbeta(n = nrow(h00), 
                                  shape1 = add00[1], shape2 = add00[2]))
  network[h11] <- rbinom(n = nrow(h11), size = 1,
                     prob = rbeta(n = nrow(h11),
                                  shape1 = lose11[1], shape2 = lose11[2]))
  network[h01] <- rbinom(n = nrow(h01), size = 1,
                     prob = lose01)
  network[h10] <- rbinom(n = nrow(h10), size = 1,
                     prob = add10)
  
  # Symmetrize the node-removed network
  network <- as.matrix(sna::symmetrize(network, rule = "upper")) # copy the upper triangle over the lower triangle
  
  # Remove node
  network[del,] <- NA # rows 
  network[,del] <- NA # columns
  N <- N-1 # update the population size.
  
  # Calculate extent of bereavement--how many friends did each individual lose?
  nFriendsLost <- colSums(edges)
  nFriendsHad <- colSums(previous) # how many friends they had
  propFriendsLost <- nFriendsLost/nFriendsHad
  propFriendsLost[is.nan(propFriendsLost)] <- 0
  propFriendsLostDF <- data.frame(node = 1:length(propFriendsLost),
                                  propFriendsLost = propFriendsLost)
  

  # SAVE AFTER-LOSS NETWORK
  afterLoss <- network # save this as a time step--the network after a loss but before rewiring.
  
  # REWIRING
  prev <- afterLoss
  prevPrev <- previous
  new <- prev
  
  # Calculate baseline probabilities
  h00 <- as.data.frame(dedup(which(prev == prevPrev & prev == 0, arr.ind = T), "upper")) %>%
    setNames(., c("ind1", "ind2"))
  h11 <- as.data.frame(dedup(which(prev == prevPrev & prev == 1, arr.ind = T), "upper")) %>%
    setNames(., c("ind1", "ind2"))
  h01 <- as.data.frame(dedup(which(prevPrev < prev, arr.ind = T), "upper")) %>%
    setNames(., c("ind1", "ind2"))
  h10 <- as.data.frame(dedup(which(prevPrev > prev, arr.ind = T), "upper")) %>%
    setNames(., c("ind1", "ind2"))
  
  h00$baselineProb <- rbeta(n = nrow(h00), shape1 = add00[1], shape2 = add00[2])
  h11$baselineProb <- rbeta(n = nrow(h11), shape1 = lose11[1], shape2 = lose11[2])
  h01$baselineProb <- lose01
  h10$baselineProb <- add10
  
  # Modify baseline probabilities based on friends lost
  h00 <- left_join(h00, propFriendsLostDF, by = c("ind1" = "node")) %>%
    rename("mod1" = propFriendsLost) %>%
    left_join(., propFriendsLostDF, by = c("ind2" = "node")) %>%
    rename("mod2" = propFriendsLost) %>%
    mutate(newProb = baselineProb + ((mod1+mod2)*coefAdd*baselineProb)) # XXX coefAdd
  
  h10 <- left_join(h10, propFriendsLostDF, by = c("ind1" = "node")) %>%
    rename("mod1" = propFriendsLost) %>%
    left_join(., propFriendsLostDF, by = c("ind2" = "node")) %>%
    rename("mod2" = propFriendsLost) %>%
    mutate(newProb = baselineProb + ((mod1+mod2)*coefAdd*baselineProb)) # XXX coefAdd
  
  h11 <- left_join(h11, propFriendsLostDF, by = c("ind1" = "node")) %>%
    rename("mod1" = propFriendsLost) %>%
    left_join(., propFriendsLostDF, by = c("ind2" = "node")) %>%
    rename("mod2" = propFriendsLost) %>%
    mutate(newProb = baselineProb + ((mod1+mod2)*coefLose*baselineProb)) # XXX coefLose
  
  h01 <- left_join(h01, propFriendsLostDF, by = c("ind1" = "node")) %>%
    rename("mod1" = propFriendsLost) %>%
    left_join(., propFriendsLostDF, by = c("ind2" = "node")) %>%
    rename("mod2" = propFriendsLost) %>%
    mutate(newProb = baselineProb + ((mod1+mod2)*coefLose*baselineProb)) # XXX coefLose
  
  # Do the rewiring--calculate new edges
  new[as.matrix(h00[,1:2])] <- rbinom(n = nrow(h00), size = 1, 
                         prob = h00$newProb)
  new[as.matrix(h11[,1:2])] <- rbinom(n = nrow(h11), size = 1,
                         prob = h11$newProb)
  new[as.matrix(h01[,1:2])] <- rbinom(n = nrow(h01), size = 1,
                         prob = h01$newProb)
  new[as.matrix(h10[,1:2])] <- rbinom(n = nrow(h10), size = 1,
                         prob = h10$newProb)
  
  # Symmetrize rewired network
  new <- as.matrix(sna::symmetrize(new, rule = "upper")) # copy the upper triangle over the lower triangle
  
  # SAVE REWIRED NETWORK
  rewired <- new # save this as a time step--the network after rewiring
  
  return(list(networks = list("afterLoss" = afterLoss,
                              "rewired" = rewired),
              del = del))
}