computeProbs <- function(graphList){
  # Get a complete list of all possible edges
  completeGraph <- do.call(igraph::union, graphList)
  complete_edgelist <- expand.grid(names(igraph::V(completeGraph)),
                                   names(igraph::V(completeGraph))) %>%
    dplyr::filter(as.character(.data$Var1) < as.character(.data$Var2)) %>%
    as.matrix()
  
  # get edge list for each element of the graph list
  els <- lapply(graphList, igraph::get.edgelist)
  
  # for each graph, check whether each edge is present or not
  tf <- lapply(els, function(x){
    complete_edgelist %in% x
  })
  
  # bind into a data frame showing presence/absence of edges over time.
  overTime <- do.call(cbind, tf) %>% as.data.frame()
  
  # add two blank steps before (all FALSE), and add the data for the individuals that make up each edge.
  beforeSteps <- data.frame(stepPrevPrev = rep(FALSE, nrow(overTime)),
                            stepPrev = rep(FALSE, nrow(overTime)))
  overTime <- cbind(stats::setNames(as.data.frame(complete_edgelist),
                                    c("ID1", "ID2")), beforeSteps, overTime) # this is our final history data frame
  
  # Create a data frame of probabilities based on overTime
  histdf <- data.frame("add00" = NA, "add10" = NA, "lose01" = NA, "lose11" = NA)
  for(i in 5:ncol(overTime)){
    vec <- vector(mode = "character", nrow(overTime))
    vec[which(!overTime[,i-2] & !overTime[,i-1])] <- "hist00"
    vec[which(!overTime[,i-2] & overTime[,i-1])] <- "hist01"
    vec[which(overTime[,i-2] & !overTime[,i-1])] <- "hist10"
    vec[which(overTime[,i-2] & overTime[,i-1])] <- "hist11"
    
    # compute the probabilities
    add00 <- sum(overTime[i] & vec == "hist00")/sum(vec == "hist00")
    add10 <- sum(overTime[i] & vec == "hist10")/sum(vec == "hist10")
    lose01 <- sum(!overTime[i] & vec == "hist01")/sum(vec == "hist01")
    lose11 <- sum(!overTime[i] & vec == "hist11")/sum(vec == "hist11")
    
    histdf[i-4,] <- c("add00" = add00, "add10" = add10, "lose01" = lose01, "lose11" = lose11)
  }
  
  histdfLong <- histdf %>%
    dplyr::mutate(earlyDate = names(overTime)[-1:-4]) %>%
    tidyr::pivot_longer(cols = -.data$earlyDate, names_to = "type", values_to = "prob")
  
  return(histdfLong)
}
