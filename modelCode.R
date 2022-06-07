library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(data.table)

### Define utility functions for later use
# FOR USE IN THE MODEL:
# Function that searches for any nodes that are not connected at 
# all to other nodes, and attaches them to one other node at random.
connectIsolatedNodes <- function(g, allow){
  if(any(degree(g) == 0) & allow == FALSE){
    isolatedNodes <- which(degree(g) == 0)
    gOut <- g
    for(focalNode in isolatedNodes){
      gOut <- add_edges(gOut, 
                        edges = c(focalNode, #from the focal node...
                                  #...to any other node
                                  # (but not self!)
                                  sample((1:length(V(gOut)))[-focalNode], 
                                         size = 1)))
    }
    return(gOut)
  }else{
    return(g) # pass the graph through as it is, if there are no isolated nodes
  }
}

# Get only unique edges (since we're dealing with an undirected graph)
uniqueEdges <- function(n){
  df <- expand.grid(from = 1:n, to = 1:n) %>%
    mutate(inOrder = case_when(from == to ~ NA_character_, 
                               from > to ~ paste(to, from),
                               from < to ~ paste(from, to))) %>%
    filter(!is.na(inOrder)) %>% # remove self edges
    group_by(inOrder) %>%
    slice(1) %>% # take only one edge for each
    ungroup() %>%
    select(-inOrder)
  return(df)
}

# FOR USE IN ANALYZING THE RESULTS:
# Make a bunch of plots to show change in the network over the model run
plotSim <- function(modelOutput, pointsize = 5, edgewidth = 0.2){
  # Create coords for plotting based on the first day's network
  coords <- layout_with_fr(modelOutput[[1]])
  
  return(lapply(modelOutput, function(x){
    x %>%
      ggraph(layout = coords)+
      geom_edge_link(edge_width = edgewidth)+
      geom_node_point(size = pointsize)
  }))
}

# Get node-level stats for the model simulation output
getNodeStats <- function(modelOutput, type = "df"){
  # Check to make sure the "type" argument is valid
  if(!(type %in% c("df", "graphs", "list"))){
    stop("Argument 'type' must be 'df', 'graphs', or 'list'.")
  }
  
  # Calculate stats
  stats <- lapply(modelOutput, function(x){
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

### The main modeling function
# MODEL:
# Returns a list of igraph network objects (graphs), one for each day in 1:tmax. 
# Graphs are BINARY and UNDIRECTED, based on the UPPER TRIANGLE of the adjacency matrix. 
runSim <- function(tmax = 10, # length of time over which to run the simulation
                   n = 50, # number of individuals
                   # allow nodes to be unconnected? If F, randomly connects 
                   # each unconnected node with one other node.
                   allowIsolated = FALSE, 
                   # p gain an edge given not connected in either of the 
                   # previous two time steps
                   add00 = 0.1, 
                   # p lose an edge given connected in prev time step but 
                   # not prev prev
                   lose01 = 0.3, 
                   # p gain an edge given connected in prev prev time step 
                   # but not prev
                   add10 = 0.3, 
                   # p lose edge given connected in previous 2 time steps
                   lose11 = 0.1,
                   verbose = TRUE){ 
  # STORAGE
  # Initialize lists to store graphs and ajacency matrices for each time step
  # storage for graphs for each day
  gs <- vector(mode = "list", length = tmax) 
  # store edge info for each day
  ei <- vector(mode = "list", length = tmax) 
  
  # Make a list of all the unique node pairs
  edges <- uniqueEdges(n)
  
  # SETUP
  # Create the graph for the first day
  startingEdges <- edges %>%
    # sample some edges based on the prob of creating edges out of nowhere
    sample_n(size = nrow(edges)*add00) 
  day1G <- graph_from_data_frame(startingEdges, directed = FALSE, vertices = 1:n)
  
  # Deal with isolated nodes
  # The connectIsolatedNodes function has a built-in switch based on the 
  # allowIsolated parameter.
  day1G <- connectIsolatedNodes(g = day1G, allow = allowIsolated)
  gs[[1]] <- day1G
  
  # RUN SIMULATION
  # Loop through the days
  progress <- 0 # initialize progress indicator
  #cat("0%\n")
  for(day in 2:tmax){
    if(verbose){
      # minimal feedback while running the model, for sanity check on speed
      cat(paste("day", day, "\n")) 
    }
    
    prevG <- gs[[day-1]] # previous day's matrix to operate on
    if(day < 3){
      # if we don't have enough info to get the am from two days ago, 
      # make a zeroes graph
      prevprevG <- graph_from_adjacency_matrix(matrix(0, n, n),
                                               mode = "upper", diag = FALSE) 
    }else{
      prevprevG <- gs[[day-2]]
    }
    
    # Use get.edge.ids to search for the edge, using directed = FALSE
    # Annoyingly, have to use a for loop with get.edge.ids--returns mysterious 
    # weird results if I try to use it with mutate. Not sure why. I wonder if 
    # there's a vectorized alternative in tidygraph. XXX Need to investigate 
    # this further.
    edgesInfo <- edges
    edgesInfo$idInPrevPrev <- NA
    edgesInfo$idInPrev <- NA
    for(edge in 1:nrow(edges)){
      from <- edges[edge, "from"]
      to <- edges[edge, "to"]
      edgesInfo$idInPrevPrev[edge] <- get.edge.ids(prevprevG, 
                                                   vp = c(from, to), directed = F)
      edgesInfo$idInPrev[edge] <- get.edge.ids(prevG, 
                                               vp = c(from, to), directed = F)
    }
    # draw random number for each individual edge
    edgesInfo <- edgesInfo %>% mutate(rand = runif(n = nrow(edges), 
                                                   min = 0, max = 1)) %>% 
      # Characterize the type of relationship between each pair of nodes
      mutate(case = case_when(idInPrev == 0 & idInPrevPrev == 0 ~ "h00",
                              idInPrev != 0 & idInPrevPrev == 0 ~ "h01",
                              idInPrev == 0 & idInPrevPrev != 0 ~ "h10",
                              idInPrev != 0 & idInPrevPrev != 0 ~ "h11",
                              TRUE ~ "error"))
    
    # Error message just in case!
    if(any(edgesInfo$case == "error")){
      stop("error computing edge cases!")
    }
    
    # Rules (note that "stay" is comprehensible but not actionable, 
    # so I've translated it into ADD/LOSE)
    edgesInfo <- edgesInfo %>%
      mutate(action = case_when(case == "h00" & rand <= add00 ~ "add",
                                case == "h00" & rand > add00 ~ "none",
                                case == "h01" & rand <= lose01 ~ "lose",
                                case == "h01" & rand > lose01 ~ "none",
                                case == "h10" & rand <= add10 ~ "add",
                                case == "h10" & rand > add10 ~ "none",
                                case == "h11" & rand <= lose11 ~ "lose",
                                case == "h11" & rand > lose11 ~ "none",
                                TRUE ~ "error"))
    
    # Error message just in case!
    if(any(edgesInfo$action == "error")){
      stop("error computing edge actions!")
    }
    
    # Loop through the edges and operate on them
    # Idea for later: df[,c("from", "to")] %>% as.matrix() %>% t() %>% 
    # as.vector() to be fed into add_edges so I can do them all at once 
    # (assuming that df is edgesInfo %>% filter(action == "add)). But no 
    # point in refactoring the code to do this right now since it turns 
    # out it doesn't work for delete_edges.
    newG <- prevG # initialize h01 graph
    for(edge in 1:nrow(edgesInfo)){
      pair <- c(edgesInfo[edge, "from"], edgesInfo[edge, "to"])
      if(edgesInfo[edge, "action"] == "add"){
        newG <- add_edges(graph = newG, edges = pair)
      }else if(edgesInfo[edge, "action"] == "lose"){
        edgeID <- get.edge.ids(newG, vp = pair, directed = F) # XXXXXXXXXXXX
        newG <- delete_edges(graph = newG, edgeID) 
        # XXX UNBELIEVABLE--WHY DOES DELETE_EDGES WORK DIFFERENTLY THAN ADD_EDGES??????
        # In the cold light of day I see that this is probably because igraph 
        # allows for networks with more than one edge between nodes. Okay, 
        # I concede that this is logical. But they should be clearer about it 
        # in the documentation, and they should allow two options for 
        # specifying this.
      }else{
        newG <- newG
      }
      # print(paste0("iteration ", edge, " complete")) # for debug
    }
    
    # If we're not allowing isolated nodes, randomly connect each node one time.
    newG <- connectIsolatedNodes(g = newG, allow = allowIsolated)
    
    # Save this day's graph
    gs[[day]] <- newG
    # Save this day's edge info
    ei[[day]] <- edgesInfo
    if((day/tmax) >= progress + 0.1){
      progress <- day/tmax
      #cat(paste0(progress*100, "%\n"))
    }
  } # close day
  
  return(list("gs" = gs, "ei" = ei))
}