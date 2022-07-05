# Sandbox script to play around with the feeding network data

# load packages
library(vultureUtils)
library(tidyverse)

# load data
load("data/feedingEdges2021.Rda")
load("data/feedingPoints2021.Rda")

# Filter to include only southern individuals
indivs <- feedingPoints %>%
  sf::st_drop_geometry() %>%
  dplyr::select(trackId, location_long.1, location_lat.1) %>%
  dplyr::group_by(trackId) %>%
  dplyr::summarize(mnlat = mean(location_lat.1))

southern <- indivs %>%
  dplyr::filter(mnlat < 32) %>%
  dplyr::pull(trackId)

southernEdges2021 <- feedingEdges2021 %>%
  dplyr::filter(ID1 %in% southern & ID2 %in% southern)
  
unweighted <- vultureUtils::makeGraphs(edges = southernEdges2021, interval = "20 days", dateTimeStart = "2021-01-01 00:00:00", dateTimeEnd = "2021-12-31 11:59:59", weighted = FALSE)$graphs

# Plot the unweighted networks
plots <- vultureUtils::plotGraphs(unweighted)

# Make an animated gif
vultureUtils::makeGIF(plots, fileName = "2021_20Days.gif", interval = 0.2)

# Examine network metrics vs. time window ---------------------------------
# Let's look at e.g. network density
timeWindows <- seq(1, 50, by = 10)
intervals <- paste(timeWindows, "days")
dts <- "2021-01-01 00:00:00"
dte <- "2021-12-31 11:59:59"
e <- feedingEdges2021
f <- feedingPoints2021

dat <- vector(mode = "list", length = length(intervals))
names(dat) <- intervals
for(i in 1:length(dat)){
  cat(paste("Computing graphs for an interval of", intervals[i], "\n"))
  dat[[i]] <- makeGraphs(edges = e, interval = intervals[i], dateTimeStart = dts,
                    dateTimeEnd = dte, weighted = FALSE)
}

graphLists <- lapply(dat, function(x){x[["graphs"]]})
names(graphLists) <- intervals

densities <- lapply(graphLists, function(x){
  unlist(lapply(x, igraph::edge_density)) %>%
    as.data.frame() %>%
    setNames("networkDensity")
})
names(densities) <- intervals

# Still need to do some work to be able to make graphs.

# Make gifs with different time windows -----------------------------------
plots <- lapply(graphLists, function(x){
  vultureUtils::plotGraphs(graphList = x, coords = "fixed")
})

map2(.x = plots, .y = intervals, 
     .f = function(.x, .y){
  makeGIF(plotList = .x, fileName = paste0("gif_", stringr::str_replace(.y, " ", ""), ".gif"))
})

# Do some calculations ----------------------------------------------------
dayGraphs <- makeGraphs(edges = southernEdges2021, interval = "1 day", dateTimeStart = "2021-01-01 00:00:00", dateTimeEnd = "2021-12-31 11:59:59", weighted = FALSE)$graphs

# Need to now add all vertices to these graphs.
dayGraphs_allVertices <- makeGraphs(edges = southernEdges2021, interval = "1 days", dateTimeStart = "2021-01-01 00:00:00", dateTimeEnd = "2021-12-31 11:59:59", weighted = FALSE, allVertices = TRUE)$graphs
complete_edgelist <- do.call(igraph::union, dayGraphs_allVertices) %>%
  igraph::get.edgelist()
els <- lapply(dayGraphs_allVertices, igraph::get.edgelist)

# get trues and falses for each edge
test <- lapply(els, function(x){
  complete_edgelist %in% x
})

# bind into a data frame showing presence/absence of edges over time.
overTime <- do.call(cbind, test) %>% as.data.frame()

# name the columns
names(overTime) <- paste0("step", stringr::str_pad(1:length(els), width = 3, side = "left", pad = "0"))

# add two blank steps before, and the edges
beforeSteps <- data.frame(stepPrevPrev = rep(FALSE, nrow(overTime)), 
                             stepPrev = rep(FALSE, nrow(overTime)))
overTime <- cbind(setNames(as.data.frame(complete_edgelist), c("ID1", "ID2")), beforeSteps, overTime)

# Okay now we have complete history. 
# Let's write a function to turn that history into a matrix of h00, h01, etc.
computeProbs <- function(df, startCol = 5){
  histdf <- data.frame("add00" = NA, "add10" = NA, "lose01" = NA, "lose11" = NA)
  for(i in startCol:ncol(df)){
    vec <- vector(mode = "character", nrow(df))
    vec[which(!df[,i-2] & !df[,i-1])] <- "hist00"
    vec[which(!df[,i-2] & df[,i-1])] <- "hist01"
    vec[which(df[,i-2] & !df[,i-1])] <- "hist10"
    vec[which(df[,i-2] & df[,i-1])] <- "hist11"
    
    add00 <- sum(df[i] & vec == "hist00")/sum(vec == "hist00")
    add10 <- sum(df[i] & vec == "hist10")/sum(vec == "hist10")
    lose01 <- sum(!df[i] & vec == "hist01")/sum(vec == "hist01")
    lose11 <- sum(!df[i] & vec == "hist11")/sum(vec == "hist11")
    
    histdf[i-(startCol-1),] <- c("add00" = add00, "add10" = add10, "lose01" = lose01, "lose11" = lose11)
  }
  
  histdfLong <- histdf %>%
    mutate(earlyDate = names(df)[-1:-4]) %>%
    pivot_longer(cols = -earlyDate, names_to = "type", values_to = "prob")
  
  return(histdfLong)
}

test <- computeProbs(df = overTime)

# probability distributions
test %>%
  ggplot(aes(x = prob))+
  geom_density()+
  facet_wrap(~type)+
  theme_minimal()

# Before putting these distributions into the model, need to do a sensitivity analysis. Let's examine the network in increments of 5 days, ranging from 1 day to 30 days.
interval.num <- seq(1, 31, by = 5)
interval <- paste(interval.num, "days")

# run the loop
histdfs <- vector(mode = "list", length = length(interval))
for(i in 1:length(interval)){
  # make the graphs
  graphs <- makeGraphs(edges = southernEdges2021, interval = interval[i], 
                       dateTimeStart = "2021-01-01 00:00:00",
                       dateTimeEnd = "2021-12-31 11:59:00",
                       weighted = FALSE, allVertices = TRUE)$graphs
  complete_edgelist <- do.call(igraph::union, graphs) %>%
    igraph::get.edgelist()
  els <- lapply(graphs, igraph::get.edgelist)
  
  # get trues and falses for each edge
  tf <- lapply(els, function(x){
    complete_edgelist %in% x
  })
  
  # bind into a data frame showing presence/absence of edges over time.
  overTime <- do.call(cbind, tf) %>% as.data.frame()
  
  # add two blank steps before, and the edges
  beforeSteps <- data.frame(stepPrevPrev = rep(FALSE, nrow(overTime)), 
                            stepPrev = rep(FALSE, nrow(overTime)))
  overTime <- cbind(setNames(as.data.frame(complete_edgelist), c("ID1", "ID2")), beforeSteps, overTime)
  
  # compute the probabilities and save them.
  probs <- computeProbs(df = overTime)
  
  histdfs[[i]] <- probs
}

# Add the intervals
histdfs <- map2(.x = histdfs, .y = interval, .f = function(.x, .y){
  .x$interval = factor(.y, levels = .y)
  return(.x)
})

# Make a data frame for plotting
sensData <- data.table::rbindlist(histdfs) %>%
  mutate(earlyDate = lubridate::ymd(earlyDate))

# Make two plots: distribution of probabilities, and probabilities over time.
sensData %>%
  ggplot(aes(x = prob, col = interval))+
  geom_density(size = 1)+
  facet_wrap(~type)+
  theme_minimal()+
  scale_color_viridis_d()

# okay looks like we can probably safely model `add00` and `lose11` as exponential distributions and `add10` and `lose01` as uniform distributions.
# Going to use the 6-day window:
sixdays_add00 <- sensData %>%
  filter(interval == "6 days",
         type == "add00")
r_6days_add00 <- MASS::fitdistr(na.omit(sixdays_add00$prob), densfun = "exponential")$estimate

sixdays_lose11 <- sensData %>%
  filter(interval == "6 days",
         type == "lose11")
r_6days_lose11 <- MASS::fitdistr(na.omit(sixdays_lose11$prob), densfun = "exponential")$estimate

sensData %>%
  ggplot(aes(x = earlyDate, y = prob, col = interval))+
  geom_smooth(se = FALSE)+
  facet_wrap(~type)+
  theme_minimal()+
  scale_color_viridis_d()

# Now let's take a look at the network densities and see how they change.
graphs <- lapply(interval, function(x){
  graphs <- makeGraphs(edges = southernEdges2021, interval = x, 
                       dateTimeStart = "2021-01-01 00:00:00",
                       dateTimeEnd = "2021-12-31 11:59:00",
                       weighted = FALSE, allVertices = TRUE)$graphs
  
})

# compile the density information
densities <- map2(.x = graphs, .y = interval, .f = function(.x, .y){
  lapply(.x, igraph::edge_density) %>% 
    unlist() %>% 
    as.data.frame() %>%
    setNames(., "density") %>%
    mutate(earlyDate = row.names(.),
           interval = factor(.y, levels = .y),
           earlyDate = lubridate::ymd(earlyDate))
}) %>% 
  data.table::rbindlist()

# plot the density information
densities %>%
  ggplot(aes(x = earlyDate, y = density, col = interval))+
  geom_smooth(se = FALSE)+
  theme_minimal()+
  scale_color_viridis_d()

# okay, this graph is a bit hard to interpret. But in general it looks like we get a lot more noise as the interval goes up. Best interval is probably 1 day or 6 days. I'm going to choose 1 day because it's biologically meaningful.

# As for what the density actually IS, let's see the distribution of densities for 1 day
densities %>%
  filter(interval == "1 days") %>%
  ggplot(aes(x = density))+
  geom_density()+
  theme_minimal() # huh, this complicates things. When we have a 1-day time window, the reason the density is so consistent is that it's so close to zero, because almost no edges are present at any given time. It's just small groups of individuals feeding together each day.

# What is the density distribution for a 6 day interval?
densities %>%
  filter(interval == "6 days") %>%
  ggplot(aes(x = density))+
  geom_density()+
  theme_minimal() # 6 days is looking like a more reasonable distribution. BUT, all of these assume that we're allowing isolated nodes and that most individuals aren't connected on any given day. That's a fundamental difference from my model as written. Need to go back and figure that out.

# So if we are including all of the individuals in the network, we should set the density pretty low and use a window of 5 or 6 days (maybe 5 for simplicity?)
# e.g. for 6 days:
(mn6days <- densities %>%
    filter(interval == "6 days") %>%
    pull(density) %>%
    mean()) # 0.0425.
