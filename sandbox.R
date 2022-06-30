# Sandbox script to play around with the feeding network data

# load packages
library(vultureUtils)
library(tidyverse)

# load data
load("data/feedingEdges2021.Rda")
load("data/feedingPoints2021.Rda")

unweighted <- vultureUtils::makeGraphs(edges = feedingEdges2021, fullData = feedingPoints2021, interval = "20 days", dateTimeStart = "2021-01-01 00:00:00", dateTimeEnd = "2021-12-31 11:59:59", weighted = FALSE)$graphs

# Plot the unweighted networks
plots <- vultureUtils::plotGraphs(unweighted)

# Make an animated gif
vultureUtils::makeGIF(plots, fileName = "2021_20Days.gif", interval = 0.2)


# Examine network metrics vs. time window ---------------------------------
# Let's look at e.g. network density
timeWindows <- seq(1, 365/2, by = 10)
intervals <- paste(timeWindows, "days")
dts <- "2021-01-01 00:00:00"
dte <- "2021-12-31 11:59:59"
e <- feedingEdges2021
f <- feedingPoints2021

dat <- vector(mode = "list", length = length(intervals))
names(dat) <- intervals
for(i in 1:length(graphLists)){
  cat(paste("Computing graphs for an interval of", intervals[i], "\n"))
  dat[[i]] <- makeGraphs(edges = e, fullData = f,
                    interval = intervals[i], dateTimeStart = dts,
                    dateTimeEnd = dte, weighted = FALSE)
}

graphLists <- lapply(dat, function(x){x[["graphs"]]})
breaks <- lapply(dat, function(x){x[["breaks"]]})
names(breaks) <- intervals
names(graphLists) <- intervals

densities <- lapply(graphLists, function(x){
  lapply(x, igraph::edge_density)
})
names(densities) <- intervals

# Make a data frame for plotting
together <- map2(.x = breaks, .y = densities, .f = function(.x, .y){
  data.frame("earlyDate" = .x,
             "networkDensity" = .y)
})




