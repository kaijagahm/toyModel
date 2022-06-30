# Sandbox script to play around with the feeding network data

# load packages
library(vultureUtils)
library(tidyverse)

# load data
load("data/feedingEdges2021.Rda")
load("data/feedingPoints2021.Rda")

unweighted <- makeGraphs(edges = feedingEdges2021, fullData = feedingPoints2021, interval = "20 days", dateTimeStart = "2021-01-01 00:00:00", dateTimeEnd = "2021-12-31 11:59:59", weighted = FALSE)$graphs

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
for(i in 1:length(dat)){
  cat(paste("Computing graphs for an interval of", intervals[i], "\n"))
  dat[[i]] <- makeGraphs(edges = e, fullData = f,
                    interval = intervals[i], dateTimeStart = dts,
                    dateTimeEnd = dte, weighted = FALSE)
}

graphLists <- lapply(dat, function(x){x[["graphs"]]})
names(graphLists) <- intervals

densities <- lapply(graphLists, function(x){
  unlist(lapply(x, igraph::edge_density)) %>%
    as.data.frame() %>%
    setNames("networkDensity") %>%
    mutate(earlyDate = row.names(.))
})
names(densities) <- intervals

df <- data.table::rbindlist(densities, idcol = TRUE) %>%
  as.data.frame() %>%
  rename("interval" = .id) %>%
  mutate(interval = stringr::str_replace(interval, " days", ""),
         interval = as.numeric(interval)) %>%
  mutate(temp = lag(earlyDate)) %>%
  mutate(earlyDate = case_when(is.na(earlyDate) ~ as.character(as.Date(temp) + as.numeric(interval)),
                               TRUE ~ earlyDate)) %>%
  select(-temp) # XXX SHOULD DO THIS FILLING IN IN THE FUNCTION ITSELF-- ADD ONE MORE TO BREAKS.

df %>%
  ggplot(aes(x = as.Date(earlyDate), y = networkDensity, col = as.factor(interval)))+
  geom_smooth(se = F, method = "lm")
# this pattern could be an artifact of the last group being fewer days than the rest. Remove the last group.

dfRemoved <- df %>%
  group_by(interval) %>%
  slice(1:(n()-1))

dfRemoved %>%
  ggplot(aes(x = as.Date(earlyDate), y = networkDensity, col = as.factor(interval)))+
  geom_smooth(se = F, method = "lm")
# Huh, still see a steep downward slope for a few of them. But the bigger picture is that there is not much pattern here--except for the last few, it's quite flat.


