# Creating a co-feeding network--test

# Load packages
library(move)#for downloading data
library(spatsoc)
library(data.table) #for the manual section where i build the SN myself
library(moveVis)
library(dplyr)
library(igraph)
library(vultureUtils)
library(lubridate)

# key parameter values
DistThreshold=50 #---at what distance two indi are considered interacting
TimeThreshold='10 minutes' # timegroups - temporally overlapping
coOccurrenceThreshold <- 2 #minimal number of coocurences for considering a viable pair- 
inIsraelThreshold <- 0.33 # proportion of days tracked that must fall in Israel

############# other variables specified##########

roostBuffer <- 50 # buffer around roosting polygons (in metres)
feedingBuffer <- 100 # buffer around feeding stations (in metres)

################### reading data from movebank ###################

load("movebankCredentials/pw.Rda")
MB.LoginObject=movebankLogin(username='kaijagahm',password=pw)
rm(pw)

# Get the 2021 data and save it XXX will need to make this a separate script.
# data2021 <- vultureUtils::downloadVultures(loginObject = MB.LoginObject, extraSensors = F, removeDup = T, dateTimeStartUTC = as.POSIXct("2021-01-01 00:00"), dateTimeEndUTC = as.POSIXct("2021-12-31 11:59"))
#testData2021 <- vultureUtils::downloadVultures(loginObject = MB.LoginObject, extraSensors = F, removeDup = T, dateTimeStartUTC = as.POSIXct("2021-01-01 00:00"), dateTimeEndUTC = as.POSIXct("2021-03-30 11:59"))
#save(testData2021, file = "data/testData2021.Rda")
# save to data/ folder as a .Rda
# save(data2021, file = "data/data2021.Rda")
# Load the pre-downloaded 2021 data
#load("data/testData2021.Rda")
load("data/data2021.Rda")

# convert to a data frame
datDF <- methods::as(data2021, "data.frame")

# remove unnecessary columns
datDF <- vultureUtils::removeUnnecessaryVars(datDF)

datDF$dateOnly <- as.Date(as.character(datDF$timestamp))

# Remove out of Israel locations: select only points falling in Israel
datDFIsrael <- vultureUtils::maskIsrael(dataset = datDF, crs = "WGS84")

# Remove vultures that have fewer than 1/3 of their duration recorded inside Israel
longEnoughIndivs <- vultureUtils::mostlyInIsrael(dataset = datDF, israelDataset = datDFIsrael, thresh = inIsraelThreshold, dateCol = "dateOnly")
datDF <- datDF %>% # using datDF because we don't want to actually restrict it to Israel yet
  filter(trackId %in% longEnoughIndivs)

# Clean the data and extract metadata
## by setting speedUpper to 5, we are restricting this to non-flight interactions.
cleaned <- vultureUtils::cleanAndMetadata(df = datDFIsrael, speedUpper = 5)
filteredData <- cleaned$filteredData
tagsMetadata <- cleaned$tagsMetadata

cleanedIsrael <- vultureUtils::maskIsrael(dataset = filteredData, longCol = "location_long.1",
                                          latCol = "location_lat.1", crs = "WGS84")

# Import feeding site locations and buffer them
feedingSites <- read.csv("data/FeedingSites_AllActiveSouthNorth.csv")
feedingSites <- vultureUtils::bufferFeedingSites(feedingSites, crsToSet = "WGS84", crsToReturn = "WGS84")

# Import roost polygons
roostPolygons <- sf::st_read("data/AllRoostPolygons.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")

# Exclude any points that fall within a roost polygon
feedingPoints <- cleanedIsrael[lengths(sf::st_intersects(cleanedIsrael, roostPolygons)) == 0,]

# Convert to UTM
feedingPoints <- feedingPoints %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1], # save lat/long coordinates, just in case we need them later.
                lat = sf::st_coordinates(.)[,2]) %>%
  # convert to UTM (this will be useful for calculating distance locally.)
  sf::st_transform(32636) %>%
  dplyr::mutate(utmE = sf::st_coordinates(.)[,1], # get utm coords as separate columns.
                utmN = sf::st_coordinates(.)[,2]) %>%
  sf::st_drop_geometry() # drop geometry, making this just a data frame, so that spatsoc will work.

# Convert the `timestamp` column to POSIXct
feedingPoints <- feedingPoints %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))
data.table::setDT(feedingPoints) # convert to a data.table so we can proceed with spatsoc

# Group these points into timegroups using `spatsoc`
spatsoc::group_times(feedingPoints, datetime = 'timestamp', threshold = TimeThreshold) 

# Group into point groups
spatsoc::group_pts(feedingPoints, threshold = DistThreshold, id = 'trackId', 
          coords = c('utmE', 'utmN'), timegroup = 'timegroup')

# Generate edge lists by timegroup
edges <- edge_dist(DT = feedingPoints, threshold = DistThreshold, id = "trackId",
                   coords = c('utmE', "utmN"), timegroup = "timegroup",
                   returnDist = TRUE, fillNA = FALSE)

# Remove duplicates
edges <- edges %>%
  filter(ID1 < ID2)

# Now create a list where the edge only stays if it occurred in at least `coOccurrenceThreshold` consecutive time steps
edges <- vultureUtils::consecEdges(edgeList = edges, consecThreshold = coOccurrenceThreshold) 
# XXX this still doesn't work to remove the grouping variable, but I just can't be bothered to fix it right now because it's such a pain. Come back to this.

save(edges, file = "data/feedingEdges2021.Rda")