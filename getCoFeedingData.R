# Creating a co-feeding network--test

# Load packages
library(move)#for downloading data
library(spatsoc)
library(data.table) #for the manual section where i build the SN myself
library(moveVis)
library(dplyr)
library(igraph)
library(vultureUtils) # self-written package. Can be installed using devtools::install_github("kaijagahm/vultureUtils").
library(lubridate)

# key parameter values
distThreshold <- 50 # distance at which vultures are considered interacting
consecThreshold <- 2 #minimal number of coocurences for considering a viable pair- 
inIsraelThreshold <- 0.33 # proportion of days tracked that must fall in Israel

############# other variables specified##########

roostBuffer <- 50 # buffer around roosting polygons (in metres)
feedingBuffer <- 100 # buffer around feeding stations (in metres)

################### reading data from movebank ###################

load("movebankCredentials/pw.Rda")
MB.LoginObject=movebankLogin(username='kaijagahm',password=pw)
rm(pw)

# Get the data and save it.
#data_20200101_20220430 <- vultureUtils::downloadVultures(loginObject = MB.LoginObject, extraSensors = F, removeDup = T, dateTimeStartUTC = as.POSIXct("2020-01-01 00:00"), dateTimeEndUTC = as.POSIXct("2022-04-30 11:59"))
# save to data/ folder as a .Rda
#save(data_20200101_20220430, file = "data/data_20200101_20220430.Rda")
# Load the pre-downloaded data
load("data/data_20200101_20220430.Rda"); beepr::beep()

# convert to a data frame
datDF <- methods::as(data_20200101_20220430, "data.frame")

# remove unnecessary columns
datDF <- vultureUtils::removeUnnecessaryVars(datDF)

datDF$dateOnly <- as.Date(as.character(datDF$timestamp))

# Get Israel mask
israelMask <- sf::st_read("data/maskIsrael.kml")

# Remove out of Israel locations: select only points falling in Israel
datDFIsrael <- vultureUtils::maskData(dataset = datDF, mask = israelMask, longCol = "location_long.1", latCol = "location_lat.1", crs = "WGS84")

# Remove vultures that have less than 1/3 of their duration recorded inside Israel
longEnoughIndivs <- vultureUtils::mostlyInMask(dataset = datDF, maskedDataset = datDFIsrael, thresh = inIsraelThreshold, dateCol = "dateOnly")
datDF <- datDF %>% # using datDF because we don't want to actually restrict it to Israel yet
  filter(trackId %in% longEnoughIndivs)

# Filter
## by setting speedUpper to 5, we are restricting this to non-flight interactions.
filteredData <- vultureUtils::filterLocs(df = datDF, speedThreshUpper = 5)

# Now mask again to remove the out-of-Israel points.
cleanedIsrael <- vultureUtils::maskData(dataset = filteredData, mask = israelMask, longCol = "location_long.1", latCol = "location_lat.1", crs = "WGS84")

# Import feeding site locations and buffer them
feedingSites <- read.csv("data/FeedingSites_AllActiveSouthNorth.csv")
feedingSites <- vultureUtils::bufferFeedingSites(feedingSites, crsToSet = "WGS84", crsToReturn = "WGS84")

# Import roost polygons
roostPolygons <- sf::st_read("data/AllRoostPolygons.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")

# Exclude any points that fall within a roost polygon
feedingPoints <- cleanedIsrael[lengths(sf::st_intersects(cleanedIsrael, roostPolygons)) == 0,]

feedingEdges_20200101_20220430 <- spaceTimeGroups(dataset = feedingPoints, distThreshold = distThreshold, consecThreshold = consecThreshold)

feedingPoints_20200101_20220430 <- feedingPoints

save(feedingEdges_20200101_20220430, file = "data/feedingEdges_20200101_20220430.Rda")
save(feedingPoints_20200101_20220430, file = "data/feedingPoints_20200101_20220430.Rda")

# Subset the data to only southern individuals. Properly, I should compute home ranges or something to determine where each individual hangs out. But I'm going to just take individuals' mean latitude for now, and use that as a proxy for their space use, since we know the southern population is generally distinct from the two northern populations. Note that this metric will confound vultures in the Golan and in the Carmel, but since both of those are northern populations that I don't need to deal with, I won't worry about that for now.

# Filter to include only southern individuals
# get mean latitude for each individual
indivs <- feedingPoints_20200101_20220430 %>%
  sf::st_drop_geometry() %>%
  dplyr::select(trackId, location_long.1, location_lat.1) %>%
  dplyr::group_by(trackId) %>%
  dplyr::summarize(mnlat = mean(location_lat.1))

# get southern individuals
southern <- indivs %>%
  dplyr::filter(mnlat < 32) %>%
  dplyr::pull(trackId)

# restrict the edgelist to only southern individuals
southernEdges_20200101_20220430 <- feedingEdges_20200101_20220430 %>%
  dplyr::filter(ID1 %in% southern & ID2 %in% southern)

# restrict the points to only southern individuals
southernPoints_20200101_20220430 <- feedingPoints_20200101_20220430 %>%
  dplyr::filter(trackId %in% southern)

save(southernEdges_20200101_20220430, file = "data/southernEdges_20200101_20220430.Rda")
save(southernPoints_20200101_20220430, file = "data/southernPoints_20200101_20220430.Rda")
beepr::beep()
