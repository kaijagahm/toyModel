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

feedingEdges2021 <- spaceTimeGroups(dataset = feedingPoints, distThreshold = distThreshold, consecThreshold = consecThreshold)

feedingPoints2021 <- feedingPoints

save(feedingEdges2021, file = "data/feedingEdges2021.Rda")
save(feedingPoints2021, file = "data/feedingPoints2021.Rda")
