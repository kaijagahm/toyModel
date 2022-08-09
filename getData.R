# Script to download and format/subset data, using functions from vultureUtils.
# Created on 2022-08-09 as a replacement for getCoFeedingData.R, getCoFlightData.R, and getCoRoostingData.R. No need for three separate scripts.

# Load packages -----------------------------------------------------------
library(move)#for downloading data
library(spatsoc)
library(data.table) #for the manual section where i build the SN myself
library(moveVis)
library(dplyr)
library(igraph)
library(vultureUtils) # self-written package. Can be installed using devtools::install_github("kaijagahm/vultureUtils").
library(lubridate)

# Movebank login information
load("movebankCredentials/pw.Rda")
MB.LoginObject <- movebankLogin(username = 'kaijagahm', password = pw)

# Download data for the right dates, and save it so we can re-load it later.
#data_20200101_20220430 <- vultureUtils::downloadVultures(loginObject = MB.LoginObject, extraSensors = F, removeDup = T, dateTimeStartUTC = as.POSIXct("2020-01-01 00:00"), dateTimeEndUTC = as.POSIXct("2022-04-30 11:59"))
# save to data/ folder as a .Rda
#save(data_20200101_20220430, file = "data/data_20200101_20220430.Rda")
# Load the pre-downloaded data
load("data/data_20200101_20220430.Rda")

# Co-feeding --------------------------------------------------------------
# convert to a data frame
datDF <- methods::as(data_20200101_20220430, "data.frame")

# remove unnecessary columns
datDF <- vultureUtils::removeUnnecessaryVars(datDF)

datDF$dateOnly <- as.Date(as.character(datDF$timestamp))

# Get Israel mask
israelMask <- sf::st_read("data/maskIsrael.kml")

# Import roost polygons
roostPolygons <- sf::st_read("data/AllRoostPolygons.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")

# Buffer the roosts by the buffer distance
roostPolygons <- vultureUtils::convertAndBuffer(obj = roostPolygons, dist = roostBuffer)

feedingEdges_20200101_20220430 <- vultureUtils::getFeedingEdges(dataset = datDF, mask = israelMask, roostPolygons = roostPolygons, roostBuffer = roostBuffer, inMaskThreshold = inIsraelThreshold, consecThreshold = consecThreshold, distThreshold = distThreshold)

save(flightEdges_20200101_20220430, file = "data/flightEdges_20200101_20220430.Rda")

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

# Co-flight ---------------------------------------------------------------

# convert to a data frame
datDF <- methods::as(data_20200101_20220430, "data.frame")

# remove unnecessary columns
datDF <- vultureUtils::removeUnnecessaryVars(datDF)

datDF$dateOnly <- as.Date(as.character(datDF$timestamp))

# Get Israel mask
israelMask <- sf::st_read("data/maskIsrael.kml")

# Import roost polygons
roostPolygons <- sf::st_read("data/AllRoostPolygons.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")

# Buffer the roosts by the buffer distance
roostPolygons <- vultureUtils::convertAndBuffer(obj = roostPolygons, dist = roostBuffer)

flightEdges_20200101_20220430 <- vultureUtils::getFlightEdges(dataset = datDF, mask = israelMask, roostPolygons = roostPolygons, roostBuffer = roostBuffer, inMaskThreshold = inIsraelThreshold, consecThreshold = consecThreshold, distThreshold = distThreshold)

save(flightEdges_20200101_20220430, file = "data/flightEdges_20200101_20220430.Rda")

# Subset the data to only southern individuals. Properly, I should compute home ranges or something to determine where each individual hangs out. But I'm going to just take individuals' mean latitude for now, and use that as a proxy for their space use, since we know the southern population is generally distinct from the two northern populations. Note that this metric will confound vultures in the Golan and in the Carmel, but since both of those are northern populations that I don't need to deal with, I won't worry about that for now.

# Filter to include only southern individuals
# get mean latitude for each individual
indivs <- flightPoints_20200101_20220430 %>%
  sf::st_drop_geometry() %>%
  dplyr::select(trackId, location_long.1, location_lat.1) %>%
  dplyr::group_by(trackId) %>%
  dplyr::summarize(mnlat = mean(location_lat.1))

# get southern individuals
southern <- indivs %>%
  dplyr::filter(mnlat < 32) %>%
  dplyr::pull(trackId)

# restrict the edgelist to only southern individuals
southernFlightEdges_20200101_20220430 <- flightEdges_20200101_20220430 %>%
  dplyr::filter(ID1 %in% southern & ID2 %in% southern)

# restrict the points to only southern individuals
southernFlightPoints_20200101_20220430 <- flightPoints_20200101_20220430 %>%
  dplyr::filter(trackId %in% southern)

save(southernFlightEdges_20200101_20220430, file = "data/southernFlightEdges_20200101_20220430.Rda")
save(southernFlightPoints_20200101_20220430, file = "data/southernFlightPoints_20200101_20220430.Rda")


# Co-roosting -------------------------------------------------------------
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
## by setting speedThreshUpper to 5, we are restricting this to non-flight interactions.
filteredData <- filterLocs(df = datDF, speedThreshUpper = 5)

# Now mask again to remove the out-of-Israel points.
cleanedIsrael <- vultureUtils::maskData(dataset = filteredData, mask = israelMask, longCol = "location_long.1", latCol = "location_lat.1", crs = "WGS84")

# Import roost polygons
roostPolygons <- sf::st_read("data/AllRoostPolygons.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")

# Buffer roost polygons by the buffer distance:
roostPolygons <- vultureUtils::convertAndBuffer(obj = roostPolygons, dist = roostBuffer)

# Find last locations of each vulture on each night, and first locations of each vulture on each morning.
lastLocs <- cleanedIsrael %>%
  group_by(trackId, dateOnly) %>%
  filter(timestamp == max(timestamp)) %>%
  ungroup() %>%
  dplyr::select(trackId, dateOnly, timestamp)

firstLocs <- cleanedIsrael %>%
  group_by(trackId, dateOnly) %>%
  filter(timestamp == min(timestamp)) %>%
  ungroup() %>%
  dplyr::select(trackId, dateOnly, timestamp)

sleepInRoost <- sf::st_join(lastLocs, roostPolygons, join = sf::st_within) %>%
  mutate(type = "night") %>%
  select(-Description) %>%
  mutate(nightDate = dateOnly) %>%
  group_by(trackId, nightDate) %>%
  slice(1)

wakeInRoost <- sf::st_join(firstLocs, roostPolygons, join = sf::st_within) %>%
  mutate(type = "morning") %>%
  select(-Description) %>%
  mutate(nightDate = dateOnly - 1) %>%
  group_by(trackId, nightDate) %>% # remove any duplicates--just take the first roost assignment
  slice(1)

# Bo's idea: just do all four calculations for all of the bird/nights, and then have the iterative if/elses apply to the choice of pre-calculated values, not to the calculations themselves. This is a bit less computationally efficient, but probably a lot easier to understand.
morningNight <- bind_rows(sleepInRoost, wakeInRoost) %>%
  arrange(trackId, nightDate, type) %>%
  sf::st_drop_geometry() %>%
  tidyr::pivot_wider(id_cols = c(trackId, nightDate), names_from = "type", values_from = "Name")

midpoints <- bind_rows(sleepInRoost, wakeInRoost) %>%
  select(-c(Name, type, timestamp)) %>%
  mutate(lon = sf::st_coordinates(.)[,1],
         lat = sf::st_coordinates(.)[,2]) %>%
  sf::st_drop_geometry() %>%
  group_by(trackId, nightDate) %>%
  summarize(midpointLon = mean(lon),
            midpointLat = mean(lat)) %>% 
  sf::st_as_sf(., coords = c("midpointLon", "midpointLat")) %>%
  sf::st_set_crs("WGS84") %>%
  sf::st_join(., roostPolygons, join = sf::st_within) %>%
  select(-Description) %>%
  rename("midpointRoost" = "Name") %>%
  group_by(trackId, nightDate) %>% # get rid of any duplicates
  slice(1)

nearestRoostIndex <- sf::st_nearest_feature(x = midpoints, y = roostPolygons)
nearestRoost <- roostPolygons$Name[nearestRoostIndex]
midpoints$nearestRoost <- nearestRoost

nrow(midpoints) == nrow(morningNight)

# Join them back together
roosts <- left_join(morningNight, midpoints, by = c("trackId", "nightDate")) %>%
  select(-geometry)
nrow(roosts) == nrow(midpoints) # check that the row count is the same

# Now let's assign a final roost based on the criteria
roosts <- roosts %>%
  mutate(finalRoost = coalesce(night, morning, midpointRoost, nearestRoost))

roostAssignments <- roosts %>%
  select(trackId, nightDate, "roost" = finalRoost)

# Now, convert this to an edge list. XXX CAN'T FIGURE OUT HOW TO DO THIS, GRR
head(roostAssignments)

test <- roostAssignments %>%
  mutate(trackId = as.character(trackId)) %>%
  group_by(nightDate, roost) %>%
  group_split(.keep = TRUE) %>%
  lapply(., function(x){
    indivs <- x$trackId
    nightDate <- x$nightDate[1]
    roost <- x$roost[1]
    if(length(indivs) > 1){
      el <- t(combn(indivs, 2)) %>% as.matrix() %>% 
        as.data.frame() %>%
        mutate(nightDate = nightDate,
               roost = roost)
      return(el)
    }else{
      return(NULL)
    }
  }) %>%
  data.table::rbindlist() %>%
  as.data.frame()
# XXX START HERE--NOT DONE.

