# Creating a co-flight network--test

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
consecThreshold <- 2 #minimal number of coocurences for considering a viable pair- 
inIsraelThreshold <- 0.33 # proportion of days tracked that must fall in Israel
roostBuffer <- 50 # buffer around roosting polygons (in metres)

load("movebankCredentials/pw.Rda")
MB.LoginObject=movebankLogin(username='kaijagahm',password=pw)
rm(pw)

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
## by setting speedThreshUpper to 5, we are restricting this to non-flight interactions.
filteredData <- filterLocs(df = datDF, speedThreshUpper = 5)

# Now mask again to remove the out-of-Israel points.
cleanedIsrael <- vultureUtils::maskData(dataset = filteredData, mask = israelMask, longCol = "location_long.1", latCol = "location_lat.1", crs = "WGS84")

# Import roost polygons
roostPolygons <- sf::st_read("data/AllRoostPolygons.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")

# Buffer roost polygons by the buffer distance:
# XXX do this.

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

test <- bind_rows(sleepInRoost, wakeInRoost) %>%
  arrange(trackId, nightDate, type) %>%
  rename("directAssignment" = "Name") %>%
  mutate(lon = sf::st_coordinates(.)[,1],
         lat = sf::st_coordinates(.)[,2]) %>%
  sf::st_drop_geometry() %>%
  group_by(trackId, nightDate) %>%
  mutate(midpointLon = mean(lon),
         midpointLat = mean(lat)) %>% 
  sf::st_as_sf(., coords = c("midpointLon", "midpointLat")) %>%
  sf::st_set_crs("WGS84") %>%
  sf::st_join(., roostPolygons, join = sf::st_within) %>%
  select(-Description) %>%
  rename("midpointRoost" = "Name")

nearestRoostIndex <- sf::st_nearest_feature(x = test, y = roostPolygons)
nearestRoost <- roostPolygons$Name[nearestRoostIndex]
test$nearestRoost <- nearestRoost

# MorningNight
morningNight <- test %>%
  dplyr::select(trackId, nightDate, directAssignment, type) %>%
  tidyr::pivot_wider(id_cols = c("trackId", "nightDate"), names_from = "type", values_from = "directAssignment")


 # XXX this went wrong--figure out what is wrong here!
test2 <- test %>%
  sf::st_drop_geometry() %>%
  dplyr::select(trackId, directAssignment, type, nightDate, midpointRoost, nearestRoost) %>%
  tidyr::pivot_wider(id_cols = c("trackId", "nightDate"), names_from = "type", values_from = "directAssignment", names_prefix = "roost_")
# XXXX end here

# Now either use coalesce with four columns or use case_when or whatever to make a new column for the final roost assignment.
nrow(test) # should be the same as the number of unique vulture/date combinations
bind_rows(sleepInRoost, wakeInRoost) %>%
  sf::st_drop_geometry() %>%
  select(trackId, nightDate) %>%
  distinct() %>%
  nrow()



# Let's create a data frame full of all the possible combinations of vultures and nights.
vultureNights <- bind_rows(sleepInRoost, wakeInRoost) %>%
  sf::st_drop_geometry() %>%
  select(trackId, nightDate) %>% # remove any duplicates--just take the first roost assignment
  distinct()

# Assign a roost based on night, if possible
vultureNights <- vultureNights %>%
  left_join(sleepInRoost %>% 
              sf::st_drop_geometry() %>% 
              select(trackId, nightDate, Name), 
            by = c("trackId", "nightDate"))

# Assign a roost based on morning, if it couldn't be assigned based on the night
vultureNights <- vultureNights %>%
  rows_patch(., wakeInRoost %>%
               sf::st_drop_geometry() %>%
               select(trackId, nightDate, Name),
             by = c("trackId", "nightDate"))

# For those that are still not assigned, get the data
unAssignedNights <- vultureNights %>%
  filter(is.na(Name)) %>%
  select(-Name) %>%
  left_join(sleepInRoost %>%
              select(trackId, nightDate), by = c("trackId", "nightDate")) %>%
  mutate(type = "night") %>%
  sf::st_as_sf() %>%
  filter(!sf::st_is_empty(.))

unAssignedDays <- vultureNights %>%
  filter(is.na(Name)) %>%
  select(-Name) %>%
  left_join(wakeInRoost %>%
              select(trackId, nightDate), by = c("trackId", "nightDate")) %>%
  mutate(type = "morning") %>%
  sf::st_as_sf() %>%
  filter(!sf::st_is_empty(.))

unAssigned <- bind_rows(unAssignedNights, unAssignedDays) %>%
  arrange(trackId, nightDate)

# Get midpoints
averages <- unassigned %>%
  group_by(trackId, nightDate) %>% 
  summarize(geometry = sf::st_union(geometry)) %>% 
  sf::st_centroid()

# See if the midpoints can be assigned to roosts
averagesAssigned <- sf::st_join(averages, roostPolygons, join = sf::st_within)

# Update the Name column for those that could be assigned
vultureNights <- vultureNights %>%
  rows_patch(., averagesAssigned %>%
               sf::st_drop_geometry() %>%
               select(trackId, nightDate, Name),
             by = c("trackId", "nightDate"))

# Get the averages that still can't be assigned
unassignedAverages <- averagesAssigned %>%
  filter(is.na(Name))

nearestFeatures <- sf::st_nearest_feature(x = unassignedAverages, y = roostPolygons)
roostNames <- roostPolygons$Name[nearestFeatures]
unassignedAverages$Name <- roostNames

# Now update vultureNights again:
vultureNights <- vultureNights %>%
  rows_patch(., unassignedAverages %>%
               sf::st_drop_geometry() %>%
               select(trackId, nightDate, Name),
             by = c("trackId", "nightDate"))

# Check if any are still unassigned
vultureNights %>%
  filter(is.na(Name)) # XXX why are these still unassigned???
  
#
# i. Assigned last locations of each vulture on each night to a roost polygon
# ii. If locations unassigned, found the first location for the morning after (location of where a vulture woke up the next day)
# ---> assigned to a roost polygon
# iii. If still unassigned, found the average location of the last location at night and first location the following morning and then
# ----> assigned to a roost polygon
# iv. If owing to the 50m buffer around roost polygons, a location is assigned to >1 roost, then remove duplicates and
# -----> assign the location to the closer roost polygon
#
