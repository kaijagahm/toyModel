# Creating a co-feeding network--test

# Load packages
library(move)#for downloading data
library(mapproj)
library(spatsoc);library("asnipe") # for working with the SN parts
library(reshape);library(data.table) #for the manual section where i build the SN myself
library(adehabitatLT);
library(moveVis)
library(dplyr)
library(igraph)
library(vultureUtils)
library(lubridate)

# key parameter values
MaxSpeedPermited=120 #in movebank units (m/s) anything above this will be filtered

VulturesToPlotMap=10:15 
DistThreshold=50 #---at what distance two indi are considered interacting
TimeThreshold='10 minutes' # timegroups - temporally overlapping
MinCoocurForValue=2; #minimal number of coocurences for considering a viable pair- 

############# other variables specified##########

MaxSpeedPermited=120 #in movebank units (m/s??) anything above this will be filtered
roostbuffer = 50 # buffer around roosting polygons (in metres)
feedBuff = 100 # buffer around feeding stations (in metres)

################### reading data from movebank ###################

load("movebankCredentials/pw.Rda")
MB.LoginObject=movebankLogin(username='kaijagahm',password=pw)
rm(pw)

# Get the 2021 data and save it XXX will need to make this a separate script.
# data2021 <- vultureUtils::downloadVultures(loginObject = MB.LoginObject, extraSensors = F, removeDup = T, dateTimeStartUTC = as.POSIXct("2021-01-01 00:00"), dateTimeEndUTC = as.POSIXct("2021-12-31 11:59"))
# save to data/ folder as a .Rda
# save(data2021, file = "data/data2021.Rda")
# Load the pre-downloaded 2021 data
load("data/data2021.Rda")

# convert to a data frame
datDF <- methods::as(data2021, "data.frame")

# Stage 1 filtering
tagsMetadata <- setNames(data.frame(matrix(ncol = 8, nrow = length(datDF))),
                         c("name", "initialPoint", "percentRemoved", "nLocs", "trackDurationDays", "daysBetweenStartEnd", "firstDay", "lastDay"))

datDF$dateOnly <- as.Date(as.character(datDF$timestamp))

# Remove out of Israel locations: select only points falling in Israel
datDFIsrael <- vultureUtils::maskIsrael(dataset = datDF, crs = "WGS84")

# Remove vultures that have fewer than 1/3 of their duration recorded inside Israel
longEnoughIndivs <- vultureUtils::mostlyInIsrael(dataset = datDF, israelDataset = datDFIsrael, thresh = 0.33, dateCol = "dateOnly")

datDF <- datDF %>%
  filter(trackId %in% longEnoughIndivs)

#######################################################################
# use df2move to convert the data.frame into a moveStack
projection(MoveStackDatasetOhad_2021)
Movestacked_BrCoFeed<-df2move(Dataset_BrCoFeed_MinDays,
                              proj = projection(MoveStackDatasetOhad_2021), 
                              x = "coords.x1", y = "coords.x2", time = "timestamps", track_id = "trackId",
                              data = Dataset_BrCoFeed_MinDays)


extent(Movestacked_BrCoFeed)

class(Movestacked_BrCoFeed)
str(Movestacked_BrCoFeed)
summary(Movestacked_BrCoFeed)
colnames(Dataset_BrCoFeed_MinDays)
unstacked_BrCoFeed<-split(Movestacked_BrCoFeed)

length(unstacked_BrCoFeed)

## a loop on tags. ###### #
for (indv in 1:length(unstacked_BrCoFeed) ){## loop on individuals, now in separate Move objects
  TagsMetaData$InitialPoints[indv]=  dim(unstacked_BrCoFeed[[indv]]@data)[1]  ##number of fixes for an individual
  
  TagsMetaData$name[indv]=          names(unstacked_BrCoFeed)[indv]
  plot(unstacked_BrCoFeed[[indv]],col='blue', type='b',main=paste("Indiv=",indv, ', name=',TagsMetaData$name[indv],sep=' '))#just simple plot of this individual's points
  #dim(unstacked_BrCoFeed[[indv]]@coords)
  
  ## removing unneeded columns
  VarsToRemove <- names(unstacked_BrCoFeed[[indv]]@data) %in% c("sensor_type_id","taxon_canonical_name","nick_name","earliest_date_born","sensor","optional",
                                                                "sensor_type","mw_activity_count","eobs_accelerations_raw","eobs_acceleration_sampling_frequency_per_axis",
                                                                "eobs_acceleration_axes","argos_valid_location_algorithm","argos_sensor_4","argos_sensor_3","argos_sensor_2",
                                                                "argos_sensor_1","argos_semi_minor","argos_semi_major","argos_pass_duration","argos_orientation","argos_nopc",
                                                                "argos_lat1","argos_lat2","1084088","argos_lon1","argos_lon2","argos_nb_mes","argos_nb_mes_120",
                                                                "eobs_key_bin_checksum","eobs_fix_battery_voltage","eobs_battery_voltage","eobs_status",
                                                                "eobs_start_timestamp","eobs_type_of_fix","eobs_used_time_to_get_fix","eobs_temperature",
                                                                "gps_dop","magnetic_field_raw_x","magnetic_field_raw_y","magnetic_field_raw_z","ornitela_transmission_protocol",
                                                                "tag_voltage","algorithm_marked_outlier","argos_altitude","argos_best_level","argos_lc","argos_iq",
                                                                "argos_gdop","argos_error_radius","argos_calcul_freq","location_lat.1","location_long.1","timestamps","height_raw",
                                                                "barometric_pressure","barometric_height","battery_charging_current","eobs_activity","manually_marked_outlier",
                                                                "eobs_activity_samples", "acceleration_raw_y", "battery_charge_percent", "data_decoding_software","gps_vdop","height_above_ellipsoid",
                                                                'acceleration_raw_x','acceleration_raw_z',"acceleration_raw_z","eobs_horizontal_accuracy_estimate","eobs_speed_accuracy_estimate");  
  
  
  unstacked_BrCoFeed[[indv]]@data=unstacked_BrCoFeed[[indv]]@data[!VarsToRemove]  ##all columns except above listed ones
  dim(unstacked_BrCoFeed[[indv]]@data)#checking if colunms removed
  
  ## filtering: choosing indices to keep for this individual
  
  indx=1:dim(unstacked_BrCoFeed[[indv]]@data)[1] #starting with all points a
  if(sum(unstacked_BrCoFeed[[indv]]@data$heading <= 360,na.rm=T)){#do i have heading data or this one?
    
    indx=intersect(indx,which(unstacked_BrCoFeed[[indv]]@data$heading <= 360))} #if yes, now index include only points with realistic heading 
  
  if(sum(unstacked_BrCoFeed[[indv]]@data$ground_speed<=MaxSpeedPermited,na.rm=T)){#below threshhold speed?
    indx=intersect(indx,which(unstacked_BrCoFeed[[indv]]@data$ground_speed<=120))}
  
  if(sum(unstacked_BrCoFeed[[indv]]@data$gps_satellite_count>=3,na.rm=T)){#enough satellite numbers?
    indx=intersect(indx,which(unstacked_BrCoFeed[[indv]]@data$gps_satellite_count>=3))}
  
  
  ## subsetting the different slots of this move object
  print(paste("indiv",indv,"name",TagsMetaData$name[indv],'. I throw out', TagsMetaData$InitialPoints[indv]-length(indx), 'points, out of',TagsMetaData$InitialPoints[indv]))
  TagsMetaData$PercentThrown[indv]=(TagsMetaData$InitialPoints[indv]-length(indx))/TagsMetaData$InitialPoints[indv]
  
  unstacked_BrCoFeed[[indv]]@timestamps=unstacked_BrCoFeed[[indv]]@timestamps[indx]
  unstacked_BrCoFeed[[indv]]@sensor=unstacked_BrCoFeed[[indv]]@sensor[indx]
  unstacked_BrCoFeed[[indv]]@data=unstacked_BrCoFeed[[indv]]@data[indx,]
  if(dim(unstacked_BrCoFeed[[indv]]@data)[1]>1){ ##Nitika = to avoid stalling the loop if not enough data locations
    unstacked_BrCoFeed[[indv]]@coords=unstacked_BrCoFeed[[indv]]@coords[indx,]
    unstacked_BrCoFeed[[indv]]@bbox[1,]=range(unstacked_BrCoFeed[[indv]]@coords[,1]);unstacked_BrCoFeed[[indv]]@bbox[2,]=range(unstacked_BrCoFeed[[indv]]@coords[,2])
    
    
    ## collecting metadata and plotting fitered track: 
    TagsMetaData$N_locs[indv]=  dim(unstacked_BrCoFeed[[indv]]@data)[1]
    TagsMetaData$TrackDurationDays[indv]=  length(unique(as.Date(as.character(unstacked_BrCoFeed[[indv]]@data$timestamp))))
    TagsMetaData$FirstDay[indv]=as.character(min(as.Date(as.character(unstacked_BrCoFeed[[indv]]@data$timestamp))));
    TagsMetaData$LastDay[indv]= as.character(max(as.Date(as.character(unstacked_BrCoFeed[[indv]]@data$timestamp))));
    TagsMetaData$DaysBetweenStartEnd[indv]=as.Date(TagsMetaData$LastDay[indv])-as.Date(TagsMetaData$FirstDay[indv]);
    lines(unstacked_BrCoFeed[[indv]],col='red')
    #plot(unstacked_BrCoFeed[[indv]], type="o", col=3, lwd=2, pch=20, xlab="location_long", ylab="location_lat")
    
    ##logging metadata
    
    head(timeLag(unstacked_BrCoFeed[[indv]], units="mins"))
    head(timestamps(unstacked_BrCoFeed[[indv]]))
  }
}#loop on individuals

################## making coords of feeding locations consistent with UTM   ######

BackStacked_BrCoFeed2021=moveStack(unstacked_BrCoFeed)

Dataset2021=as.data.frame(BackStacked_BrCoFeed2021) #converting to a dataset with all indis and their rows

################## getting all non-flying locations ####################

#############   only 2021 datapoints that are at a ground_speed<5 
## check for each date which non-flight datapoints fell inside the food polygons simulataneously on a daily basis
NonFlyingPts <- subset(Dataset2021, heading <= 360 & gps_time_to_fix<=89 & ground_speed<=5 & gps_satellite_count>=3  )


##############    Remove last locations that fall outside Israel cut-off boundary##########

setwd("C:/XXXXX/KML_Files")
FileName='CutOffRegion.kml'
LayerName='Regional_polygon'
Outline = readOGR(dsn="C:/XXXXX/KML_Files/CutOffRegion.kml", 
                  layer="CutOffRegion.kml")
plot(Outline)
Outline@polygons  ##  WGS84

####  convert all Last Locations to a spatial points object

xy <- NonFlyingPts[,c("location_long","location_lat")]
Sp_NonFlyingPts_df <- SpatialPointsDataFrame(coords = xy, data = NonFlyingPts,
                                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
NonFlyingPts_df_IN<-Sp_NonFlyingPts_df[complete.cases(over(Sp_NonFlyingPts_df, Outline)), ]
plot(Outline)
plot(NonFlyingPts_df_IN, col="red", add=TRUE)
#####################Import feeding locations  ####################
setwd("C:/XXXX/")
FoodSites = read.csv("FeedingSites_AllActiveSouthNorth.csv") 
summary(FoodSites)
################## making coords of feeding locations consistent with UTM   ######
BackStacked_BrCoFeed2021=moveStack(unstacked_BrCoFeed)
#preparing a dataframe for mapping:
Dataset2021=as.data.frame(BackStacked_BrCoFeed2021) #converting to a dataset with all indis and their rows
tail(Dataset2021)
Dataset2021$DateOnly=as.Date(as.character(Dataset2021$timestamp));
##for last locations of vultures on each night:
projection(BackStacked_BrCoFeed2021)#this was the original data projection from movebank
class(FoodSites)
colnames(FoodSites)[6:7]<-c("WGS_lat","WGS_long")
FoodSitesF_wgs=FoodSites;
coordinates(FoodSitesF_wgs)<-~WGS_long+WGS_lat
proj4string(FoodSitesF_wgs)<-CRS(projection("+proj=longlat +datum=WGS84 +no_defs"))
##converting to UTM
utmS <- '+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs' #south, but most points are in the north

utmN <- '+proj=utm +zone=36        +ellps=WGS84 +datum=WGS84 +units=m +no_defs'  #north

## converting to UTM36North, 
FoodSitesF_utmN_Poly <- spTransform(FoodSitesF_wgs, CRS(utmN))
head(coordinates(FoodSitesF_utmN_Poly))#now the lat long are in metric
# just plotting the new dataset to see the locations look fine: plot(Dataset2021_utmN,col='blue')
## appending the new coordinates in easting northing- for calculating distance UTM locally
FoodSites$Easting_poly=coordinates(FoodSitesF_utmN_Poly)[,1]
FoodSites$Northing_poly=coordinates(FoodSitesF_utmN_Poly)[,2]
######    creating a circular buffer of radius = 100m #############
coords <- FoodSites[ , c("Easting_poly", "Northing_poly")]   # coordinates
data   <- FoodSites[ , 1:4]          # data
# make the SpatialPointsDataFrame object
spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = data, 
                               proj4string = CRS(utmN))

## make spatial points
dat_sf <- st_as_sf(FoodSites, coords = c("Easting_poly", "Northing_poly"), crs = utmN) 
FoodPolys = st_buffer(dat_sf, dist = feedBuff)
FoodPolys$geometry
nrow(FoodSites)
##add buffer to spatialPolygonDataFrame in UTM projection:
plot(FoodPolys$geometry[c(1:2)], border='blue')
class(FoodPolys$geometry)
FoodPolys$FoodSite
##Getting coordinates of circular buffer around feeding stations 
FoodBuffCoods = as.data.frame(st_coordinates(FoodPolys))
colnames(FoodBuffCoods)<-c("Buff_long","Buff_lat","L1","Sno")
Sites = as.data.frame(cbind(FoodPolys$Sno, as.character(FoodPolys$FoodSite)))
colnames(Sites)<-c("Sno", "StationName")
allPolysFood = merge(FoodBuffCoods, Sites, by = "Sno")
###projecting food buffer boundaries to wgs just for visualizing
projection(BackStacked_BrCoFeed2021)#this was the original data projection from movebank
allPolysFood_utm=allPolysFood;
coordinates(allPolysFood_utm)<-~Buff_long+Buff_lat
proj4string(allPolysFood_utm)<-CRS(projection(utmN))
allPolysFood_wgs = spTransform(allPolysFood_utm, CRS(projection("+proj=longlat +datum=WGS84 +no_defs"))) ##transformed utm to wgs
## appending the new coordinates in easting northing- for calculating distance UTM locally
allPolysFood$wgsLong=coordinates(allPolysFood_wgs)[,1]
allPolysFood$wgsLat=coordinates(allPolysFood_wgs)[,2]


#### Importing all roost polygons ######
setwd("C:/XXXXX/KML_Files")
FileName='AllRoostPolygons.kml'
LayerName='Roosting'
Roosts = readOGR(dsn=FileName, layer=LayerName)
head(Roosts)
summary(Roosts)
allRoostCoods = data.frame()
for (i in 1:length(Roosts)){
  roost1 = as.data.frame(cbind((as.character(Roosts@data[i,1])),Roosts@polygons[[i]]@Polygons[[1]]@coords))
  allRoostCoods = rbind(allRoostCoods, roost1)
}
colnames(allRoostCoods)<-c("Site","x","y") 
allRoostCoods<-allRoostCoods[,c("x","y","Site")]
allRoostCoods$x<-as.numeric(as.character(allRoostCoods$x))
allRoostCoods$y<-as.numeric(as.character(allRoostCoods$y))

allFoodCoods<-as.data.frame(allPolysFood_wgs)[,c("Buff_long","Buff_lat","StationName")]
colnames(allFoodCoods)<-c("x","y","Site") ##with the huge polygons that should be removed

### joining roost polgons with Food buffer polygons

#############IF ONLY LOOKING AT ALL NON-FLIGHT INTERACTIONS EXCEPT WITHIN ROOSTS######
localpolydf<-as.data.frame(allRoostCoods)
localpoly <- localpolydf %>% split(localpolydf$Site) %>% 
  lapply(function(x) rbind(x,x[1,])) %>%
  lapply(function(x) x[,1:2]) %>%
  lapply(function(x) list(as.matrix(x))) %>%
  lapply(function(x) st_polygon(x))

###points
offsetdf<-as.data.frame(NonFlyingPts_df_IN)[,c("coords.x1","coords.x2","trackId")]
points <- st_as_sf(offsetdf,coords=c('coords.x1','coords.x2'),remove = F)
#convert polygons to sf object and add id column
polys <- localpoly %>% st_sfc() %>% st_sf(geom=.) %>% 
  dplyr::mutate(Stn=factor(unique(localpolydf$Site))) 
#find intersection
# Joined is the points falling inside roost polygons
library(sf)
joined <- polys  %>% st_intersection(points) 
#### selecting points falling OUTSIDE food buffers:
#st_intersects which returns a list of the same length as nc_point. If the point isn't in a polygon then it returns integer(0) elements.
#So you can any of a few ways of turning that output into TRUE/FALSE to then select your points. For example, look for length-0 elements in the list:
NonFlyOutRoostFood<-sapply(st_intersects(points, polys),function(x){length(x)==0})
FeedingOutRoostOnly<-points[NonFlyOutRoostFood,]
nrow(FeedingOutRoostOnly %>% st_drop_geometry())
nrow(FeedingOutRoostOnly %>% st_drop_geometry())-nrow(points %>% st_drop_geometry())

setwd("C:/XXXXX")
write.csv(FeedingOutRoostOnly, "Breeding2021_NonFlyingPtsOutsideRoostOnly.csv")

#######################   USING THESE outside roost and outside food points to check spatial proximity and deduce interactions##################
####Make sure to include the seasonal non-flying out-roost points only #####
setwd("C:/XXXXX")
FeedingOutRoostOnly<-read.csv("Breeding2021_NonFlyingPtsOutsideRoostOnly.csv")
Dataset<-FeedingOutRoostOnly
setDT(Dataset, keep.rownames = TRUE)[]
colnames(Dataset)<-c("Row","coords.x1","coords.x2", "trackId", "geometry1","geometry2")

NonFlyingPts_IN<-as.data.frame(NonFlyingPts_df_IN)
setDT(NonFlyingPts_IN, keep.rownames = TRUE)[]
colnames(NonFlyingPts_IN)[1]<-"Row"
head(NonFlyingPts_IN)
Dataset_NonFlightF<-dplyr::left_join(Dataset, NonFlyingPts_IN, by = c("Row"))

############    pruned dataset within Israel to move object
projection(BackStacked_BrCoFeed2021)
BackStacked_NonFlight = df2move(Dataset_NonFlightF,
                                proj = CRS("+proj=longlat +datum=WGS84 +no_defs"), 
                                x = "coords.x1.x", y = "coords.x2.x", time = "timestamps", track_id = "trackId.x")
Dataset_NonFlightF_ltraj=move2ade(BackStacked_NonFlight)
Dataset_NonFlightF=as.data.frame(BackStacked_NonFlight) #converting to a dataset F for filtered
Dataset_NonFlightF=data.table::setDT(Dataset_NonFlightF,keep.rownames=T)#and now to a data.table
## renaming\adding columns
Dataset_NonFlightF$ID=Dataset_NonFlightF$trackId
Dataset_NonFlightF$location_long1=Dataset_NonFlightF$coords.x1
Dataset_NonFlightF$location_lat1=Dataset_NonFlightF$coords.x2
## adding coordinates with a metric value ##### #
## first setting the same coordinate system as in movebank
projection(BackStacked_NonFlight)#this was the original data projection from movebank
Dataset_NonFlightF_wgs=Dataset_NonFlightF;
coordinates(Dataset_NonFlightF_wgs)<-~coords.x1+coords.x2
proj4string(Dataset_NonFlightF_wgs)<-CRS(projection(BackStacked_NonFlight))
utmS <- '+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs' #south, but most points are in the north
utmN <- '+proj=utm +zone=36        +ellps=WGS84 +datum=WGS84 +units=m +no_defs'  #north

## converting to UTM36North, 
Dataset_NonFlightF_utmN <- spTransform(Dataset_NonFlightF_wgs, CRS(utmN))
head(coordinates(Dataset_NonFlightF_utmN))#now the lat long are in metric
# just plotting the new dataset to see the locations look fine: plot(Dataset_NonFlightF_utmN,col='blue')

## appending the new coordinates in easting northing- for calculating distance UTM locally
Dataset_NonFlightF$Easting=coordinates(Dataset_NonFlightF_utmN)[,1]
Dataset_NonFlightF$Northing=coordinates(Dataset_NonFlightF_utmN)[,2]
##spatsoc- timegroups
start_time <- Sys.time()
## using spatsoc for grouping into time groups
library(spatsoc)
start_time <- Sys.time()
group_times(Dataset_NonFlightF, datetime = 'timestamps', threshold = TimeThreshold) 
##group all datapoints that fall within 10 minutes of each other 
## using spatsoc for grouping into spatial groups with the metric values --- using northing and easting
group_pts(Dataset_NonFlightF, threshold = DistThreshold, id = 'ID', coords = c('Northing', 'Easting'), timegroup = 'timegroup')
end_time <- Sys.time()
end_time - start_time

##### manual distance calculation and SN construction #####
SimlDataPntCnt   = expand.grid(unique(as.character(Dataset_NonFlightF$ID)),unique(as.character(Dataset_NonFlightF$ID)))
#a long form of all possible dyads to count interaction
SimlDataPntCnt$counter=0;names(SimlDataPntCnt)=c('ID','ID2','counter') ##setting empty matrix of dimensions pf ID1*ID2
##counter is to see if they co-occur
CoocurCountr = SimlDataPntCnt
#a long form of all possible dyads to count intervals both were at the same timegroup
#creating empty matrix for SRI
SRIlongform=CoocurCountr;names(SRIlongform)[3]='SRI'#--- exhaustive pair-wise
subset(SRIlongform, SRI != 0)
start_time <- Sys.time()
matrices = createDirectedMatrices(Dataset_NonFlightF, DistThreshold)
end_time <- Sys.time()
end_time - start_time
SimlDataPntCnt = matrices[[1]]
CoocurCountr = matrices [[2]]
CoocurCountr$counter = ifelse(CoocurCountr$counter<MinCoocurForValue, 0, CoocurCountr$counter) #eliminating dyads which co-occured less than the minimum cut-off
# Time difference of  mins
##    SRIlongform = actual co-occurence / opportunities of co-occurences
##    SimlDataPntCnt = opportunities of co-occurences
##    CoocurCountr = actual co-occurence

##order of ID1 & ID2 shouldn't matter because undirectional interactions
library(igraph)
SimlDataPntCnt_ordered = get.data.frame(graph.data.frame(SimlDataPntCnt[,c('ID','ID2')], directed=FALSE),"edges") 
SimlDataPntCnt$ID=SimlDataPntCnt_ordered[,1]
SimlDataPntCnt$ID2=SimlDataPntCnt_ordered[,2]
### I did not aggregate the values of co=occurence because in a given time when A and B co-occured it doesn't matter for us is A co-occured with B or B co-occured with A
SimlDataPntCnt=unique(SimlDataPntCnt)
CoocurCountr_ordered = get.data.frame(graph.data.frame(CoocurCountr[,c('ID','ID2')], directed=FALSE),"edges") 
CoocurCountr$ID=CoocurCountr_ordered[,1]
CoocurCountr$ID2=CoocurCountr_ordered[,2]
CoocurCountr=unique(CoocurCountr)
##order of ID1 & ID2 shouldn't matter
library(igraph)
nrow(SRIlongform)
SRIlongform_ordered = get.data.frame(graph.data.frame(SRIlongform[,c('ID','ID2')], directed=FALSE),"edges") 
SRIlongform$ID=SRIlongform_ordered[,1]
SRIlongform$ID2=SRIlongform_ordered[,2]
SRIlongform=unique(SRIlongform)
SRIlongform$SRI=as.numeric(CoocurCountr$counter/SimlDataPntCnt$counter)# ratio of number of co-occurances/number of simluatnous datapoints 
#SRIlongform has exhaustive dyads which may or may not have had the opportunity to co-occur
hist(SRIlongform$SRI);range(SRIlongform$SRI,na.rm=T)
SRI_mrtx=as.matrix(tidyr::spread(data=SRIlongform, key= ID2, value=SRI)) ##ID2 is key meaning it becomes column names in the matrix
Coocur_mrtx=as.matrix(tidyr::spread(data=CoocurCountr, key= ID2, value=counter))
NamesInMtrx=SRI_mrtx[,1];SRI_mrtx=SRI_mrtx[,-1]#just setting row names from the dataframe
M=mapply(SRI_mrtx, FUN=as.numeric)#converting to numeric matrix
SRI_mrtx<- matrix(data=M, ncol=ncol(SRI_mrtx), nrow=nrow(SRI_mrtx))
rownames(SRI_mrtx)=NamesInMtrx;colnames(SRI_mrtx)=NamesInMtrx
Diag=diag(SRI_mrtx);unique(Diag);#self, indeed always 1, lets remove them:
diag(SRI_mrtx)=NA
MaxIndSRI=apply(SRI_mrtx,2,max, na.rm=T)#the max value in each row
print(paste('even after removing diagonal self connection there are still', sum(MaxIndSRI==1),'fully connected dyads with MinCoocurForValue of',MinCoocurForValue))