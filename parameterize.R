# AUTHOR: Kaija Gahm
# CREATED: 2022-06-21

# Load packages
library(igraph)
library(tidyverse)
library(vultureUtils)
library(sna)
library(spatsoc)
library(move)

# Create movebank login object
load("movebankCredentials/pw.Rda")
MB.LoginObject <- movebankLogin(username = 'kaijagahm', 
                                password = pw)

# Download 2021 vulture data (assuming that it's in UTC already--need to check this) XXX
dat <- downloadVultures(loginObject = MB.LoginObject, dateTimeStartUTC = as.POSIXct("2021-01-01 00:00"), dateTimeEndUTC = as.POSIXct("2021-12-31 11:59"))

# Create a co-feeding network for the whole year

# Create a co-feeding network for each month

# Create a co-feeding network for each day