# ---
# title: "Translate BBS data"
# author: "Melina Houle"
# date: "July 08, 2026"
# Note on translation:
# -- The script download data and extract the ones not being uploaded in WildTrax. 

library(googlesheets4)
library(googledrive)
library(tidyr)
library(readr)
library(chron)
library(readxl)
library(plyr) #rbind.fill
library(dplyr)
library(stringr) #str_pad
 
##USGS URL
#url <- "https://www.sciencebase.gov/catalog/item/6a0b0b0ab66b0188da36aedd"

## Initialize variables
wd <- "E:/MelinaStuff/BAM/WildTrax/WT-Integration"
organization <- "BAM"
dataset_code <- "BAM-BBS2024"
setwd(file.path(wd))

WTpj_Tbl <- read_sheet("https://docs.google.com/spreadsheets/d/1fqifS_E5O_IpW1B-UG_xthr9hzY6FIek-nFjCrt1G0w", sheet = "project")
lu <- "./lookupTables"
WT_spTbl <- read.csv(file.path(lu, "species_codes.csv"))
colnames(WT_spTbl) <- c("species_common_name", "species_code", "scientific_name")
WT_durMethTbl <- read.csv(file.path(lu, "duration_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distMethTbl <- read.csv(file.path(lu, "distance_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_durBandTbl <- read.csv(file.path(lu, "duration_interval_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distBandTbl <- read.csv(file.path(lu, "distance_band_codes.csv"), fileEncoding="UTF-8-BOM")

BBS_spTbl <- read.csv(file.path(lu, "BAM_BBS_speciesList.csv"))
XYtbl <- read.csv(file.path(lu, "BAM_BBS_stopsXY_july2022.csv"))

project_dir <- file.path("./project", dataset_code)
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}
out_dir <- file.path(wd, "out", dataset_code)    # where output dataframe will be exported
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

#Download from GoogleDrive
if (length(list.files(project_dir)) ==0) {
  pid <- WTpj_Tbl %>%
    filter(dataset_code =="BAM-BBS2024") %>%
    select("GSharedDrive location")
  gd.list <- drive_ls(as.character(pid))
  detection <- gd.list %>%
    filter(name =="50-StopData.zip") %>%
    select("id")
  drive_download(as_id(as.character(detection)), path = file.path(project_dir, "50-StopData.zip"))
  unzip(file.path(project_dir,"50-StopData.zip"), exdir= project_dir)
  folderName <- sub("(.*)\\..*$","\\1", "50-StopData")
  version <- list.files(file.path(project_dir,folderName))
  flist <- list.files(file.path(project_dir,folderName,version))
  lapply(flist,  function(x) {unzip(file.path(project_dir,folderName,version,x), exdir= project_dir)})
  
  visit <- gd.list %>%
    filter(name =="Weather.csv") %>%
    select("id")
  drive_download(as_id(as.character(visit)), path = file.path(project_dir, "Weather.csv"))
  #unzip(file.path(project_dir,"weather.zip.zip"), exdir= project_dir)
}

f_weather <- read.csv(file.path(project_dir,"Weather.csv"), header = TRUE)

#--------------------------------------------------------------
#
#       TRANSLATE
#
#--------------------------------------------------------------
# Temporal range
rangeYear <- c(2024)
############################
#### LOCATION TABLE ####
############################

location_tbl <- XYtbl %>%
  mutate(organization = organization,
         project = dataset_code,
         site = sapply(str_split(location, ":"), function(x) x[2]),
         station = sapply(str_split(location, ":"), function(x) x[3]),
         buffer_m = NA,
         location_visibility = "Visible",
         true_coordinates = TRUE,
         location_comments = NA,
         internal_wildtrax_id = NA)


############################
#### VISIT TABLE ####
############################
##Load source weather file and subset temporally using year range and geographically using location 
# subset only those with lat/long
fweather <-f_weather %>%
  select(-TempScale , -StartWind , -EndWind , -StartSky, -EndSky, -StartTemp, -EndTemp) %>%
  filter(CountryNum %in% location_tbl$CountryNum & StateNum %in% location_tbl$StateNum & Route %in% location_tbl$Route,
         Year %in% rangeYear,
         RPID < 200, #Drop experimental protocol (RPID >199 are experimental)
         QualityCurrentID ==1)  #Drop survey that do not meet BBS weather, date, time, and route completion criteria (QualityCurrentID = 0) 

visit_tbl <- fweather %>%
  right_join(location_tbl, by = c("CountryNum", "StateNum", "Route"), relationship = "many-to-many") %>%
  drop_na(Year) %>%
  mutate(visitDate = as.character(as.Date(with(location_tbl, paste(Year,Month,Day, sep = "-")),format = "%Y-%m-%d")),
         bait = "None",
         data_origin = "USGS download",
         survey_year = sub("\\-.*", "", visitDate),
         StartTime = ifelse(is.na(StartTime), 0, str_pad(StartTime, width = 4, pad = 0)),
         StartTime_BAM = times(format(strptime(StartTime, format="%H%M"), format = "%H:%M:%S")),
         EndTime = ifelse(is.na(EndTime), 0, str_pad(EndTime, width = 4, pad = 0)),
         EndTime_BAM = times(format(strptime(EndTime, format="%H%M"), format = "%H:%M:%S")),
         time_int = EndTime_BAM - StartTime_BAM,
         survey_time = as.character(StartTime_BAM + ((Stop - 1)*(time_int/50))),
         surveyDateTime = paste(visitDate, survey_time),
         observer = ObsN,
         comments = paste0("route time interval: ", StartTime, "-", EndTime),
         pkey_dt = paste(location, paste0(gsub("-", "", as.character(visitDate)),"_", gsub(":", "", survey_time)), observer, sep=":"),
)

############################
#### SURVEY TABLE ####
############################
##Protocol duration  (3 min)
survey_tbl<-visit_tbl %>%
  mutate(durationMethod = "0-3min",
         distanceMethod = "0m-400m")

#Species
fifties <- c("Fifty1.csv","Fifty2.csv","Fifty3.csv","Fifty4.csv","Fifty5.csv","Fifty6.csv","Fifty7.csv","Fifty8.csv","Fifty9.csv","Fifty10.csv")
outPtCount <- lapply(fifties, function(x){
  path2fifty = file.path(project_dir,"50-StopData",x)
  fifty <- read.csv(file = path2fifty, header = TRUE)
  # Subset temporal and geographical range at the route level using RouteDataID
  fifty <-subset(fifty, fifty$RouteDataID %in% survey_tbl$RouteDataID)

  # if fifty not empty -> transpose
  if(!nrow(fifty)==0){
    fiftygather <- gather(fifty, key=Stop, value = count, Stop1:Stop50)
    fiftygather$Stop <- gsub("Stop", "", fiftygather$Stop)
    
    #Attached pkey_dt
    fiftygather <-merge(fiftygather, survey_tbl[,c("RouteDataID", "Stop", "Year", "pkey_dt")], by=c("RouteDataID", "Stop", "Year"))
    # Calculate frequency table per stop
    #fifty_sumCOunt <- fiftygather %>% 
    #  group_by(pkey_dt) %>% 
    #  summarise(count = sum(count), .groups= "keep")
    
    fifty_sumCOunt <- setNames(aggregate(fiftygather$count, by=list(fiftygather$pkey_dt), FUN=sum), c("pkey_dt", "count"))
    fifty_Obs<- fiftygather[which(fiftygather$count!=0), ] 
    # Isolate REAL no observation (stops where no birds were seen)
    fifty_noObs <- fifty_sumCOunt[fifty_sumCOunt$count ==0,]
    if(nrow(fifty_noObs)>0){
      fifty_noObs$AOU <- 0
      fifty_all<-rbind.fill(fifty_Obs, fifty_noObs)
    }else{
      fifty_all<-fifty_Obs
    }

    # Add fifty REAL no observation (fifty_noObs)
    print(paste(x, "real obs:", as.character(nrow(fifty_all))))
    #---------
    ## Merge to location_tbl 
    #---------
    fifty_BAM <-merge(survey_tbl, fifty_all[,c("pkey_dt", "count", "AOU")], by=c("pkey_dt"))
    
    #Species
    fifty_BAM$species <- BBS_spTbl$WT_species_code[match(fifty_BAM$AOU, BBS_spTbl$AOU)]
    fifty_BAM$species[fifty_BAM$count ==0] <- "NONE"
    print(unique(fifty_BAM$species[!(fifty_BAM$species %in% WT_spTbl$species_code)]))
    
    fifty_BAM$scientificname <- WT_spTbl$scientific_name[match(fifty_BAM$species, WT_spTbl$species_code)]
    fifty_BAM$original_species <- fifty_BAM$AOU
    fifty_BAM$distanceband <- "0m-400m"
    fifty_BAM$durationinterval <- "0-3min"
    fifty_BAM$ind_count <- fifty_BAM$count
    
    #HEARD
    fifty_BAM$isHeard <- "DNC"
    #SEEN
    fifty_BAM$isSeen <- "DNC"
    #comments
    fifty_BAM$rawObserver <- fifty_BAM$ObsN
    fifty_BAM$raw_distance_code <- "0.25-mile radius"
    fifty_BAM$raw_duration_code <- "3-minute point count"
    fifty_BAM$originalBehaviourData <- NA
    fifty_BAM$missingindetections <- NA
    fifty_BAM$pc_vt <- "DNC"
    fifty_BAM$pc_vt_detail <- "DNC"
    fifty_BAM$age <- "DNC"
    fifty_BAM$fm <- "DNC"
    fifty_BAM$group <- "DNC"
    fifty_BAM$flyover <- "DNC"
    fifty_BAM$displaytype <- "DNC"
    fifty_BAM$nestevidence <- "DNC"
    fifty_BAM$behaviourother <- "DNC"
    fifty_BAM$utmZone <- "DNC"
    fifty_BAM$easting <- "DNC"
    fifty_BAM$northing <- "DNC"
    fifty_BAM$missinginlocations <- "DNC"
    fifty_BAM$time_zone <- "DNC"
    fifty_BAM$missinginvisit <- "DNC"
   
    #---------------------
    print(paste(x, "final:", as.character(nrow(fifty_BAM))))
    return(fifty_BAM)
    }
  })
data_flat <-do.call(rbind, outPtCount) # 2024


############################
##EXPORT
############################
dr<- drive_get(paste0("DataTransfered/",organization), shared_drive= "BAM_AvianData")

#Set GoogleDrive id
if (nrow(drive_ls(as_id(dr), pattern = dataset_code)) == 0){
  dr_dataset_code <-drive_mkdir(dataset_code, path = as_id(dr), overwrite = NA)
} else {
  dr_dataset_code <- drive_ls(as_id(dr), pattern = dataset_code)
}
dr_ls <- drive_ls(as_id(dr), pattern = dataset_code)

#--- Output columns
#WTlocation <- c("organization", "location", "latitude", "longitude", "buffer_m", "location_visibility", "true_coordinates", "location_comments", "internal_wildtrax_id")

WTsurvey <- c("location", "surveyDateTime", "durationMethod", "distanceMethod", "observer", "species", "distanceband",
              "durationinterval", "abundance", "isHeard", "isSeen", "comments")

Extended <- c("organization", "project","location", "surveyDateTime", "species", "distanceband", "durationinterval", "site", 
              "station", "utmZone", "easting", "northing", "missinginlocations", "time_zone", "data_origin", 
              "missinginvisit", "pkey_dt", "survey_time", "survey_year", "rawObserver", "original_species", 
              "scientificname", "raw_distance_code", "raw_duration_code", "originalBehaviourData", 
              "missingindetections", "pc_vt", "pc_vt_detail", "age", "fm", "group", "flyover", 
              "displaytype", "nestevidence", "behaviourother")

#location_tbl <- data_flat[!duplicated(data_flat[,WTlocation]), WTlocation] 
#write.csv(location_tbl, file= file.path(out_dir, "BAM-BBS2024_location.csv"), row.names = FALSE, na = "")
#location_out <- file.path(out_dir,"BAM-BBS2024_location.csv")
#drive_upload(media = location_out, path = as_id(dr_dataset_code), name = "BAM-BBS2024_location.csv", overwrite = TRUE) 

survey_tbl <- data_flat %>% 
  dplyr::filter(ind_count > 0) %>% 
  group_by(location, surveyDateTime, durationMethod, distanceMethod, observer, species, distanceband, durationinterval, isHeard, isSeen, comments) %>%
  dplyr::summarise(abundance = sum(ind_count), .groups= "keep")

write.csv(survey_tbl, file= file.path(out_dir, "BAM-BBS2024_survey.csv"), na = "", row.names = FALSE)
survey_out <- file.path(out_dir, "BAM-BBS2024_survey.csv")
drive_upload(media = survey_out, path = as_id(dr_ls), name = "BAM-BBS2024_survey.csv", overwrite = TRUE) 

extended_tbl <- data_flat[!duplicated(data_flat[,Extended]), Extended] 
write.csv(extended_tbl, file.path(out_dir, "BAM-BBS2024_behavior.csv"), na = "", row.names = FALSE)
beh_out <- file.path(out_dir, "BAM-BBS2024_behavior.csv")
drive_upload(media = beh_out, path = as_id(dr_ls), name = "BAM-BBS2024_behavior.csv", overwrite = TRUE) 


#---PROCESSING STATS
write_lines(paste0("Organization: ", organization), file.path(out_dir, paste0(dataset_code, "_stats.csv")))
write_lines(paste0("Project: ", dataset_code), file.path(out_dir, paste0(dataset_code, "_stats.csv")), append= TRUE)
nrow_location <- paste0("Number of locations: ", nrow(location_tbl))
write_lines(nrow_location, file.path(out_dir, paste0(dataset_code, "_stats.csv")), append= TRUE)
nrow_visit <- paste0("Number of visit: ", nrow(visit_tbl))
write_lines(nrow_visit, file.path(out_dir, paste0(dataset_code, "_stats.csv")), append= TRUE)
nrow_survey <- paste0("Number of survey: ", nrow(survey_tbl))
write_lines(nrow_survey, file.path(out_dir, paste0(dataset_code, "_stats.csv")), append= TRUE)

