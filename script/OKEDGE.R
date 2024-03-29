# ---
# title: "Translate Okanagan Edge Effects data"
# author: "Elly Knight"
# date: "July 14, 2022"
# Note on translation:
# --  Fix distance band > 100m as UNKNOWN to fit protocol
# --- Fix distance band > 250m as UNKNOWN to fit protocol

library(tidyverse)
library(readxl)
library(sf)
library(googledrive)

############################
## IMPORT ####
############################
## Initialize variables
wd <- "E:/MelinaStuff/BAM/WildTrax/WT-Integration"
setwd(wd)

organization = "Independent-Knight"
dataset_code = "OKEDGE"
dat <- "OK Database - PC data - Elly Knight.csv"
pc <- "OK Database - PC info - Elly Knight.csv"
plot <- "OK Database - plot info - Elly Knight.csv"
site <- "OK Database - site info - Elly Knight.csv"
lu <- "./lookupTables"
WT_spTbl <- read.csv(file.path("./lookupTables/species_codes.csv"))
colnames(WT_spTbl) <- c("species_common_name", "species_code", "scientific_name")
WT_durMethTbl <- read.csv(file.path("./lookupTables/duration_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distMethTbl <- read.csv(file.path("./lookupTables/distance_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_durBandTbl <- read.csv(file.path("./lookupTables/duration_interval_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distBandTbl <- read.csv(file.path("./lookupTables/distance_band_codes.csv"), fileEncoding="UTF-8-BOM")

#set working folder
project_dir <- file.path(wd, "project", dataset_code)
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}
data_source <- file.path(project_dir, dat)
if (!file.exists(data_source)) {
  #Download from GoogleDrive
  drive_download(paste0("sourceData/", dat), path = data_source)
}
plot_source <- file.path(project_dir, plot)
if (!file.exists(plot_source)) {
  #Download from GoogleDrive
  drive_download(paste0("sourceData/", plot), path = plot_source)
}
site_source <- file.path(project_dir, site)
if (!file.exists(site_source)) {
  #Download from GoogleDrive
  drive_download(paste0("sourceData/", site), path = site_source)
}
pc_source <- file.path(project_dir, pc)
if (!file.exists(pc_source)) {
  #Download from GoogleDrive
  drive_download(paste0("sourceData/", pc), path = pc_source)
}

out_dir <- file.path("./out", dataset_code)    # where output dataframe will be exported
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

site <- read.csv(site_source)
plot <- read.csv(plot_source)
pc <- read.csv(pc_source) 
dat <- read.csv(data_source)

############################
#### LOCATION TABLE ####
############################

#Read in template
location.temp <- read_xlsx("template/Template-WT-1Location.xlsx")

#Tidy raw data
location.ok <- left_join(plot, site)

#Reproject to lat long
location.nad83 <- location.ok %>% 
  st_as_sf(coords=c("UTM.E", "UTM.N"), crs=26911) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(longitude = X, latitude = Y)

#Format
location.all <- location.ok %>% 
  cbind(location.nad83) %>% 
  mutate(location = paste0("OKEDGE:", str_sub(Name, 1, 3), ":", Plot),
         elevationMeters = NA,
         bufferRadiusMeters = NA,
         isHidden = FALSE,
         trueCoordinates = TRUE,
         comments = NA,
         internal_wildtrax_id = NA,
         internal_update_ts = NA)

############################
#### VISIT TABLE ####
############################

#Read in template
visit.temp <- read_xlsx("template/Template-WT-2Visit.xlsx", sheet=1)[,-c(12:18)]

#Tidy raw data
visit.ok <- pc %>% 
  dplyr::filter(!is.na(Site),
                !is.na(Round))

#Format
visit.all <- visit.ok %>% 
  mutate(SiteChar = case_when(nchar(Site)==1 ~ paste0("00", Site),
                          nchar(Site)==2 ~ paste0("0", Site),
                          nchar(Site)==3 ~ paste0(Site)),
         location = paste0("OKEDGE:", SiteChar, ":", Plot),
         visitDate = Date, 
         snowDepthMeters = NA,
         waterDepthMeters = NA,
         crew = NA,
         bait = "NONE",
         accessMethod = NA,
         landFeatures = NA, 
         comments = Comments, 
         wildtrax_internal_update_ts = NA,
         wildtrax_internal_lv_id = NA)
  
############################
#### SURVEY TABLE ####
############################

#Read in template
survey.temp <- read_xlsx("template/Template-WT-3Surveys.xlsx", sheet=1)

#Check species codes
lu.sp <- read.csv("lookupTables/species_codes.csv")
survey.sp <- dplyr::filter(dat, !Species %in% lu.sp$species_code)
table(survey.sp$Species)

#Tidy raw data
survey.ok<- dat %>% 
  dplyr::filter(!is.na(Site)) %>% 
  dplyr::select(-ID) %>% 
  rename(ObsTime = Time) %>% 
  full_join(visit.ok %>% 
              dplyr::select(-ID)) %>%
  mutate(abundance = case_when(Species=="BHCOx50?" ~ 50,
                               !is.na(Species)~ 1, 
                               is.na(Species) ~ 0)) %>% 
  mutate(Observer = ifelse(Observer=="eck", "ECK", Observer),
         species= case_when(Species=="CEWA" ~ "CEDW",
                            Species=="BHCOx50?" ~ "BHCO",
                            Species=="TRSW" ~ "TRES",
                            Species=="raptor" ~ "URPT",
                            Species %in% c("SPSP", "spsp", "CCSP?", "GRSP?") ~ "UNSP",
                            Species %in% c("DOSP", "HUSP", "") ~ "UNBI",
                            Species=="GRPA" ~ "GRAP",
                            Species=="ECDO" ~ "EUCD",
                            !is.na(Species) ~ Species,
                            is.na(Species) ~ "NONE"))

#Check species codes again
survey.sp <- dplyr::filter(survey.ok, !species %in% lu.sp$species_code)
table(survey.sp$species)

#Get and create lookup tables
lu.obs <- survey.ok %>% 
  dplyr::select(Observer) %>% 
  unique() %>% 
  mutate(row = row_number(),
         observer = paste0("OKEDGE_0", row))

# Set >100m and >250m as UNKNOWN
lu.db <- survey.ok %>% 
  dplyr::select(Distance) %>% 
  unique() %>% 
  arrange(Distance) %>% 
  mutate(distanceband = c("UNKNOWN",
                          "0m-10m",
                          "500m-INF",
                          "UNKNOWN",
                          "500m-INF",
                          "UNKNOWN",
                          "500m-INF",
                          "500m-INF",
                          "0m-10m",
                          "10m-20m",
                          "100m-200m",
                          "20m-50m",
                          "200m-300m",
                          "200m-300m",
                          "300m-400m",
                          "400m-500m",
                          "50m-75m",
                          "500m-INF",
                          "75m-100m",
                          "UNKNOWN"))

#Format
survey.all <- survey.ok %>% 
  left_join(lu.obs) %>% 
  left_join(lu.db) %>% 
  mutate(SiteChar = case_when(nchar(Site)==1 ~ paste0("00", Site),
                              nchar(Site)==2 ~ paste0("0", Site),
                              nchar(Site)==3 ~ paste0(Site)),
         location = paste0("OKEDGE:", SiteChar, ":", Plot),
         Time = as.POSIXct(Time, format = "%Y-%m-%d %H:%M:%S"),
         Time = as.character(format(Time, format = "%H:%M:%S")),
         Time = case_when(is.na(Time) ~ "00:00:01",
                          !is.na(Time)~ Time),
         #Time = as.character(str_sub(Time, -8, -1)),
         surveyDateTime = as.character(paste(Date, Time)),
         durationMethod = "0-10min",
         distanceMethod = "0m-10m-20m-50m-75m-100m-200m-300m-400m-500m-INF",
         durationinterval = "0-10min",
         isHeard = ifelse(Activity %in% c("Singing", "Calling", "chipping", "Defend", "Drumming", "Chick beg call", "singing", "Display", "drumming"), "Yes", "No"),
         isSeen = ifelse(!'Visual.' %in% c("", "`", NA), "Yes", "No")) %>% 
  rename(comments = Comments)


############################
#### EXTENDED TABLE ####
############################

#Read in template
extended.temp <- read_xlsx("template/Template-WT-4Extended.xlsx", sheet=1)

#Get and create lookup tables
lu.vt <- survey.all %>% 
  dplyr::select(Activity) %>% 
  unique() %>% 
  mutate(pc_vt = c("Song", "Call", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "Call", "Call", "Non-vocal", "Call", "Song", "Song", "Non-vocal", NA),
         pc_vt_detail=c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "Drumming", NA, NA, "Display", "Drumming", NA))

lu.agefm <- survey.all %>% 
  dplyr::select(Sex) %>% 
  unique() %>% 
  mutate(age = c("Adult", NA, "Adult", "Adult", "Juvenile", "Fledgling", "Adult", NA, "Adult", NA),
         fm = c("Male", NA, "Female", "Female", "Unknown", "Unknown", "Unknown", "NA", "Unknown", NA)) 

#Format
extended.all <- survey.all %>% 
  left_join(location.all %>% 
              dplyr::select(-Distance, -comments)) %>% 
  left_join(lu.sp %>% 
              rename(species = species_code)) %>% 
  left_join(lu.vt) %>% 
  left_join(lu.agefm) %>% 
  rename(site = SiteChar,
         station = Plot,
         utmZone = 11,
         easting = UTM.E,
         northing = UTM.N,
         rawObserver = Observer,
         original_species = Species,
         scientificname = scientific_name,
         rawDistanceCode = Distance,
         originalBehaviourData = Activity) %>% 
  mutate(missinginlocations = NA,
         time_zone = "MDT",
         data_origin = NA,
         missinginvisit = NA,
         pkey_dt = paste0(location,":", surveyDateTime, ":", observer),
         survey_year = as.numeric(str_sub(Date, 1, 4)),
         survey_time = str_sub(Time, 1, 5),
         rawDurationCode = NA,
         missingindetections = NA,
         group=NA,
         flyover = ifelse(originalBehaviourData=="Fly over", "Yes", "No"),
         displaytype = NA,
         nestevidence=NA,
         behaviourother = NA,
         atlas_breeding_code = NA)
         
############################
## EXPORT ####
############################
#Extract GoogleDrive id to store output
dr<- drive_get(paste0("toUpload/",organization))

if (nrow(drive_ls(as_id(dr), pattern = dataset_code)) == 0){
  dr_dataset_code <-drive_mkdir(dataset_code, path = as_id(dr), overwrite = NA)
} else {
  dr_dataset_code <- drive_ls(as_id(dr), pattern = dataset_code)
}
#Location
WTlocation <- c("location", "latitude", "longitude")

# Remove duplicated location
location_tbl <- location.all[!duplicated(location.all[,c("location")]), WTlocation] # 
write.csv(location_tbl, file= file.path(out_dir, paste0(dataset_code,"_location.csv")), row.names = FALSE, na = "")
location_out <- file.path(out_dir, paste0(dataset_code,"_location.csv"))
drive_upload(media = location_out, path = as_id(dr_dataset_code), name = paste0(dataset_code,"_location.csv"), overwrite = TRUE) 


#Visit
WTvisit <- c("location", "visitDate", "snowDepthMeters", "waterDepthMeters", "crew", "bait", "accessMethod", "landFeatures", "comments", 
             "wildtrax_internal_update_ts", "wildtrax_internal_lv_id")

#Delete duplicated based on WildtTrax attributes (double observer on the same site, same day). 
visit_tbl <- visit.all[!duplicated(visit.all[,WTvisit]), WTvisit] 
write.csv(visit_tbl, file= file.path(out_dir, paste0(dataset_code,"_visit.csv")), quote = FALSE, row.names = FALSE, na = "")
visit_out <- file.path(out_dir, paste0(dataset_code,"_visit.csv"))
drive_upload(media = visit_out, path = as_id(dr_dataset_code), name = paste0(dataset_code,"_visit.csv"), overwrite = TRUE) 



#Survey
survey_tbl <- survey.all %>% 
  group_by(location, surveyDateTime, durationMethod, distanceMethod, observer, species, distanceband,
           durationinterval, isHeard, isSeen, comments) %>% 
  dplyr::summarise(abundance = sum(abundance), .groups= "keep")
WTsurvey <- c("location", "surveyDateTime", "durationMethod", "distanceMethod", "observer", "species", "distanceband",
              "durationinterval", "abundance", "isHeard", "isSeen")
survey_tbl <- survey_tbl[!duplicated(survey_tbl[,WTsurvey]), WTsurvey] 


write.csv(survey_tbl, file= file.path(out_dir, paste0(dataset_code,"_survey.csv")), quote = FALSE, row.names = FALSE, na = "")
survey_out <- file.path(out_dir, paste0(dataset_code,"_survey.csv"))
drive_upload(media = survey_out, path = as_id(dr_dataset_code), name = paste0(dataset_code,"_survey.csv")) 





#Extend
extended.wt <- extended.all %>% 
  dplyr::select(colnames(extended.temp))
write.csv(extended.wt, "out/OKEDGE_extended.csv", row.names = FALSE)
