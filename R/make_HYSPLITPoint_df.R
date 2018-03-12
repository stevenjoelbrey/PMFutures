# make_HYSPLITPoint_df.R 


y1 <- 2007
y2 <- 2017

library(stringr)

# Read daily Hysplit point text file

# The header looks like the next line
# "Lon,    Lat, YearMmDd, HhMm,  Dur"
expectedHeaders <- c("Lon", "Lat", "YearMmDd", "HhMm", "Dur")

# Create an empty dataframe to rbind daily files to. 
merged_df <-  data.frame(matrix(ncol = length(expectedHeaders), nrow = 0))
names(merged_df) <- expectedHeaders

# Establish the location of the data 
dataDir <- "../HMS_ARCHIVE/"
years   <- y1:y2
nYears  <- length(years)

# Loop through each year of interest
for (y in 1:nYears){
  
  analysisYear <- years[y]
  print(paste("Working on year:", analysisYear))
  
  # Loop through each year in the HMS Archive
  yearDir <- paste0(dataDir, analysisYear, "/TEXT/")
  
  # get all files for year. List only text files that include "hmshysplit".
  # We do not care about HMS hotspots 
  yearFiles <- list.files(path=yearDir, pattern="hmshysplit")
  nFiles    <- length(yearFiles)
  
  # Handle individual daily files loop
  for (f in 1:nFiles){
  
    # read each file and check for wierdness
    dailyFile <- paste0(yearDir, yearFiles[f])
    df <- read.csv(dailyFile, stringsAsFactors = FALSE)
    
    # Check for desired headers
    if(any(names(df) != expectedHeaders)){
      stop(paste("One or more of the expected headers for:", 
                 dailyFile, 
                 "are missing."))
    } else{
      
      # Append the daily data. This will grow merged_df
      merged_df <- rbind(merged_df, df, stringsAsFactors=FALSE)
      
    }
  
  } # end of daily loop

} # end of year loop

################################################################################
# Begin quality control check on merged_df features. Remove all rows with 
# anything wierd. Create all masks such that TRUE is what we keep. 
################################################################################

# See if any YearMmDd are insane...
# aka, read date as R object, if it fails, exclude
yearString <- as.character(merged_df$YearMmDd)
DATES      <- as.POSIXct(yearString, format="%Y%m%d", tz="UTC")
dateMask   <- !is.na(DATES)

# Now place the R-formatted version of the dates into the dataframe
merged_df$DATE <- DATES

# Check for "****" in Dur and in "HhMm" columns. First white space must be removed
DUR <- str_replace(merged_df$Dur, " ", "") 
HhMm <- str_replace(merged_df$HhMm, " ", "") 
# quickly place version that has no spaces 
merged_df$Dur <- DUR
merged_df$HhMm<- HhMm

DURMask <- DUR != "****"
HhMmMask<- HhMm!= "****"

# It is easier to mask date range after subsetting the data by missing dates
# first
m1 <- dateMask & DURMask & HhMmMask
print(paste("Percent of rows retained after masking:", sum(m1)/length(m1)*100))
merged_df_subset <- merged_df[m1, ]

# Date outside realistic range
lowerMask <- merged_df_subset$DATE >= as.POSIXct(paste0(y1, "-01-01"), tz="UTC")
upperMask <- merged_df_subset$DATE <= as.POSIXct(paste0(y2, "-12-31"), tz="UTC")
m2 <- lowerMask & upperMask

hysplitPoints <- merged_df_subset[m2, ]

# Dur is meant to be HHMM but often leading 0 is left off. Where not length 4,
# append 0 up front 
DUR <- hysplitPoints$Dur
DurLength <- str_length(DUR)

DUR[DurLength==1] <- paste0("000", DUR[DurLength==1])
DUR[DurLength==2] <- paste0("00", DUR[DurLength==2])
DUR[DurLength==3] <- paste0("0", DUR[DurLength==3])

# Now get just the hourly quantity from the HHMM field
DurHours <- as.numeric(str_sub(DUR, 1, 2))
hysplitPoints$DurHours <- DurHours

# save away
save(hysplitPoints, file=paste0("Data/HMS/hysplitPoints_", y1,"_", y2,".RData"))




