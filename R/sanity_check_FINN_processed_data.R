# sanity_check_FINN_processed_data.R

# ------------------------- Description ----------------------------------------
# This script is used to perform basic checks on the FPA FOD wildfires that 
# were put through the FINN framework to have PM25 emissions estimated. It will
# save out the processed data with FPA FOD attributes appended. 

# ------------------------- Arguments ------------------------------------------
processData <- TRUE
mergeYears  <- TRUE

years <- 2002:2015

library(rgdal)
library(maps)
library(ggplot2)

if(processData){
  for (year in years){

  # Read in the exact data that was given to be processed
  data_projected <- readOGR(dsn="Data/FINN/shape/", layer=paste0("FPA_FOD_", year))
  df_in <- data_projected@data
  
  # Now load the data that was subset to be processed by FINN (all columns)
  load(paste0("Data/FINN/FPA_FOD_forFINN_", year, ".RData"))
  
  # See if the unique IDs match for these dataframes
  if( !(sum(df_in$FOD_ID == df$FOD_ID)/length(df$FOD_ID)==1) ){
    stop("The data do not have the same number of rows, attributes lost")
  } else{
    # If the code is still running it means all the rows matched. That means we
    # can use df instead of df_in, since they match were they have the same
    # data and df has many more attributes (columns).
    rm(df_in, data_projected)
  }
  
  # Read in the data processed by Christine 
  df_processed <- read.csv(paste0("Data/FINN/FINNv1.5_BREY_",year,"_061220182.txt"))
  print(names(df_processed))
  
  # Do the dimensions (rows) match the file you gave? Different columns by design
  print(paste("Original dimensions", dim(df)))
  print(paste("Processed dimensions", dim(df_processed)))
  if(dim(df)[1] != dim(df_processed)[1]){
    stop("Fires were lost or added in the FINN processing.")
  }else{
    print("Same number of wildfires before and after processing.")
  }
  
  # Match the rows of these data VERY carefully. match(): "match returns a vector 
  # of the positions of (first) matches of its first argument in its second."
  ProcessedLocationIndf<- match(df_processed$FOD_ID, df$FOD_ID)
  
  # Since match returns the first instance of a match, make sure every match
  # indicie is unique, otherwise it could mean something strange is happenning. 
  if(length(unique(ProcessedLocationIndf)) != dim(df)[1]){
    stop("There is an issue with the data being duplicated or repeated")
  } else{
    print("All FOD_ID in processed data have exactly one match in df")
  }
  
  # Make sure the match aligns
  if(unique(df_processed$FOD_ID == df$FOD_ID[ProcessedLocationIndf])){
    print("All data aligned")
  }else{
    stop("The data are not alignged. Something is backwards.")
  }
  
  # Combine the data by re-ordering the original data to match the processed data
  # We know the seperate FOD_IDs match, rename the second one, so that functions
  # still work. 
  names(df_processed)[names(df_processed) == "FOD_ID"] <- "FOD_ID_copy"
  df_combined <- cbind(df[ProcessedLocationIndf,], df_processed)
  print(names(df_combined))
  
  # Make sure attributes that were a part of both match
  dLat <- df_combined$lat - df_combined$LATITUDE
  dLon <- df_combined$long - df_combined$LONGITUDE
  
  if( max(abs(dLat)) > 1e-5 | max(abs(dLon)) > 1e-5){
    stop("A significant different in latitude or longitude was found.")
  }
  
  png(filename=paste0("Data/FINN/Figures/locations_",year,".png"), 
      width=3000, height=1900, res=250)
  plot(df_combined$LONGITUDE, df_combined$LATITUDE, pch=".", col="orange",
       xlab="Longitude", ylab="Latitude")
  points(df_combined$long, df_combined$lat, pch=".", col=df_combined$genLC)
  map("state", add=T)
  title("Orange dots should be hidden by range of genLC colors")
  legend("bottomleft",
         legend=unique(df_combined$genLC),
         fill=unique(df_combined$genLC),
         title="gen LC #")
  dev.off()
  
  # Map and color code the points by land cover
  png(filename=paste0("Data/FINN/Figures/genLC_dist_",year,".png"))
  hist(df_combined$genLC, main="distribution of genLC")
  dev.off()
  
  # Plot emissions [PM2.5] vs. fire size for each land cover type
  png(filename=paste0("Data/FINN/Figures/PM_v_size_",year,".png"))
  plot(df_combined$PM25, df_combined$AREA_m2, col=df_combined$genLC, pch=19)
  title("PM25 emitted vs. fire size colored by genLC")
  dev.off()
  
  # Look at human and lightning emission totals
  lightningMask <- df_combined$STAT_CAUSE_DESCR == "Lightning"
  
  # Total emissions of each. 
  lightningPM25 <- sum(df_combined$PM25[lightningMask])
  humanPM25     <- sum(df_combined$PM25[!lightningMask])
  
  print(lightningPM25)
  print(humanPM25)
  
  png(filename=paste0("Data/FINN/Figures/PM_total_",year,".png"))
  barplot(height=c(lightningPM25, humanPM25), col=c("blue", "orange"),
          ylab="kg")
  text(x=c(0.5,1.5), y=c(-100000, -100000), labels=c(lightningPM25, humanPM25))
  title(paste("CONUS wildfire PM25 total emissions for", year))
  dev.off()
  
  # Save the file 
  FPA_FOD <- df_combined
  rm(df_combined)
  save(FPA_FOD, file=paste0("Data/FINN/processed/FPA_FOD_wFINN-PM25_",year,".RData"))
  
  print(paste0("Finished processing ", year, " FPA FOD w/FINN fires."))

}
}

# Merge the yearly processed data into a single file
if(mergeYears){
  nRow <- 0
  for (year in years){
    
    if(year == years[1]){
      df_base <- get(load(paste0("Data/FINN/processed/FPA_FOD_wFINN-PM25_",year,".RData")))
      rm(FPA_FOD)
      nRow <- nRow + dim(df_base)[1]
    } else{
      # Loads FPA_FOD object
      load(paste0("Data/FINN/processed/FPA_FOD_wFINN-PM25_",year,".RData"))
      nRow <- nRow + dim(FPA_FOD)[1]
      df_base <- rbind(df_base, FPA_FOD)
      rm(FPA_FOD)
    }
      
  }
  if(nRow != dim(df_base)[1]){
    stop("Data was lost in the merge")
  }
  # Rename final merged data
  FPA_FOD <- df_base
  # Get rid of duplicate column with name FPA_FOD
  
  save(FPA_FOD, file=paste0("Data/FINN/processed/FPA_FOD_wFINN-PM25_2002-2015.RData"))
}