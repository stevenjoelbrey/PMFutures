# plot_BA_met_r_values.R 

# This script makes plots that display the r (cor) values of regional burn area
# correlation with meteorology values. These r values are calculated using
# R/plot_BA_Met_relationbs.R. These values are pearson correlation coefs. 

# TODO: indicate which correlations are statistically different from zero. 

library(ggplot2)
library(ggthemes)

regionName <- "southeast_and_west"
lightningCol <- adjustcolor("cyan4", alpha.f = 1)
humanCol     <- adjustcolor("orange", alpha.f = 1)

# where the r (coef) values live. 
dataDir <- "Data/correlations/"

if(regionName == "west"){
  
  df_6.2 <- get(load(paste0(dataDir, "ecoregion_6.2.RData")))
  df_10.1 <- get(load(paste0(dataDir, "ecoregion_10.1.RData")))
  df_11.1 <- get(load(paste0(dataDir, "ecoregion_11.1.RData")))
  
  df <- rbind(df_6.2, df_10.1, df_11.1)
  df$r <- as.numeric(df$r)
  df$ecoregion <- factor(df$ecoregion, levels=c("6.2", "10.1", "11.1"))
  # Give nice names to ecoregions 
  levels(df$ecoregion) <- c("Forested Mountains", "High Deserts", "Mediterranean California")

}else if(regionName == "southeast"){
  
  df_8.3 <- get(load(paste0(dataDir, "ecoregion_8.3.RData")))
  df_8.4 <- get(load(paste0(dataDir, "ecoregion_8.4.RData")))
  df_8.5 <- get(load(paste0(dataDir, "ecoregion_8.5.RData")))
  df_15.4 <- get(load(paste0(dataDir, "ecoregion_15.4.RData")))
  
  df <- rbind(df_8.3, df_8.4, df_8.5, df_15.4)
  df$r <- as.numeric(df$r)
  df$ecoregion <- factor(df$ecoregion, levels=c("8.3", "8.4", "8.5", "15.4"))
  # Give nice names to ecoregions 
  levels(df$ecoregion) <- c("Southeast Plains", "Ozark", 
                            "Mississippi Alluvial", "Everglades")
} else {
  
  # both southeast and west together requested
  df_6.2 <- get(load(paste0(dataDir, "ecoregion_6.2.RData")))
  df_10.1 <- get(load(paste0(dataDir, "ecoregion_10.1.RData")))
  df_11.1 <- get(load(paste0(dataDir, "ecoregion_11.1.RData")))
  
  df_8.3 <- get(load(paste0(dataDir, "ecoregion_8.3.RData")))
  df_8.4 <- get(load(paste0(dataDir, "ecoregion_8.4.RData")))
  df_8.5 <- get(load(paste0(dataDir, "ecoregion_8.5.RData")))
  df_15.4 <- get(load(paste0(dataDir, "ecoregion_15.4.RData")))
  
  # combine each together
  df <- rbind(df_6.2, df_10.1, df_11.1, df_8.3, df_8.4, df_8.5, df_15.4)
  
  df$r <- as.numeric(df$r)
  df$ecoregion <- factor(df$ecoregion, levels=c("6.2", "10.1", "11.1", "8.3", "8.4", "8.5", "15.4"))
  # Give nice names to ecoregions 
  levels(df$ecoregion) <- c("Forested Mt.", "High Deserts", "Med. Cal.",
                            "S.E. Plains", "Ozark", "MS. Alluvial", "Everglades")
  
  
}

rm(r_df)

df$ignitionType <- factor(df$ignitionType, levels=c("Lightning", "Human"))
df$variable     <- factor(df$variable, levels=c("T2M","TP","RH"))

# Make the plot 
p <- ggplot(df, aes(x=variable, y=r, fill=ignitionType))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  scale_x_discrete("Meteorology variable",     
                   labels=c("T","P", "RH"))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  facet_grid(~ecoregion)+ 
  scale_fill_manual(values=c(lightningCol, humanCol))+
  scale_y_continuous("Pearson correlation")+
  guides(fill=guide_legend(title="Ignition"))+
  # geom_text(aes(label = paste(round(r, 2)), # for labeling the bars 
  #               vjust = ifelse(r >= 0, 0, 1)),
  #           fontface = "bold") +
  theme_tufte(ticks=T, base_size = 15)

# Save the plot
png(filename=paste0("Figures/regional_met_relations/correlation_summary_",
                    regionName,".png"),
    width=2400, height=1000, res=260)
print(p) 
dev.off()