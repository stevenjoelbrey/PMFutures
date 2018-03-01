# plot_BA_met_r_values.R 

# This script makes plots that display the r values of regional burn area
# correlation with meteorology values. These r values are calculated using
# R/plot_BA_Met_relationbs.R. These are values are pearson correlation coefs. 

library(ggplot2)
library(ggthemes)

# where the r (ceof) values live. 
dataDir <- "Data/correlations/"

df_6.2 <- get(load(paste0(dataDir, "ecoregion_6.2.RData")))
df_10.1 <- get(load(paste0(dataDir, "ecoregion_10.1.RData")))
df_11.1 <- get(load(paste0(dataDir, "ecoregion_11.1.RData")))
rm(r_df)

df <- rbind(df_6.2, df_10.1, df_11.1)
df$r <- as.numeric(df$r)
df$ecoregion <- factor(df$ecoregion, levels=c("6.2", "10.1", "11.1"))
df$ignitionType <- factor(df$ignitionType, levels=c("Lightning", "Human"))
df$variable <- factor(df$variable, levels=c("T2M","TP","RH"))

# Give nice names to ecoregions 
levels(df$ecoregion) <- c("Forested Mountains", "High Deserts", "Mediterranean California")

# Make the plot 
p <- ggplot(df, aes(x=variable, y=r, fill=ignitionType))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  scale_x_discrete("Meteorology variable",     
                   labels=c("T","Precip", "RH%"))+
  
  facet_wrap(~ecoregion)+ 
  scale_fill_manual(values=c("gray", "orange"))+
  scale_y_continuous("Pearson correlation")+
  guides(fill=guide_legend(title="Ignition"))+
  # geom_text(aes(label = paste(round(r, 2)), # for labeling the bars 
  #               vjust = ifelse(r >= 0, 0, 1)),
  #           fontface = "bold") +
  theme_tufte(ticks=F, base_size = 18)

# Save the plot
png(filename="Figures/regional_met_relations/correlation_summaruy.png",
    width=2300, height=1100, res=260)
p 
dev.off()