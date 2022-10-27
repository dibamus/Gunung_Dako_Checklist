#Gunung Dako Analysis
setwd("path to your working directory with the scripts and data")

# Load in these libraries
library("tidyverse")
library("readxl")
library("vegan")
library("elevatr")
library("reshape2")
library("leaflet")
library("cowplot")
library("sf")

#source in the functions necessary for analysis
source("checklistFunctions2.R")

##SETUP - color palette for plotting
colors <- setNames(c("#B5BA72","#a8ccde","#3b5374","#586A6A"),
                   nm = c("frog","lizard","snake","cumulative"))

#### STEP 1 - Data Import/Cleanup ####
df <- read_excel("Supplimentary Table S1 -DakoCatalog.xlsx",sheet = 1) %>% mutate(binom = "Genus_")

colnames(df) <- sub(" ","_",colnames(df))
colnames(df) <- gsub("[()]","",colnames(df))

df <- df %>% mutate(binom = Genus_species) 

#Remove specimens with uncertain affinities
df <- df %>%
  filter(Genus_species != "Sphenomorphus \"nigrilabris\" indet.") %>% #exclude indeterminate Sphenomorphus specimens
  filter(JAM_Number != 16552) %>% #exlude unwanted specimen in city
  filter(JAM_Number != 16535) #exclude turtle (Leucocephalon yunwoi)

#### STEP 2 - Add groups/Look up Elevation####

cutoffs <- c(400,850,1500) #- cutoff values for elevational bands

df <- df %>% addGroups() %>% elevBands(cutoffs)

#put specimens missing coords in the middle "upland" elevational band (found near 1000m camp)
df$eband[which(df$JAM_Number %in% c(15998,16201,16409))] <- "850-1500m" 

#### STEP 3 - Generate Plots ####
accumulationPlots <- accCurve(df,colors)
elevationPlot <- plotElev(df,colors)

alignedAccumulation <- align_plots(plotlist = list(accumulationPlots$taxonacc + theme(legend.position = "none"),
                                                   accumulationPlots$bandacc + theme(legend.position = "none")),
                                   align = "hv",
                                   axis = "blr")

accumulation <- plot_grid(alignedAccumulation[[1]], alignedAccumulation[[2]], 
                          rel_widths = c(3,2),
                          rel_heights = c(2,1))

#### STEP 4 - Generate Tables ####

#Table 1
t1 <- speciesTable(filter(df))

#Table 2
t2 <- missingSpecies(df) 

#### STEP 5 - Calculate beta diversity for elevational bands #####

beta <- betadiver(veganize(df, elevation = TRUE))

#### STEP 6 - Leaflet Interactive Map #####
elevationcolors <- colorFactor(c("#A8CCDE","#4E819A","#3B5374"), df$eband)

dakoPoly <- st_read("shp_0/WDPA_WDOECM_Sep2022_Public_555571285_shp-polygons.shp")
# 

map <- df %>% 
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldTopo) %>% # using ESRI World Topo for the background map tiles
  addPolygons(data = dakoPoly$geometry,
              color = "#586A6A", weight = 2,
              fillColor = "#B5BA72",
              smoothFactor = 0.2,
              opacity = 1.0, fillOpacity = 0.5) %>%
  
  addCircleMarkers(label = paste(df$JAM_Number,"-",df$Genus_species),
                   color = ~elevationcolors(df$eband),
                   popup = paste(df$JAM_Number,"-",df$Genus_species,"<br>",
                                 df$Collection_Date,"<br>",
                                 "Lat =",df$Latitude,"  ","Lon =",df$Longitude,"  ","Elev =",df$Elevation,"<br>",
                                 df$Habitat),
                   radius=4) %>% #this adds in markers. on mouse-over, it will say their JAM ID and the binomial
  
  setView(lat = df$Latitude[10], lng = df$Longitude[10],zoom = 11) %>%
  addScaleBar(position = "bottomleft",
              options= scaleBarOptions(metric = TRUE))

#### Push files to gDrive ####
ggsave(filename = paste0(mtn, "_AccumulationPlot.png"),
       plot = accumulation,
       width = 8, height = 6, units = "in",
       bg = "white")

ggsave(filename = paste0(mtn, "_ElevationPlot.png"),
       plot = elevationPlot[[1]],
       width = 8, height = 6.2, units = "in",
       bg = "white")
source('gDrive.R') #not included in published code
