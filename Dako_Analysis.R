#Gunung Dako Analysis
#setwd("path to your working directory with the scripts and data")

# Load in these libraries
require("tidyverse")
require("readxl")
require("vegan")
require("elevatr")
require("reshape2")
require("leaflet")
require("cowplot")
require("sf")
require('ggplot2')
require('cowplot')
require('divDyn')

#source in the functions necessary for analysis
source("checklistFunctions2.R")

##SETUP - color palette for plotting
colors <- setNames(c("#B5BA72","#a8ccde","#3b5374","#586A6A"),
                   nm = c("frog","lizard","snake","cumulative"))

#### STEP 1 - Data Import/Cleanup ####
df <- read_excel("Supplimentary Table S1 -DakoCatalog.xlsx",sheet = 1)

colnames(df) <- sub(" ","_",colnames(df))
colnames(df) <- gsub("[()]","",colnames(df))

df <- df %>% mutate(binom = Genus_species) 

#Remove specimens with uncertain affinities
df <- df %>%
  filter(JAM_Number != 16552) %>% #exlude unwanted specimen in city
  filter(JAM_Number != 16535) %>% #exclude turtle (Leucocephalon yunwoi)
  filter(JAM_Number != 16045) #exclude missing specimen

#all the S. "nigrilabris indeterminate" are above 1000m, so put them in the "high" species
df$binom[which(df$binom == "Sphenomorphus \"nigrilabris\" indet.")] <- "Sphenomorphus nigrilabris \"high elevation\""

#lump the indet. O. semipalmata in with the high ones
df$binom[which(df$binom == "Occidozyga semipalmata indet.")] <- "Occidozyga semipalmata \"high elevation\""

#C. "aspinosus" to C.sp. "aspinosus"
df$binom[which(df$binom == "Occidozyga semipalmata")] <- "Occidozyga semipalmata \"low elevation\""

#### STEP 2 - Add groups/Look up Elevation####

#### OPTIONAL: Check elevatr accuracy ####
elevs <- filter(df, !is.na(Latitude))[,c("Longitude","Latitude")] %>% 
  as.data.frame()%>%
  get_elev_point(prj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                 src = 'aws',
                 z = 14)
df$elevatr <- NA

df$elevatr[which(!is.na(df$Latitude))] <- elevs@data$elevation

df$elev_diff <- df$Elevation-df$elevatr


ggplot(df)+
  geom_density(aes(x = elev_diff))

ggplot(df)+
  geom_point(aes(x = Elevation, y= elevatr))

elev_corr <- lm(elevatr ~Elevation, df, na.action = na.omit)

#elevatr is biased to +2.8m vs our gps points
correctionfactor <- mean(df$elev_diff%>% na.exclude())

checkelevs <- filter(df, elev_diff < -50 | elev_diff >50) #r-squared of 99.989

#### continue ####
cutoffs <- c(700,1400) #- cutoff values for elevational bands

df <- df %>% addGroups() %>% elevBands(cutoffs,cf=2.8)

#put specimens missing coords in the 700-1400 elevational band (found near 1000m camp)
df$eband[which(df$JAM_Number %in% c(15998,16201,16409))] <- levels(df$eband)[2]

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
plot_grid(elevationPlot$rtDiv + 
            coord_cartesian() + 
            ylab("Species Richness") + 
            theme(axis.text.y = element_text(color = "black")),
       elevationPlot$scatter + 
         coord_flip() + 
         scale_x_discrete(limits = rev) + 
         theme(axis.text.x = element_text(face = "plain", angle = 0,vjust = 0.5, hjust=0.5),
               axis.text.y = element_text (face = "italic")) +
         facet_grid(rows = vars(group), scales = "free", space = "free"),
        nrow = 2, rel_heights = c(1,4), align = "v", axis = "rl")


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

df$Habitat[which(is.na(df$Habitat))] <- "not recorded"

map <- df %>% 
  leaflet(options = leafletOptions(
    attributionControl=FALSE,zoomControl = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldTopo) %>% # using ESRI World Topo for the background map tiles
  addPolygons(data = dakoPoly$geometry,
              color = "#586A6A", weight = 2,
              fillColor = "#B5BA72",
              smoothFactor = 0.2,
              opacity = 1.0, fillOpacity = 0.5) %>%
  
  addCircleMarkers(label = paste(df$JAM_Number,"-",df$binom),
                   color = ~elevationcolors(df$eband),
                   popup = paste("<center>",df$JAM_Number, "-", df$binom, "<br>",
                                 "Collected ",df$Collection_Date, "</center>","<br>",
                                 "Elevation: ", df$Elevation, "m", "<br>",
                                 "Habitat: ", df$Habitat),
                   radius=4) %>% #this adds in markers. on mouse-over, it will say their JAM ID and the binomial
  
  setView(lat = df$Latitude[10], lng = df$Longitude[10],zoom = 11) %>%
  addScaleBar(position = "bottomleft",
              options= scaleBarOptions(metric = TRUE))


#### Push files to gDrive ####

mtn <- "Dako"
source('gDrive.R') #not included in published code; pushes results to shared folder
