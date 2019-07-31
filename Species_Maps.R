# This is the code to produce a species map

# Libraries ---------------------------------------------------------------

library(ggplot2)
library(rgdal) # reads in shape files (requires GDAL)
library(broom) # to manage shapefiles


# Data --------------------------------------------------------------------

# Loading data:
data <- read.csv(file = "ESM_5.csv", header = T)
# Note: these maps use the breeding probabilities calculated using the data from each site from the years in which each site was visited. See the Compiling_Results script for some notes on creating different variables from the data (such as probability of presence) which could be mapped instead.



# read the shapefile (currently using GADM: https://gadm.org/ map of Antarctica, which has been cut to only include the Antarctic peninsula):
peninsula <- readOGR(dsn = "./Shapefiles", layer = "GADM_peninsula")
# Turn the shapefile into an object that ggplot2 can use:
peninsula <- fortify(peninsula)


# Mapping -----------------------------------------------------------------

# Make Basic Map
basic <- ggplot() +
  geom_polygon(data = peninsula, aes(x = long, y = lat, group = group),
               size = 0.25, fill = "gray40", alpha = 1) +
  
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_line(colour = "white", linetype = "dotted")) +
  coord_map(project="stereographic", orientation=c(-90,-55,0),
            ylim = c(-68, -60), xlim = c(-70, -42)) +
  labs(x = "", y = "")


# Function that adds data of a particular species to the basic map:
SpMap <- function(base, dat, species, let) {
  frame <- dat[,c(1:6,which(colnames(dat)==species))] # subsets the data
  colnames(frame)[7] <- "Probability" # standardizes the name of the data column
  frame <- frame[order(frame$Probability),] # Orders the sites
  
  base +
    geom_point(data = frame,
               aes(x = Lon, y = Lat, colour = Probability),
               size = 4,
               shape = 19, alpha = 0.8) +
    theme(plot.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12)) +
    scale_color_gradient(low = "#FFFFFF", high = "#990056") +
    ggtitle(paste("Breeding Probability: ", species, sep = ""))
}


# Then this 
SpMap(base = basic, dat = data, species = "CAPE")







