# Code to produce a species accumulation curve for any site


# Libraries ---------------------------------------------------------------

library(vegan) # for rarefaction techniques
library(HDInterval) # for the hdi function



# Data --------------------------------------------------------------------

# Loading the species accumulation results:
load(file = "accum.results.Rdata")

# Loading the edited model results results:
load(file = "breed.results.edited.Rdata")

# Loading Effort Data:
load(file = "EffortData.RData")

# Loading the raw visit data:
visitdata <- read.csv("Raw_ASI_Data.csv", header = TRUE)



# Preparing Model Data ----------------------------------------------------

# This summarizes the species richness by site x year x iteration
sa.array <- apply(accum.results, MARGIN = c(1,2,4), FUN = sum)




# Preparing Raw Data ------------------------------------------------------

# We first need to summarize the raw data so that we can view raw species richness in the same way as the model output.

# Subset the raw data to create a dataframe to manipulate the raw data:
raw.data.subset <- visitdata[,c("Code", "Season",
                               "ADPE",
                               "ANSH",
                               "ANTE",
                               "AnySkua",
                               "BBSP",
                               "BRSK",
                               "CAPE",
                               "CHPE",
                               "GEPE",
                               "KEGU",
                               "MACP",
                               "SGPE",
                               "SNPE",
                               "SNSH",
                               "SOFU",
                               "SPSK",
                               "WISP")]

# Determining breeding:
for (i in 1:nrow(raw.data.subset)) {
  for (j in 3:ncol(raw.data.subset)) {
    if (raw.data.subset[i,j] != 2 | is.na(raw.data.subset[i,j])) {
      raw.data.subset[i,j] <- 0
    } else {raw.data.subset[i,j] <- 1}
  }
}
# Reformatting season:
for (i in 1:nrow(raw.data.subset)) {
  raw.data.subset$NumSeason[i] <- as.numeric(
    substring(raw.data.subset$Season[i],first = 1, last = 4)
  )
}
# Create Total Skua category for raw data:
for (i in 1:nrow(raw.data.subset)) {
  temp <- sum(raw.data.subset[i,c("BRSK","SPSK")])
  raw.data.subset$TotSkua[i] <- max(temp, raw.data.subset$AnySkua[i])
}
# Re-subsetting the data:
raw.data.subset <- raw.data.subset[,c("Code","NumSeason",
                                      "ADPE",
                                      "ANSH",
                                      "ANTE",
                                      "BBSP",
                                      "CAPE",
                                      "CHPE",
                                      "GEPE",
                                      "KEGU",
                                      "MACP",
                                      "SGPE",
                                      "SNPE",
                                      "SNSH",
                                      "SOFU",
                                      "WISP",
                                      "TotSkua")]
colnames(raw.data.subset)[2] <- "Season" # Renaming the Season column
# Number of years visited from the most-visted site:
max.visits.yr <- max(visit.tot.yr.tab)
# List of sites:
site.list <- dimnames(visit.tot.yr.tab)[[1]]
# Create an empty array for site x species x year
raw.array.yr <- array(NA, dim = c(nrow(visit.tot.yr.tab),
                                  ncol(raw.data.subset)-2,
                                  max.visits.yr),
                      dimnames = list(site.list,
                                      colnames(raw.data.subset)[-c(1:2)],
                                      min(raw.data.subset$Season):
                                        max(raw.data.subset$Season)))
# Fill that array:
for (i in 1:dim(raw.array.yr)[1]) {
  temp <- raw.data.subset[which(raw.data.subset$Code == site.list[i]),]
  temp.table <- table(temp$Season)
  for (k in 1:length(temp.table)) { # year
    for (j in 3:ncol(temp)) { # species
      raw.array.yr[i,j-2,names(temp.table[k])] <- 
        max(temp[which(temp$Season <= names(temp.table[k])),j])
    }
  }
}
# Raw richness by year (NA for years without visits):
raw.sa.yr.mat <- apply(raw.array.yr, MARGIN = c(1,3), FUN = sum, na.rm = F)




# Projecting diversity ----------------------------------------------------

# This function projects species accumulation using a variety of indices at a site (see the vegan package and the function specpool() for more details. The skuas are again split into two columns, but in this case represent "at least one species of skua", and "a second species of skua". It would be impossible to get the true pattern of each species, since so many of the original records are unidentified to species.
Specpool <- function(data, site) {
  dat <- data[site,,]
  dat <- dat[rowSums(is.na(dat)) != ncol(dat), ]
  Skua1 <- as.numeric(dat[,"TotSkua"] >= 1)
  Skua2 <- as.numeric(dat[,"TotSkua"] == 2)
  dat <- cbind(dat, Skua1, Skua2)
  dat <- dat[,which(colnames(dat) != "TotSkua")]
  res <- specpool(dat)
  return(res)
}

# For the function to work with the raw data, they need to be reorganized:
# Create an empty array for raw projections (site x year x species)
proj_dat_raw <- array(NA, dim = c(length(site.list),
                                  max.visits.yr,
                                  ncol(raw.data.subset)-2),
                      dimnames = list(site.list,
                                      min(raw.data.subset$Season):
                                        max(raw.data.subset$Season),
                                      colnames(raw.data.subset)[-c(1:2)]))
# Fill that array:
for (i in 1:dim(proj_dat_raw)[1]) {
  temp <- raw.data.subset[which(raw.data.subset$Code == site.list[i]),]
  temp.table <- table(temp$Season)
  for (j in 1:length(temp.table)) { # year
    for (k in 3:ncol(temp)) { # species
      proj_dat_raw[i,names(temp.table[j]),k-2] <- 
        max(temp[which(temp$Season == names(temp.table[j])),k])
    }
  }
}
rm(temp.table, i, j, k) # cleaning up


# List of sites with more than one visit (the projections will only work for those sites):
specpool_list <- as.character(site.list[which(visit.tot.yr.tab[,"2016-17"] > 1)])


# This object will be a framework for both raw and model projection results:
specpool_res <- as.data.frame(site.list)
names(specpool_res) <- "Code"
specpool_res$Years_Visited <- visit.tot.yr.tab[,"2016-17"]

# Projections for raw data:
temp.frm <- data.frame()
for (i in 1:length(specpool_list)) {
  temp <- Specpool(data = proj_dat_raw, site = specpool_list[i])
  temp$Code <- specpool_list[i]
  temp.frm <- rbind(temp.frm, temp)
}
specpool_res_raw <- merge(specpool_res, temp.frm, by = "Code", all.x = TRUE)

# Projections for model data (separately for each iteration)
# This step takes a long time:
specpool_res_mod <- array(NA, dim = c(length(site.list),
                                      dim(specpool_res_raw)[[2]]-1,
                                      dim(breed.results.edit)[[4]]),
                          dimnames = list(site.list,
                                          dimnames(specpool_res_raw)[[2]][-1],
                                          NULL))
for (j in 1:dim(specpool_res_mod)[[3]]) {
  temp.frm <- data.frame()
  for (i in 1:length(specpool_list)) {
    temp <- Specpool(data = breed.results.edit[,,,j], site = specpool_list[i])
    temp$Code <- specpool_list[i]
    temp.frm <- rbind(temp.frm, temp)
  }
  specpool_res_mod_tmp <- merge(specpool_res, temp.frm, by = "Code", all.x = TRUE)
  specpool_res_mod[,,j] <- as.matrix(specpool_res_mod_tmp[,-1])
}

rm(specpool_res_mod_tmp, temp.frm, temp, i, j) # Cleaning up



# Accumulation Plots ------------------------------------------------------

# The function for making the plot:
PlotAccumYrProj <- function(site) {
  dat <- sa.array[site,,] # extract accumulation data
  proj_raw <- specpool_res_raw[which(specpool_res_raw$Code == site),]
  proj_mod <- specpool_res_mod[site,,]
  vis <- visit.tot.yr.tab[site,] # extract number of years visited over time
  raw <- as.numeric(raw.sa.yr.mat[site,])
  yrs <- which(visit.yr.tab[site,] == 1)
  
  if (vis["2016-17"] <= 10) {
    lim <- 10} else {lim <- 22}
  
  # Extracting model accumulation data:
  datf <- as.data.frame(apply(dat, MARGIN = 1, FUN = median))
  colnames(datf) <- "median"
  datf$upper95 <- apply(dat, MARGIN = 1, FUN = hdi, credMass = 0.95)[2,]
  datf$lower95 <- apply(dat, MARGIN = 1, FUN = hdi, credMass = 0.95)[1,]
  datf$vis <- as.numeric(vis)
  
  # Basic Plot:
  plot(NULL,
       ylim = c(0,16),
       xlim = c(1,lim),
       xlab = NA, ylab = NA)
  title(main = site, ylab = "Richness", xlab = "Cumulative Years Visited")
  
  # Model accumulation data:
  polygon(x = c(datf$vis[yrs], rev(datf$vis[yrs])),
          y = c(datf$lower95[yrs], rev(datf$upper95[yrs])),
          col = rgb(0.7,0.7,0.7,0.5),
          border = F)
  lines(median ~ vis, data = datf[yrs,], lwd = 2)
  lines(x = 1:sum(!is.na(raw)), y = raw[!is.na(raw)], lwd = 2, col = "red")
  
  # Projected richness:
  # Raw:
  lines(x = rep(lim-0.25, 2), y = c(proj_raw$chao-proj_raw$chao.se, proj_raw$chao+proj_raw$chao.se),
        type = "l", col = "red", lwd = 2)
  points(x = lim-0.25, y = proj_raw$chao, pch = 16, col = "red", cex = 1.25)
  # Model:
  pm <- median(proj_mod["chao",])
  psem <- median(proj_mod["chao.se",])
  lines(x = rep(lim-0.5, 2), y = c(pm-psem, pm+psem),
        type = "l", col = "black", lwd = 2)
  points(x = lim-0.5, y = pm, pch = 16, col = "black", cex = 1.25)  
  
}


# Then you can draw the species accumulation plot for any site that has more than one visit:
# The list of those sites is here (see the ESM_5 supplement for full names and locations of each of those sites:
specpool_list

# Plotting the first site on the list (simply replace the site code with the desired site):
PlotAccumYrProj(site = "ACTI")


