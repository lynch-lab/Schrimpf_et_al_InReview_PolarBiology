# Code to combine results from all different species into the objects necessary to make and analyze figures


# Libraries ---------------------------------------------------------------

library(abind) # for the abind function
library(reshape2) # for the dcast function
library(vegan) # for rarefaction techniques
library(HDInterval) # for the hdi function


# Basic Data --------------------------------------------------------------

# First load the site info:
metadata <- read.csv("ESM_5.csv", header = TRUE) # Loads site metadata, including sea-ice concentration.
metadata <- metadata[,1:6] # While this csv also includes a summary of the results for each species at each site, we're only using the sea-ice data here to run the model.



# Effort Data -------------------------------------------------------------

visit.list <- read.csv("Raw_ASI_Data.csv", header = T)

# Create a table of visits per site per year:
visit.tab <- dcast(visit.list, Code ~ Season, fun.aggregate = length)
# Ignore the message about using one column as value.var (it does what it's supposed to)

# Takes away the Code column and makes it the rownames:
rownames(visit.tab) <- visit.tab$Code
visit.tab <- visit.tab[,-1]

# Create table of yes/no visitation per year:
visit.yr.tab <- visit.tab > 0
visit.yr.tab <- as.data.frame(apply(visit.yr.tab, MARGIN = c(1,2), FUN = as.numeric))

# Create cummulative visitation per year table:
visit.tot.yr.tab <- visit.yr.tab #replicate structure
for (i in 1:nrow(visit.tot.yr.tab)) {
  for (j in 2:ncol(visit.tot.yr.tab)) {
    visit.tot.yr.tab[i,j] <- sum(visit.yr.tab[i,1:j])
  }
}



save(visit.tab, visit.yr.tab, visit.tot.yr.tab, file = "EffortData.RData")







# Compiling Parameter Results ---------------------------------------------

# Object to store the results:
param.results <- array(dim = c(3000, 17, 10), dimnames = list(NULL, c("ADPE",
                                                                      "ANSH",
                                                                      "ANTE",
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
                                                                      "WISP",
                                                                      "AnySkua"),
                                                              c("alpha.psi",
                                                                "alpha.r",
                                                                "beta.psi",
                                                                "beta.r",
                                                                "gamma.psi",
                                                                "gamma.r",
                                                                "p2",
                                                                "p3_1",
                                                                "p3_2",
                                                                "p3_3")))
# Function to extract results:
ParamCompile <- function(species, file) {
  file[,species,"alpha.psi"] <- alpha.psi.results
  file[,species,"beta.psi"] <- beta.psi.results
  file[,species,"gamma.psi"] <- gamma.psi.results
  file[,species,"alpha.r"] <- alpha.r.results
  file[,species,"beta.r"] <- beta.r.results
  file[,species,"gamma.r"] <- gamma.r.results
  file[,species,"p2"] <- p2.results
  file[,species,"p3_1"] <- p3.1.results
  file[,species,"p3_2"] <- p3.2.results
  file[,species,"p3_3"] <- p3.3.results
  return(file)
}

# Then load all species:
load("Model_Output/ADPE.out.param.Rdata")
param.results <- ParamCompile(species = "ADPE", file = param.results)

load("Model_Output/ANSH.out.param.Rdata")
param.results <- ParamCompile(species = "ANSH", file = param.results)

load("Model_Output/ANTE.out.param.Rdata")
param.results <- ParamCompile(species = "ANTE", file = param.results)

load("Model_Output/BBSP.out.param.Rdata")
param.results <- ParamCompile(species = "BBSP", file = param.results)

load("Model_Output/BRSK.out.param.Rdata")
param.results <- ParamCompile(species = "BRSK", file = param.results)

load("Model_Output/CAPE.out.param.Rdata")
param.results <- ParamCompile(species = "CAPE", file = param.results)

load("Model_Output/CHPE.out.param.Rdata")
param.results <- ParamCompile(species = "CHPE", file = param.results)

load("Model_Output/GEPE.out.param.Rdata")
param.results <- ParamCompile(species = "GEPE", file = param.results)

load("Model_Output/KEGU.out.param.Rdata")
param.results <- ParamCompile(species = "KEGU", file = param.results)

load("Model_Output/MACP.out.param.Rdata")
param.results <- ParamCompile(species = "MACP", file = param.results)

load("Model_Output/SGPE.out.param.Rdata")
param.results <- ParamCompile(species = "SGPE", file = param.results)

load("Model_Output/SNPE.out.param.Rdata")
param.results <- ParamCompile(species = "SNPE", file = param.results)

load("Model_Output/SNSH.out.param.Rdata")
param.results <- ParamCompile(species = "SNSH", file = param.results)

load("Model_Output/SOFU.out.param.Rdata")
param.results <- ParamCompile(species = "SOFU", file = param.results)

load("Model_Output/SPSK.out.param.Rdata")
param.results <- ParamCompile(species = "SPSK", file = param.results)

load("Model_Output/WISP.out.param.Rdata")
param.results <- ParamCompile(species = "WISP", file = param.results)

load("Model_Output/AnySkua.out.param.Rdata")
param.results <- ParamCompile(species = "AnySkua", file = param.results)


# Saving
save(param.results, file = "param.results.RData")



# Compiling Latent State --------------------------------------------------

# Then load all species:
load("Model_Output/ADPE.out.z.Rdata")
ADPE.z <- z.sp.results

load("Model_Output/ANSH.out.z.Rdata")
ANSH.z <- z.sp.results

load("Model_Output/ANTE.out.z.Rdata")
ANTE.z <- z.sp.results

load("Model_Output/BBSP.out.z.Rdata")
BBSP.z <- z.sp.results

load("Model_Output/BRSK.out.z.Rdata")
BRSK.z <- z.sp.results

load("Model_Output/CAPE.out.z.Rdata")
CAPE.z <- z.sp.results

load("Model_Output/CHPE.out.z.Rdata")
CHPE.z <- z.sp.results

load("Model_Output/GEPE.out.z.Rdata")
GEPE.z <- z.sp.results

load("Model_Output/KEGU.out.z.Rdata")
KEGU.z <- z.sp.results

load("Model_Output/MACP.out.z.Rdata")
MACP.z <- z.sp.results

load("Model_Output/SGPE.out.z.Rdata")
SGPE.z <- z.sp.results

load("Model_Output/SNPE.out.z.Rdata")
SNPE.z <- z.sp.results

load("Model_Output/SNSH.out.z.Rdata")
SNSH.z <- z.sp.results

load("Model_Output/SOFU.out.z.Rdata")
SOFU.z <- z.sp.results

load("Model_Output/SPSK.out.z.Rdata")
SPSK.z <- z.sp.results

load("Model_Output/WISP.out.z.Rdata")
WISP.z <- z.sp.results

load("Model_Output/AnySkua.out.z.Rdata")
AnySkua.z <- z.sp.results


lc <- dim(ADPE.z)[1] # length of all three chains
S <- dim(ADPE.z)[2] # number of sites
Y <- dim(ADPE.z)[3] # number of years

# Making an array to hold the results:
z.results <- array(dim = c(S, Y, 17, lc),
                   dimnames = list(metadata$Code,
                                   c("95-96", "96-97",
                                     "97-98", "98-99",
                                     "99-00", "00-01",
                                     "01-02", "02-03",
                                     "03-04", "04-05",
                                     "05-06", "06-07",
                                     "07-08", "08-09",
                                     "09-10", "10-11",
                                     "11-12", "12-13",
                                     "13-14", "14-15",
                                     "15-16", "16-17"),
                                   c("ADPE",
                                     "ANSH",
                                     "ANTE",
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
                                     "WISP",
                                     "AnySkua"), NULL))
# Filling the array:
for (i in 1:S) {
  for (j in 1:Y) {
    for (k in 1:lc) {
      z.results[i,j,"ADPE",k] <- ADPE.z[k,i,j]
      z.results[i,j,"ANSH",k] <- ANSH.z[k,i,j]
      z.results[i,j,"ANTE",k] <- ANTE.z[k,i,j]
      z.results[i,j,"BBSP",k] <- BBSP.z[k,i,j]
      z.results[i,j,"BRSK",k] <- BRSK.z[k,i,j]
      z.results[i,j,"CAPE",k] <- CAPE.z[k,i,j]
      z.results[i,j,"CHPE",k] <- CHPE.z[k,i,j]
      z.results[i,j,"GEPE",k] <- GEPE.z[k,i,j]
      z.results[i,j,"KEGU",k] <- KEGU.z[k,i,j]
      z.results[i,j,"MACP",k] <- MACP.z[k,i,j]
      z.results[i,j,"SGPE",k] <- SGPE.z[k,i,j]
      z.results[i,j,"SNPE",k] <- SNPE.z[k,i,j]
      z.results[i,j,"SNSH",k] <- SNSH.z[k,i,j]
      z.results[i,j,"SOFU",k] <- SOFU.z[k,i,j]
      z.results[i,j,"SPSK",k] <- SPSK.z[k,i,j]
      z.results[i,j,"WISP",k] <- WISP.z[k,i,j]
      z.results[i,j,"AnySkua",k] <- AnySkua.z[k,i,j]
    }
  }
}


# Results that represent a 1 or 0 for that state (present, or breeding):
pres.results <- (z.results > 1) * 1
breed.results <- (z.results > 2) * 1

save(pres.results, file = "pres.results.RData")
save(breed.results, file = "breed.results.RData")

# Note, since many of these files are large, you may need to clear things out the Environment, and start the next sections by loading relevant data back in.



# Editing Results ---------------------------------------------------------

# To load back the results (may need to do one at a time)
load("pres.results.RData")
load("breed.results.RData")

# Making versions that combine the the skuas into a single category (from 0 to 2)
# Presence:
tot.sk.pres <- array(NA, dim = c(dim(pres.results)[1:2], 1, dim(pres.results)[4]))
dimnames(tot.sk.pres) <- list(dimnames(pres.results)[[1]],
                               dimnames(pres.results)[[2]],
                               "TotSkua",
                               c(dimnames(pres.results)[[4]]))
for (i in 1:dim(tot.sk.pres)[1]) {
  for (j in 1:dim(tot.sk.pres)[2]) {
    for (k in 1:dim(tot.sk.pres)[4]) {
      temp <- sum(pres.results[i,j,"BRSK",k], pres.results[i,j,"SPSK",k])
      tot.sk.pres[i,j,1,k] <- max(temp, pres.results[i,j,"AnySkua",k])
    }
  }
}
pres.results.edit <- pres.results[,,!dimnames(pres.results)[[3]] %in% c("BRSK","SPSK","AnySkua"),]
pres.results.edit <- abind(pres.results.edit, tot.sk.pres, along = 3)

save(pres.results.edit, file="pres.results.edited.Rdata")



# Breeding:
tot.sk.breed <- array(NA, dim = c(dim(breed.results)[1:2], 1, dim(breed.results)[4]))
dimnames(tot.sk.breed) <- list(dimnames(breed.results)[[1]],
                               dimnames(breed.results)[[2]],
                               "TotSkua",
                               c(dimnames(breed.results)[[4]]))
for (i in 1:dim(tot.sk.breed)[1]) {
  for (j in 1:dim(tot.sk.breed)[2]) {
    for (k in 1:dim(tot.sk.breed)[4]) {
      temp <- sum(breed.results[i,j,"BRSK",k], breed.results[i,j,"SPSK",k])
      tot.sk.breed[i,j,1,k] <- max(temp, breed.results[i,j,"AnySkua",k])
    }
  }
}
breed.results.edit <- breed.results[,,!dimnames(breed.results)[[3]] %in% c("BRSK","SPSK","AnySkua"),]
breed.results.edit <- abind(breed.results.edit, tot.sk.breed, along = 3)
dim(tot.sk.breed)

save(breed.results.edit, file="breed.results.edited.Rdata")




# Summarizing Breeding Probability ----------------------------------------

# This summarizes the breeding probabiity using only those years during which each site was visited (i.e. the results that can already be found in the species columns of ESM_5).
# Note that we use the edited data, leaving the skuas as a combined category
# Although we don't do it here, this same process could be completed for the presence results.

# Load the edited data
load(file="breed.results.edited.Rdata")
# And the effort data
load(file = "EffortData.RData")

# Function to create a dataframe for any species:
SpBreedRes <- function(species) {
  # Create results object:
  res <- as.data.frame(rep(NA, times = dim(breed.results.edit)[[1]]))
  names(res) <- "Breed"
  res$Code <- dimnames(breed.results.edit)[[1]] # Site names
  # Fill the object:
  for (i in 1:length(res$Code)) {
    yr_i <- which(visit.yr.tab[i,] == 1) # year index
    post <- breed.results.edit[i,yr_i,species,] # posterior
    res$Breed[i] <- round(sum(post)/length(post), digits = 4)
  }
  colnames(res)[1] <- species
  return(res)
}

# Run it for each species:
breed.ADPE <- SpBreedRes(species = "ADPE")
breed.ANSH <- SpBreedRes(species = "ANSH")
breed.ANTE <- SpBreedRes(species = "ANTE")
breed.BBSP <- SpBreedRes(species = "BBSP")
breed.CAPE <- SpBreedRes(species = "CAPE")
breed.CHPE <- SpBreedRes(species = "CHPE")
breed.GEPE <- SpBreedRes(species = "GEPE")
breed.KEGU <- SpBreedRes(species = "KEGU")
breed.MACP <- SpBreedRes(species = "MACP")
breed.SGPE <- SpBreedRes(species = "SGPE")
breed.SNPE <- SpBreedRes(species = "SNPE")
breed.SNSH <- SpBreedRes(species = "SNSH")
breed.SOFU <- SpBreedRes(species = "SOFU")
breed.WISP <- SpBreedRes(species = "WISP")
breed.TotSkua <- SpBreedRes(species = "TotSkua") # Recall that this is between 0 and 2.

# We created ESM_5 by combining these species columns with the metadata for the sites.



# This version of the function does the same thing as above, but with the original skua categories:
SpBreedResSkua <- function(species) {
  # Create results object:
  res <- as.data.frame(rep(NA, times = dim(breed.results)[[1]]))
  names(res) <- "Breed"
  res$Code <- dimnames(breed.results)[[1]] # Site names
  # Fill the object:
  for (i in 1:length(res$Code)) {
    yr_i <- which(visit.yr.tab[i,] == 1) # year index
    post <- breed.results[i,yr_i,species,] # posterior
    res$Breed[i] <- round(sum(post)/length(post), digits = 4)
  }
  colnames(res)[1] <- species
  return(res)
}

breed.BRSK <- SpBreedResSkua(species = "BRSK")
breed.SPSK <- SpBreedResSkua(species = "SPSK")
breed.AnySkua <- SpBreedResSkua(species = "AnySkua")



# Species Accumulation ----------------------------------------------------

# This code uses the breed.results to calculate species accumulation over the course of the model history. For each site, only the years in which that site was visited are used. Note that this uses the original, not the edited version (see the section about combining skuas below).

# The first step is to load the breed.results back in, if they haven't been entered already:
load("breed.results.RData")

accum.results <- breed.results # replicate data structure (with names)

# This takes a while:
for (i in 1:dim(accum.results)[1]) { # Site
  yrs <- which(visit.tab[i,] > 0) # index of first visit year
  for (j in 1:dim(accum.results)[2]) { # Year
    for (k in 1:dim(accum.results)[3]) { # Species
      for (l in 1:dim(accum.results)[4]) { # Iteration
        if (!(j %in% yrs)) { # when this year is not a visited year
          accum.results[i,j,k,l] <- 0
        } else {# By choosing 0 instead of NA, I make the code easier later, but also make it so that 0 can't be interpreted as evidence for absence
        m <- which(yrs == j)
        if (any(breed.results[i,yrs[1:m],k,l]==1)) { # if any of the visited years showed that species as breeding (in that iteration)
          accum.results[i,j,k,l] <- 1
        } else {
          accum.results[i,j,k,l] <- 0
        }
        }

      }
    }
  }
}

# Combining Skuas
# We chose to combine skuas after calculating the species accumulation because if we did it before, skuas could conceivably show up in alternating years, giving us an underestimate of the total skua species in some years.
tot.sk <- array(NA, dim = c(dim(accum.results)[1:2], 1, dim(accum.results)[4]))
dimnames(tot.sk) <- list(dimnames(accum.results)[[1]], dimnames(accum.results)[[2]], "TotSkua", c(dimnames(accum.results)[[4]]))
for (i in 1:dim(tot.sk)[1]) {
  for (j in 1:dim(tot.sk)[2]) {
    for (k in 1:dim(tot.sk)[4]) {
      temp <- sum(accum.results[i,j,"BRSK",k], accum.results[i,j,"SPSK",k])
      tot.sk[i,j,1,k] <- max(temp, accum.results[i,j,"AnySkua",k])
    }
  }
}

accum.results <- accum.results[,,!dimnames(accum.results)[[3]] %in% c("BRSK","SPSK","AnySkua"),]
accum.results <- abind(accum.results, tot.sk, along = 3)

# Save the results:
save(accum.results,file="accum.results.Rdata")


