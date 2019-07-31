# Code to run the occupancy model for a single species


# Libraries ---------------------------------------------------------------

jags.path="/usr/local/bin/jags/"
library('R2jags') # To run JAGS
library('coda') # To manipulate JAGS output
library('gtools')
library('boot')
library('reshape')
library('reshape2') # used for acast function
library('snowfall') # used to run in parallel
library('parallel') # used to run in parallel
library('SpatialEpi') # used to run in parallel



# Load Data ---------------------------------------------------------------

metadata <- read.csv("ESM_5.csv", header = TRUE) # Loads site metadata, including sea-ice concentration.
metadata <- metadata[,1:6] # NWhile this csv also includes a summary of the results for each species at each site, we're only using the sea-ice data here to run the model.

visitdata <- read.csv("Raw_ASI_Data.csv", header = TRUE) # Loads raw ASI visit data

fulldata <- merge(visitdata, metadata, by="Code", all.x=TRUE)
fulldata <- fulldata[order(fulldata$Code, fulldata$Season),]


# Preparing Data ----------------------------------------------------------

# To change the species that the model uses:
species <- "GEPE"
# Species Codes:
# Adélie penguin: ADPE
# Antarctic shag: ANSH
# Antarctic tern: ANTE
# black-bellied storm-petrel: BBSP
# brown skua: BRSK
# cape petrel: CAPE
# chinstrap penguin: CHPE
# gentoo penguin: GEPE
# kelp gull: KEGU
# macaroni penguin: MACP
# snow petrel: SNPE
# snowy sheathbill: SNSH
# south polar skua: SPSK
# southern fulmar: SOFU
# southern giant petrel: SGPE
# Wilson's storm-petrel: WISP


# Creating the data frame to be used:
data <- data.frame(fulldata$Code) # Site Codes
colnames(data) <- "Code"
data$Season <- fulldata$Season # Field season
data$state <- as.numeric(fulldata[,species]+1) # codes the data as: 1 (not detected), 2 (presence), 3 (presence and breeding)
data$seaice <- as.numeric(fulldata$AvgNovSIC) # Covariate (sea ice concentration)

# Labels the number of repeat visits in a year
data$visit <- 1
for (i in 2:dim(data)[1]){
  if (data$Code[i] == data$Code[i-1] &
      data$Season[i] == data$Season[i-1])
    data$visit[i] <- data$visit[i-1]+1
}

# turning variables into numeric for ease in modeling
data$site <- as.numeric(as.factor(data$Code))
data$season <- as.numeric(as.factor(data$Season))
final_data <- subset(data, select=c(site,season,visit,state))

# Creating sea-ice object (single row for each site)
unique_seaice <- unique(subset(data, select=c(site,seaice)))
seaice <- as.numeric(scale(unique_seaice$seaice)) # Normalizes the covariate


V <- max(final_data$visit) # Highest number of visits
S <- length(unique(final_data$site)) # Number of sites
Y <- length(unique(final_data$season)) # Number of Years
y <- acast(final_data, site ~ season ~ visit) # Data array for model




# Model Code --------------------------------------------------------------

# Creates the .jags file for JAGS to read:

sink("Occupancy.jags")
cat("
    
    model {
    
    #### Priors ####
    
    beta.psi ~ dnorm(0,0.386)
    beta.r ~ dnorm(0,0.386)
    alpha.psi ~ dnorm(0,0.386)
    alpha.r ~ dnorm(0,0.386)
    gamma.psi ~ dnorm(0,0.386)
    gamma.r ~ dnorm(0,0.386)
    
    p2 ~ dunif(0,1)
    p3[1:3] ~ ddirch(alpha.p3[1:3])
    for (i in 1:3){
    alpha.p3[i] <- 1
    }
    
    
    #### Likelihood ####
    
    ## states-space model ##
    
    # First year:
    for (i in 1:S){
    logit(mu.psi[i,1]) <- alpha.psi + beta.psi*seaice[i]
    logit(mu.r[i,1]) <- alpha.r + beta.r*seaice[i]
    z1[i,1] ~ dbin(mu.psi[i,1], 1)
    z2[i,1] ~ dbin(mu.r[i,1]*z1[i,1], 1)
    z[i,1] <- c(1-max(z1[i,1], z2[i,1]),
    z1[i,1]-z2[i,1],
    z2[i,1]
    ) %*% c(1,2,3)
    }
    
    # Subsequent years:
    for (i in 1:S){
    for (j in 2:Y){
    logit(mu.psi[i,j]) <- alpha.psi +
    beta.psi*seaice[i] +
    gamma.psi*z1[i,j-1]
    logit(mu.r[i,j]) <- alpha.r + 
    beta.r*seaice[i] + 
    gamma.r*z2[i,j-1]
    z1[i,j] ~ dbin(mu.psi[i,j], 1)
    z2[i,j] ~ dbin(mu.r[i,j]*z1[i,j], 1)
    z[i,j] <- c(1-max(z1[i,j], z2[i,j]),
    z1[i,j]-z2[i,j],
    z2[i,j]
    ) %*% c(1,2,3)
    }
    }
    
    
    
    ## observation model ##
    
    # Define observation matrix
    
    p[1,1] <- 1
    p[1,2] <- 0
    p[1,3] <- 0
    p[2,1] <- 1-p2
    p[2,2] <- p2
    p[2,3] <- 0
    p[3,1] <- p3[1]
    p[3,2] <- p3[2]
    p[3,3] <- p3[3]
    
    
    for (i in 1:S){
    for (j in 1:Y){
    for (k in 1:V){
    y[i,j,k] ~ dcat(p[z[i,j],])
    y.new[i,j,k] ~ dcat(p[z[i,j],]) # for posterior predictive check
    }
    }
    }
    
    
    #### Derived values ####
    
    ## phi = the probability of being in each state
    for (i in 1:S) {
    for (j in 1:Y) {
    phi[i,j,1] <- 1 - mu.psi[i,j]
    phi[i,j,2] <- mu.psi[i,j] * (1 - mu.r[i,j])
    phi[i,j,3] <- mu.psi[i,j] * mu.r[i,j]
    }
    }
    
    ## eval = the probability of recording each state
    for (i in 1:S) {
    for (j in 1:Y) {
    eval[i,j,1] <- (phi[i,j,1] * 1) + 
    (phi[i,j,2] * (1 - p2)) + 
    (phi[i,j,3] * p3[1])
    eval[i,j,2] <- (phi[i,j,1] * 0) + 
    (phi[i,j,2] * p2) + 
    (phi[i,j,3] * p3[2])
    eval[i,j,3] <- (phi[i,j,1] * 0) + 
    (phi[i,j,2] * 0) + 
    (phi[i,j,3] * p3[3])
    }
    }
    
    
    
    }",fill = TRUE)
sink()


# Running -----------------------------------------------------------------

#####  Bundle the data for JAGS #####

Dat <- list(y = y, S = dim(y)[1], Y = dim(y)[2], V = dim(y)[3], seaice = seaice)



##### Select Initial Values #####

InitStage <- function()
{
  list(
    z1 = matrix(1, nrow = dim(y)[1], ncol = dim(y)[2]),
    z2 = matrix(1, nrow = dim(y)[1], ncol = dim(y)[2]),
    beta.psi = runif(1,0.2,0.2),
    beta.r = runif(1,0.2,0.2),
    alpha.psi = runif(1, 0.2, 0.3),
    alpha.r = runif(1, 0.2, 0.3),
    gamma.psi = 0,
    gamma.r = 0,
    p3 = as.vector(rdirichlet(1,c(1,1,1))),
    p2 = runif(1,.1,.9)
  )
}



##### Posteriors to return #####

# Note: it may be necessary to limit this list to limit object size
# (recommend only returning "y.new" OR "z")

ParsStage <- c("beta.psi","beta.r",
               "alpha.psi","alpha.r",
               "gamma.psi","gamma.r",
               "p2","p3",
               "mu.psi","mu.r",
               "eval","y.new",
               "z")



##### MCMC settings #####

ni <- 300000 # Total draws
nt <- 50 # thinning rate
nb <- 250000 # burn-in
nc <- 3 # number of chains
n.adapt <- 10000



##### Parallel Execution #####

cl <- makeCluster(3)

clusterExport(cl, c("Dat", "y", "InitStage", "ParsStage", "n.adapt", "nb", "ni", "nt"))

system.time({
  out1 <- clusterEvalQ(cl, {
    library(rjags); library(SpatialEpi); library('gtools')
    jm <- jags.model("Occupancy.jags",
                     data = Dat, InitStage,
                     n.chains=1, n.adapt=n.adapt)
    update(jm, n.iter=nb, thin=nt)
    samples = coda.samples(jm, n.iter=ni, variable.names=ParsStage, thin=nt)
    return(as.mcmc(samples))
  })
})

stopCluster(cl)


zm <- mcmc.list(out1)

# # If you want to save the whole model output, use this, but usually it is too large to effectively save:
# save(zm,file = paste(species,".results.Rdata",sep = ""))


# Packaging Results -------------------------------------------------------

lc <- (ni - nb) / nt * nc # total number of draws from the posterior (i.e. the length of all three chains put together).
# # That should be the same as:
# length(zm$JAGSoutput$sims.array[,1:3,1])


##### The latent-tate results #####

z.sp.results <- array(dim = c(lc, S, Y))
for (i in 1:S) {
  for (j in 1:Y) {
    z.sp.results[,i,j] <- c(zm$JAGSoutput$sims.array[,1,which(names(zm$JAGSoutput$sims.array[1,1,]) == paste("z[",as.character(i),",",as.character(j),"]",sep=""))],
                            zm$JAGSoutput$sims.array[,2,which(names(zm$JAGSoutput$sims.array[1,2,]) == paste("z[",as.character(i),",",as.character(j),"]",sep=""))],
                            zm$JAGSoutput$sims.array[,3,which(names(zm$JAGSoutput$sims.array[1,3,]) == paste("z[",as.character(i),",",as.character(j),"]",sep=""))])
  }
}

save(z.sp.results, file = paste(species,".out.z.Rdata",sep=""))



##### Parameter results #####

alpha.psi.results <- c(zm$JAGSoutput$sims.array[,1,"alpha.psi"],
                       zm$JAGSoutput$sims.array[,2,"alpha.psi"],
                       zm$JAGSoutput$sims.array[,3,"alpha.psi"])
alpha.r.results <- c(zm$JAGSoutput$sims.array[,1,"alpha.r"],
                     zm$JAGSoutput$sims.array[,2,"alpha.r"],
                     zm$JAGSoutput$sims.array[,3,"alpha.r"])
beta.psi.results <- c(zm$JAGSoutput$sims.array[,1,"beta.psi"],
                      zm$JAGSoutput$sims.array[,2,"beta.psi"],
                      zm$JAGSoutput$sims.array[,3,"beta.psi"])
beta.r.results <- c(zm$JAGSoutput$sims.array[,1,"beta.r"],
                    zm$JAGSoutput$sims.array[,2,"beta.r"],
                    zm$JAGSoutput$sims.array[,3,"beta.r"])
gamma.psi.results <- c(zm$JAGSoutput$sims.array[,1,"gamma.psi"],
                       zm$JAGSoutput$sims.array[,2,"gamma.psi"],
                       zm$JAGSoutput$sims.array[,3,"gamma.psi"])
gamma.r.results <- c(zm$JAGSoutput$sims.array[,1,"gamma.r"],
                     zm$JAGSoutput$sims.array[,2,"gamma.r"],
                     zm$JAGSoutput$sims.array[,3,"gamma.r"])
p2.results <- c(zm$JAGSoutput$sims.array[,1,"p2"],
                zm$JAGSoutput$sims.array[,2,"p2"],
                zm$JAGSoutput$sims.array[,3,"p2"])
p3.1.results <- c(zm$JAGSoutput$sims.array[,1,"p3[1]"],
                  zm$JAGSoutput$sims.array[,2,"p3[1]"],
                  zm$JAGSoutput$sims.array[,3,"p3[1]"])
p3.2.results <- c(zm$JAGSoutput$sims.array[,1,"p3[2]"],
                  zm$JAGSoutput$sims.array[,2,"p3[2]"],
                  zm$JAGSoutput$sims.array[,3,"p3[2]"])
p3.3.results <- c(zm$JAGSoutput$sims.array[,1,"p3[3]"],
                  zm$JAGSoutput$sims.array[,2,"p3[3]"],
                  zm$JAGSoutput$sims.array[,3,"p3[3]"])

save(alpha.psi.results,
     alpha.r.results,
     beta.psi.results,
     beta.r.results,
     gamma.psi.results,
     gamma.r.results,
     p2.results,
     p3.1.results,
     p3.2.results,
     p3.3.results,
     file = paste(species,".out.param.Rdata",sep=""))



##### Info for the Posterior Predictive Check #####

y.new.results <- array(dim = c(lc, S, Y, V))
for (i in 1:S) {
  for (j in 1:Y) {
    for (k in 1:V) {
      y.new.results[,i,j,k] <- c(zm$JAGSoutput$sims.array[,1,which(names(zm$JAGSoutput$sims.array[1,1,]) == (paste("y.new[",as.character(i),",",as.character(j),",",as.character(k),"]",sep="")))],
                                 zm$JAGSoutput$sims.array[,2,which(names(zm$JAGSoutput$sims.array[1,2,]) == (paste("y.new[",as.character(i),",",as.character(j),",",as.character(k),"]",sep="")))],
                                 zm$JAGSoutput$sims.array[,3,which(names(zm$JAGSoutput$sims.array[1,3,]) == (paste("y.new[",as.character(i),",",as.character(j),",",as.character(k),"]",sep="")))])
    }
  }
}

for (i in 1:S) {
  for (j in 1:Y) {
    for (k in 1:V) {
      if (is.na(y[i,j,k])) y.new.results[,i,j,k] <- NA
    }
  }
}


eval.results <- array(dim = c(lc, S, Y, 3))
for (i in 1:S) {
  for (j in 1:Y) {
    for (k in 1:3) {
      eval.results[,i,j,k] <- c(zm$JAGSoutput$sims.array[,1,which(names(zm$JAGSoutput$sims.array[1,1,]) == (paste("eval[",as.character(i),",",as.character(j),",",as.character(k),"]",sep="")))],
                                zm$JAGSoutput$sims.array[,2,which(names(zm$JAGSoutput$sims.array[1,2,]) == (paste("eval[",as.character(i),",",as.character(j),",",as.character(k),"]",sep="")))],
                                zm$JAGSoutput$sims.array[,3,which(names(zm$JAGSoutput$sims.array[1,3,]) == (paste("eval[",as.character(i),",",as.character(j),",",as.character(k),"]",sep="")))])
    }
  }
}

# Note, these objects don't need to be saved after calculating t.pred and t.actual below, but they can be as a back-up:
save(y.new.results,
     eval.results,
     file = paste(species,".out.ppc.Rdata",sep=""))


##### Posterior Predicitive Check test-statistics #####
# Test statistic for predicted data:
t.pred <- rep(NA, times = lc)
for (i in 1:lc) {
  for (j in 1:S) {
    for (k in 1:Y) {
      for (l in 1:V) {
        temp.y.vect <- 
          c(as.numeric(y.new.results[i,j,k,l]==1),
            as.numeric(y.new.results[i,j,k,l]==2),
            as.numeric(y.new.results[i,j,k,l]==3))
        t.pred[i] <- sum(t.pred[i], sum((temp.y.vect - eval.results[i,j,k,]) ^2 /
                                          eval.results[i,j,k,], na.rm=T), na.rm=T)
      }
    }
  }
}

# Test statistic for actual data:
t.actual <- rep(NA, times = lc)
for (i in 1:lc) {
  for (j in 1:S) {
    for (k in 1:Y) {
      for (l in 1:V) {
        temp.y.vect <- c(as.numeric(y[j,k,l]==1),
                         as.numeric(y[j,k,l]==2),
                         as.numeric(y[j,k,l]==3))
        t.actual[i] <- sum(t.actual[i], sum((temp.y.vect - eval.results[i,j,k,]) ^2 /
                                              eval.results[i,j,k,], na.rm=T), na.rm=T)
      }
    }
  }
}


save(t.pred, t.actual, file=paste(species,".out.ppc2.Rdata",sep=""))













