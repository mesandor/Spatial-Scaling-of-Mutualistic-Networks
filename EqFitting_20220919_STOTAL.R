##### Model comparison for IAR (Equation 9) and other NARs (Equation 10) ######

# Load packages
library(R2jags)
library(R2WinBUGS)
library(MCMCvis)

# Read data
AllNetworks <- read.csv("AllStudiesNetworkData.csv")
colnames(AllNetworks) <- c("Patch_ID", "Patch_Area_m2", "Lower_Spp_Rich", "Upper_Spp_Rich", "No_Int", "Connectance", "Dataset", "Type", "Type2", "Connectivity")


# Subset Datasets
Seeds <- subset(AllNetworks, AllNetworks$Dataset=="Sandor")
Grass1 <- subset(AllNetworks, AllNetworks$Dataset=="Grass_pollin")
Sugiura <- subset(AllNetworks, AllNetworks$Dataset=="Sugiura")


#################### Write models #################### 
# A first model that fits Connectance-area relationship model and calculates parameters
CARs <- function(){
  alpha ~ dnorm(0, 0.01) 
  delta ~ dnorm(0, 0.01) 
  sigma1 ~ dunif(0, 100)
  tau1 <- 1/(sigma1*sigma1)
  
  for(i in 1:n.points){    
    mu1[i] <- alpha + delta * log.Area[i] #Connectance-area relationship
    log.connectance[i] ~ dnorm(mu1[i], tau1)
  }
}

# A second model that calculates SAR parameters
SARs <- function(){
  
  log.c1 ~ dnorm(0, 0.01) 
  log.c2 ~ dnorm(0, 0.01) 
  log.ct ~ dnorm(0, 0.01) 
  z1 ~ dnorm(0, 0.01) 
  z2 ~ dnorm(0, 0.01) 
  zt ~ dnorm(0, 0.01) 
  # u is given as data as assumption of both models
  sigma1 ~ dunif(0, 100)
  tau1 <- 1/(sigma1*sigma1)
  sigma2 ~ dunif(0, 100)
  tau2 <- 1/(sigma2*sigma2)
  sigma3 ~ dunif(0, 100)
  tau3 <- 1/(sigma3*sigma3)
  
  for(i in 1:n.points){
    mu1[i] <- log.c1 + z1 * log.Area[i] #SAR for spp group 1
    log.Spp1[i] ~ dnorm(mu1[i], tau1)
    mu2[i] <- log.c2 + z2 * log.Area[i] #SAR for spp group 2
    log.Spp2[i] ~ dnorm(mu2[i], tau2)
    mu3[i] <- log.ct + zt * log.Area[i] #SAR for total spp
    log.Sppt[i] ~ dnorm(mu3[i], tau3)
  }
}



#################### Fit CAR models to all 3 datasets #################### 
parameters.CAR <- c("alpha", "delta") 

seeds.data.car <- list(
  n.points = dim(Seeds)[1],
  log.connectance = log10(Seeds$Connectance/100),
  log.Area = log10(Seeds$Patch_Area_m2)
)
grass.data.car <- list(
  n.points = dim(Grass1)[1],
  log.connectance = log10(Grass1$Connectance/100),
  log.Area = log10(Grass1$Patch_Area_m2)
)
ant.data.car <- list(
  n.points = dim(Sugiura)[1],
  log.connectance = log10(Sugiura$Connectance/100),
  log.Area = log10(Sugiura$Patch_Area_m2)
)

CAR.S <- jags(data=seeds.data.car, inits = list(list(alpha=1, delta=0), list(alpha=-1, delta=1), list(alpha=3, delta=-1)), parameters.to.save= parameters.CAR, CARs, n.burnin=480000, n.iter=500000, n.thin=100, n.chains=3)
CAR.G <- jags(data=grass.data.car, inits = list(list(alpha=1, delta=0), list(alpha=-1, delta=1), list(alpha=3, delta=-1)), parameters.to.save= parameters.CAR, CARs, n.burnin=480000, n.iter=500000, n.thin=100, n.chains=3)
CAR.A <- jags(data=ant.data.car, inits = list(list(alpha=1, delta=0), list(alpha=1, delta=1), list(alpha=3, delta=-1)), parameters.to.save= parameters.CAR, CARs, n.burnin=480000, n.iter=500000, n.thin=100, n.chains=3)
# Remember to check fits, Rhat, N.eff
CAR.S
CAR.G
CAR.A


#################### Fit SAR models to all 3 datasets #################### 
parameters.SAR <- c("log.c1", "log.c2", "log.ct", "z1", "z2", "zt") 

seeds.data.sar <- list(
  n.points = dim(Seeds)[1],
  log.Spp1 = log10(Seeds$Lower_Spp_Rich),
  log.Spp2 = log10(Seeds$Upper_Spp_Rich),
  log.Sppt = log10((Seeds$Upper_Spp_Rich + Seeds$Lower_Spp_Rich)),
  log.Area = log10(Seeds$Patch_Area_m2)
)
grass.data.sar <- list(
  n.points = dim(Grass1)[1],
  log.Spp1 = log10(Grass1$Lower_Spp_Rich),
  log.Spp2 = log10(Grass1$Upper_Spp_Rich),
  log.Sppt = log10((Grass1$Upper_Spp_Rich + Grass1$Lower_Spp_Rich)),
  log.Area = log10(Grass1$Patch_Area_m2)
)
ant.data.sar <- list(
  n.points = dim(Sugiura)[1],
  log.Spp1 = log10(Sugiura$Lower_Spp_Rich),
  log.Spp2 = log10(Sugiura$Upper_Spp_Rich),
  log.Sppt = log10((Sugiura$Upper_Spp_Rich + Sugiura$Lower_Spp_Rich)),
  log.Area = log10(Sugiura$Patch_Area_m2)
)

SAR.S <- jags(data=seeds.data.sar, inits = list(list(log.c1=3, log.c2=3, log.ct=3, zt=0.33, z1=0.3, z2=0.35), list(log.c1=1, log.c2=1, log.ct=1, zt=0.23, z1=0.2, z2=0.25), list(log.c1=-1, log.c2=-1, log.ct=-1, zt=0.13, z1=0.1, z2=0.15)), parameters.to.save= parameters.SAR, SARs, n.burnin=480000, n.iter=500000, n.thin=100, n.chains=3)
SAR.G <- jags(data=grass.data.sar, inits = list(list(log.c1=3, log.c2=3, log.ct=3, zt=0.33, z1=0.3, z2=0.35), list(log.c1=1, log.c2=1, log.ct=1, zt=0.23, z1=0.2, z2=0.25), list(log.c1=-1, log.c2=-1, log.ct=-1, zt=0.13, z1=0.1, z2=0.15)), parameters.to.save= parameters.SAR, SARs, n.burnin=480000, n.iter=500000, n.thin=100, n.chains=3)
SAR.A <- jags(data=ant.data.sar, inits = list(list(log.c1=3, log.c2=3, log.ct=3, zt=0.33, z1=0.3, z2=0.35), list(log.c1=1, log.c2=1, log.ct=1, zt=0.23, z1=0.2, z2=0.25), list(log.c1=-1, log.c2=-1, log.ct=-1, zt=0.13, z1=0.1, z2=0.15)), parameters.to.save= parameters.SAR, SARs, n.burnin=480000, n.iter=500000, n.thin=100, n.chains=3)
# Remember to check fits, Rhat, N.eff
SAR.S
SAR.G
SAR.A

################################# Create a fig to compare fit SARs to data ########################
x.S <- log10(Seeds$Patch_Area_m2)
x.G <- log10(Grass1$Patch_Area_m2)
x.A <- log10(Sugiura$Patch_Area_m2)

y.S <- log10((Seeds$Upper_Spp_Rich+Seeds$Lower_Spp_Rich))
y.G <- log10((Grass1$Upper_Spp_Rich+Grass1$Lower_Spp_Rich))
y.A <- log10((Sugiura$Upper_Spp_Rich+Sugiura$Lower_Spp_Rich))

pdf(file = paste0("FigureS3_STOTAL_",Sys.Date(),".pdf"), width = 3.25, height = 6)
par(las = 1, mar = c(2, 3, 0.5, 0.5), oma=c(1,0,0,0),
    tcl = -0.2, mgp = c(1.8,0.25,0), mfrow = c(3, 1))

plot(x.S, y.S, ylim = c(0,2),
     xlab = "", ylab = "log(Species)", 
     pch=21, col="black", bg="#00000080", cex = 1.5)
abline(a=mean(MCMCpstr(SAR.S, params = parameters.SAR, type = "chains")$log.ct), b=mean(MCMCpstr(SAR.S, params = parameters.SAR, type = "chains")$zt), lwd=2)
title(main = "A. Seed dispersal", adj = 0.01, line = -1.2)

plot(x.G, y.G, ylim = c(0.5,2.5),
     xlab = "", ylab = "log(Species)", 
     pch=21, col="black", bg="#00000080", cex = 1.5)
abline(a=mean(MCMCpstr(SAR.G, params = parameters.SAR, type = "chains")$log.ct), b=mean(MCMCpstr(SAR.G, params = parameters.SAR, type = "chains")$zt), lwd=2)
title(main = "B. Pollination", adj = 0.01, line = -1.2)

plot(x.A, y.A, ylim = c(0.5,2.5),
     xlab = "", ylab = "log(Species)", 
     pch=21, col="black", bg="#00000080", cex = 1.5)
abline(a=mean(MCMCpstr(SAR.A, params = parameters.SAR, type = "chains")$log.ct), b=mean(MCMCpstr(SAR.A, params = parameters.SAR, type = "chains")$zt), lwd=2)
title(main = "C. Ant-Plant", adj = 0.01, line = -1.2)
mtext(text = expression(paste("log(area) ", (m^2))), side = 1, line = 1.5, cex = 2/3)

dev.off()



################################## Save interactions for later #####################################
seeds.data.int <- list(
  log.interacts = log10(Seeds$No_Int)
)

grass.data.int <- list(
  log.interacts = log10(Grass1$No_Int)
)

ant.data.int <- list(
  log.interacts = log10(Sugiura$No_Int)
)




#################### Calculate IAR for all 3 datasets #################### 
# x-value sequences
x.area.S <- seq(min(seeds.data.sar$log.Area), max(seeds.data.sar$log.Area), length.out = 100)
x.area.G <- seq(min(grass.data.sar$log.Area), max(grass.data.sar$log.Area), length.out = 100)
x.area.A <- seq(min(ant.data.sar$log.Area), max(ant.data.sar$log.Area), length.out = 100)

# extract parameters
param.S.SAR <- MCMCpstr(SAR.S, params = parameters.SAR, type = "chains")
param.G.SAR <- MCMCpstr(SAR.G, params = parameters.SAR, type = "chains")
param.A.SAR <- MCMCpstr(SAR.A, params = parameters.SAR, type = "chains")

# extract parameters
param.S.CAR <- MCMCpstr(CAR.S, params = parameters.CAR, type = "chains")
param.G.CAR <- MCMCpstr(CAR.G, params = parameters.CAR, type = "chains")
param.A.CAR <- MCMCpstr(CAR.A, params = parameters.CAR, type = "chains")

# create parameters
set.seed(412)
param.S.LSCC <- as.data.frame(cbind(log10(rnorm(length(param.S.SAR$log.c2), 2, sd(Seeds$Connectance/100))), log10(rnorm(length(param.S.SAR$log.c2), mean(Seeds$Connectance/100), sd(Seeds$Connectance/100)))))
colnames(param.S.LSCC) <- c("log.b1", "log.b2")
param.G.LSCC <- as.data.frame(cbind(log10(rnorm(length(param.G.SAR$log.c2), 2, sd(Grass1$Connectance/100))), log10(rnorm(length(param.G.SAR$log.c2), mean(Grass1$Connectance/100), sd(Grass1$Connectance/100)))))
colnames(param.G.LSCC) <- c("log.b1", "log.b2")
param.A.LSCC <- as.data.frame(cbind(log10(rnorm(length(param.A.SAR$log.c2), 2, sd(Sugiura$Connectance/100))), log10(rnorm(length(param.A.SAR$log.c2), mean(Sugiura$Connectance/100), sd(Sugiura$Connectance/100)))))
colnames(param.A.LSCC) <- c("log.b1", "log.b2")

# create a storage object that is the dimensions of x-values by posterior
pred.IAR.S <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.S)))
pred.IAR.G <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.G)))
pred.IAR.A <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.A)))

pred.LSSL.S <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.S)))
pred.LSSL.G <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.G)))
pred.LSSL.A <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.A)))

pred.CCH.S <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.S)))
pred.CCH.G <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.G)))
pred.CCH.A <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.A)))

# make predicted line for each posterior draw
#IAR = alpha + log.c1 + log.c2 + (delta + z1 + z2)*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.IAR.S[i, ] <- param.S.CAR$alpha[i] + param.S.SAR$log.c1[i] + param.S.SAR$log.c2[i] + (param.S.CAR$delta[i] + param.S.SAR$z1[i] + param.S.SAR$z2[i]) * x.area.S
  pred.IAR.G[i, ] <- param.G.CAR$alpha[i] + param.G.SAR$log.c1[i] + param.G.SAR$log.c2[i] + (param.G.CAR$delta[i] + param.G.SAR$z1[i] + param.G.SAR$z2[i]) * x.area.G
  pred.IAR.A[i, ] <- param.A.CAR$alpha[i] + param.A.SAR$log.c1[i] + param.A.SAR$log.c2[i] + (param.A.CAR$delta[i] + param.A.SAR$z1[i] + param.A.SAR$z2[i]) * x.area.A
}

# make predicted line for each posterior draw
#Eq 10, b ~2, u=1, LSSL
#LSSL = log(b) + u*log(c2) + u*z2*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.LSSL.S[i, ] <- param.S.LSCC$log.b1[i] + 1 * param.S.SAR$log.ct[i]  + 1 * param.S.SAR$zt[i] * x.area.S
  pred.LSSL.G[i, ] <- param.G.LSCC$log.b1[i] + 1 * param.G.SAR$log.ct[i]  + 1 * param.G.SAR$zt[i] * x.area.G
  pred.LSSL.A[i, ] <- param.A.LSCC$log.b1[i] + 1 * param.A.SAR$log.ct[i]  + 1 * param.A.SAR$zt[i] * x.area.A
}

# make predicted line for each posterior draw
#Eq 10, 0 < 1 < b, u=2, CCH
#CCH = log(b) + u*log(c2) + u*z2*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.CCH.S[i, ] <- param.S.LSCC$log.b2[i] + 2 * param.S.SAR$log.ct[i]  + 2 * param.S.SAR$zt[i] * x.area.S
  pred.CCH.G[i, ] <- param.G.LSCC$log.b2[i] + 2 * param.G.SAR$log.ct[i]  + 2 * param.G.SAR$zt[i] * x.area.G
  pred.CCH.A[i, ] <- param.A.LSCC$log.b2[i] + 2 * param.A.SAR$log.ct[i]  + 2 * param.A.SAR$zt[i] * x.area.A
}


#################### Plot Figure S2 #################### 

pdf(file = paste0("Figure2_STOTAL_",Sys.Date(),".pdf"), width = 3.25, height = 6)
par(las = 1, mar = c(2, 3, 0.5, 0.5), oma=c(1,0,0,0),
    tcl = -0.2, mgp = c(1.8,0.25,0), mfrow = c(3, 1))

plot(seeds.data.sar$log.Area, seeds.data.int$log.interacts, ylim = c(0.3, 2),
     xlab = "", ylab = "log(Interactions)", 
     pch=21, col="black", bg="#00000080", cex = 1.5)
lines(x.area.S, apply(pred.LSSL.S, 2, mean), col = "#929084", lwd = 2, lty = 4)
lines(x.area.S, apply(pred.CCH.S, 2, mean), col = "#2e4052", lwd = 2, lty = 2)
lines(x.area.S, apply(pred.IAR.S, 2, mean), col = "#e5323b", lwd = 2, lty = 1)
title(main = "A. Seed dispersal", adj = 0.01, line = -1.2)

plot(grass.data.sar$log.Area, grass.data.int$log.interacts, ylim = c(1, 3),
     xlab = "", ylab = "log(Interactions)", 
     pch=22, col="black", bg="#00000080", cex = 1.5)
lines(x.area.G, apply(pred.LSSL.G, 2, mean), col = "#929084",  lwd = 2, lty = 4)
lines(x.area.G, apply(pred.CCH.G, 2, mean), col = "#2e4052", lwd = 2, lty = 2)
lines(x.area.G, apply(pred.IAR.G, 2, mean), col = "#e5323b", lwd = 2, lty = 1)
title(main = "B. Pollination", adj = 0.01, line = -1.2)

plot(ant.data.sar$log.Area, ant.data.int$log.interacts, ylim = c(0.5, 2.6),
     xlab = "", ylab = "log(Interactions)", 
     pch=23, col="black", bg="#00000080", cex = 1.5)
lines(x.area.A, apply(pred.LSSL.A, 2, mean), col = "#929084", lwd = 2, lty = 4)
lines(x.area.A, apply(pred.CCH.A, 2, mean), col = "#2e4052", lwd = 2, lty = 2)
lines(x.area.A, apply(pred.IAR.A, 2, mean), col = "#e5323b", lwd = 2, lty = 1)
title(main = "C. Ant-Plant", adj = 0.01, line = -1.2)
legend(x = "bottomright", lty = c(4, 2, 1), lwd = 2, col = c("#929084", "#2e4052", "#e5323b"),
       legend = c("Link-species scaling law", "Constant connectance hypothesis", "Interaction-area relationship"))
mtext(text = expression(paste("log(area) ", (m^2))), side = 1, line = 1.5, cex = 2/3)

dev.off()



######################### CALCULATE RMSE ###############################
NS <- length(seeds.data.int$log.interacts)
NG <- length(grass.data.int$log.interacts)
NAP <- length(ant.data.int$log.interacts)

# x-value sequences
x.area.S.real <- seeds.data.sar$log.Area
x.area.G.real <- grass.data.sar$log.Area
x.area.A.real <- ant.data.sar$log.Area

# create a storage object that is the dimensions of x-values by posterior
pred.IAR.S.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.S.real)))
pred.IAR.G.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.G.real)))
pred.IAR.A.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.A.real)))

pred.LSSL.S.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.S.real)))
pred.LSSL.G.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.G.real)))
pred.LSSL.A.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.A.real)))

pred.CCH.S.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.S.real)))
pred.CCH.G.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.G.real)))
pred.CCH.A.real <- array(dim = c(length(param.S.SAR$log.c1), length(x.area.A.real)))


# make predicted line for each posterior draw
#IAR = alpha + log.c1 + log.c2 + (delta + z1 + z2)*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.IAR.S.real[i, ] <- param.S.CAR$alpha[i] + param.S.SAR$log.c1[i] + param.S.SAR$log.c2[i] + (param.S.CAR$delta[i] + param.S.SAR$z1[i] + param.S.SAR$z2[i]) * x.area.S.real
  pred.IAR.G.real[i, ] <- param.G.CAR$alpha[i] + param.G.SAR$log.c1[i] + param.G.SAR$log.c2[i] + (param.G.CAR$delta[i] + param.G.SAR$z1[i] + param.G.SAR$z2[i]) * x.area.G.real
  pred.IAR.A.real[i, ] <- param.A.CAR$alpha[i] + param.A.SAR$log.c1[i] + param.A.SAR$log.c2[i] + (param.A.CAR$delta[i] + param.A.SAR$z1[i] + param.A.SAR$z2[i]) * x.area.A.real
}

# make predicted line for each posterior draw
#Eq 10, u=1, LSSL
#LSSL = log(b) + u*log(ct) + u*zt*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.LSSL.S.real[i, ] <- param.S.LSCC$log.b1[i] + 1 * param.S.SAR$log.ct[i]  + 1 * param.S.SAR$zt[i] * x.area.S.real
  pred.LSSL.G.real[i, ] <- param.G.LSCC$log.b1[i] + 1 * param.G.SAR$log.ct[i]  + 1 * param.G.SAR$zt[i] * x.area.G.real
  pred.LSSL.A.real[i, ] <- param.A.LSCC$log.b1[i] + 1 * param.A.SAR$log.ct[i]  + 1 * param.A.SAR$zt[i] * x.area.A.real
}

# make predicted line for each posterior draw
#Eq 10, u=2, CCH
#CCH = log(b) + u*log(ct) + u*zt*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.CCH.S.real[i, ] <- param.S.LSCC$log.b2[i] + 2 * param.S.SAR$log.ct[i]  + 2 * param.S.SAR$zt[i] * x.area.S.real
  pred.CCH.G.real[i, ] <- param.G.LSCC$log.b2[i] + 2 * param.G.SAR$log.ct[i]  + 2 * param.G.SAR$zt[i] * x.area.G.real
  pred.CCH.A.real[i, ] <- param.A.LSCC$log.b2[i] + 2 * param.A.SAR$log.ct[i]  + 2 * param.A.SAR$zt[i] * x.area.A.real
}


#RMSE for IAR
SquaresIAR.S <- mat.or.vec(600,NS)
SumSquaresIAR.S <- mat.or.vec(600,1)
RMSE.IAR.S <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NS){
    SquaresIAR.S[i,j] <- (pred.IAR.S.real[i,j] - seeds.data.int$log.interacts[j])^2
  }
  SumSquaresIAR.S[i] <- sum(SquaresIAR.S[i,])
  RMSE.IAR.S[i] <- sqrt(SumSquaresIAR.S[i]/NS)
}
SquaresIAR.G <- mat.or.vec(600,NG)
SumSquaresIAR.G <- mat.or.vec(600,1)
RMSE.IAR.G <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NG){
    SquaresIAR.G[i,j] <- (pred.IAR.G.real[i,j] - grass.data.int$log.interacts[j])^2
  }
  SumSquaresIAR.G[i] <- sum(SquaresIAR.G[i,])
  RMSE.IAR.G[i] <- sqrt(SumSquaresIAR.G[i]/NG)
}
SquaresIAR.A <- mat.or.vec(600,NAP)
SumSquaresIAR.A <- mat.or.vec(600,1)
RMSE.IAR.A <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NAP){
    SquaresIAR.A[i,j] <- (pred.IAR.A.real[i,j] - ant.data.int$log.interacts[j])^2
  }
  SumSquaresIAR.A[i] <- sum(SquaresIAR.A[i,])
  RMSE.IAR.A[i] <- sqrt(SumSquaresIAR.A[i]/NAP)
}
#RMSE for LSSL
SquaresLSSL.S <- mat.or.vec(600,NS)
SumSquaresLSSL.S <- mat.or.vec(600,1)
RMSE.LSSL.S <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NS){
    SquaresLSSL.S[i,j] <- (pred.LSSL.S.real[i,j] - seeds.data.int$log.interacts[j])^2
  }
  SumSquaresLSSL.S[i] <- sum(SquaresLSSL.S[i,])
  RMSE.LSSL.S[i] <- sqrt(SumSquaresLSSL.S[i]/NS)
}
SquaresLSSL.G <- mat.or.vec(600,NG)
SumSquaresLSSL.G <- mat.or.vec(600,1)
RMSE.LSSL.G <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NG){
    SquaresLSSL.G[i,j] <- (pred.LSSL.G.real[i,j] - grass.data.int$log.interacts[j])^2
  }
  SumSquaresLSSL.G[i] <- sum(SquaresLSSL.G[i,])
  RMSE.LSSL.G[i] <- sqrt(SumSquaresLSSL.G[i]/NG)
}
SquaresLSSL.A <- mat.or.vec(600,NAP)
SumSquaresLSSL.A <- mat.or.vec(600,1)
RMSE.LSSL.A <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NAP){
    SquaresLSSL.A[i,j] <- (pred.LSSL.A.real[i,j] - ant.data.int$log.interacts[j])^2
  }
  SumSquaresLSSL.A[i] <- sum(SquaresLSSL.A[i,])
  RMSE.LSSL.A[i] <- sqrt(SumSquaresLSSL.A[i]/NAP)
}
#RMSE for CCH
SquaresCCH.S <- mat.or.vec(600,NS)
SumSquaresCCH.S <- mat.or.vec(600,1)
RMSE.CCH.S <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NS){
    SquaresCCH.S[i,j] <- (pred.CCH.S.real[i,j] - seeds.data.int$log.interacts[j])^2
  }
  SumSquaresCCH.S[i] <- sum(SquaresCCH.S[i,])
  RMSE.CCH.S[i] <- sqrt(SumSquaresCCH.S[i]/NS)
}
SquaresCCH.G <- mat.or.vec(600,NG)
SumSquaresCCH.G <- mat.or.vec(600,1)
RMSE.CCH.G <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NG){
    SquaresCCH.G[i,j] <- (pred.CCH.G.real[i,j] - grass.data.int$log.interacts[j])^2
  }
  SumSquaresCCH.G[i] <- sum(SquaresCCH.G[i,])
  RMSE.CCH.G[i] <- sqrt(SumSquaresCCH.G[i]/NG)
}
SquaresCCH.A <- mat.or.vec(600,NAP)
SumSquaresCCH.A <- mat.or.vec(600,1)
RMSE.CCH.A <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NAP){
    SquaresCCH.A[i,j] <- (pred.CCH.A.real[i,j] - ant.data.int$log.interacts[j])^2
  }
  SumSquaresCCH.A[i] <- sum(SquaresCCH.A[i,])
  RMSE.CCH.A[i] <- sqrt(SumSquaresCCH.A[i]/NAP)
}

rmse.table <- data.frame(cbind(model = c("IAR", "LSSL", "CCH")), 
                        rbind(mean(RMSE.IAR.S), mean(RMSE.LSSL.S), mean(RMSE.CCH.S)), 
                        rbind(quantile(RMSE.IAR.S, 0.025), quantile(RMSE.LSSL.S, 0.025), quantile(RMSE.CCH.S, 0.025)),
                        rbind(quantile(RMSE.IAR.S, 0.975), quantile(RMSE.LSSL.S, 0.975), quantile(RMSE.CCH.S, 0.975)),
                        rbind(mean(RMSE.IAR.G), mean(RMSE.LSSL.G), mean(RMSE.CCH.G)), 
                        rbind(quantile(RMSE.IAR.G, 0.025), quantile(RMSE.LSSL.G, 0.025), quantile(RMSE.CCH.G, 0.025)),
                        rbind(quantile(RMSE.IAR.G, 0.975), quantile(RMSE.LSSL.G, 0.975), quantile(RMSE.CCH.G, 0.975)),
                        rbind(mean(RMSE.IAR.A), mean(RMSE.LSSL.A), mean(RMSE.CCH.A)), 
                        rbind(quantile(RMSE.IAR.A, 0.025), quantile(RMSE.LSSL.A, 0.025), quantile(RMSE.CCH.A, 0.025)),
                        rbind(quantile(RMSE.IAR.A, 0.975), quantile(RMSE.LSSL.A, 0.975), quantile(RMSE.CCH.A, 0.975)))
colnames(rmse.table) <- c("model", "seeds.mean", "seeds.l95", "seeds.u95", "grass.mean", "grass.l95", "grass.u95", "ant.mean", "ant.l95", "ant.u95")
write.csv(rmse.table, file = paste0("Table1_STOTAL_",Sys.Date(),".csv"))


############################ PREDICTED VS ACTUAL DIFFERENCE ###################################
#take sqrt of Squares to get absolute value of difference between predicted and actual, take the mean across all outputs, take the mean of all data points, un-log the mean
#meanDiff.IAR.S <- 10^mean(apply(sqrt(SquaresIAR.S), 2, mean)) #Note: this is an order of magnitude less than the below (if use squares)


#IAR
Unlog.DiffIAR.S <- mat.or.vec(nr=600, nc=NS)
for(i in 1:600){
  for(j in 1:NS){
Unlog.DiffIAR.S[i,j] <- 10^(pred.IAR.S.real[i,j]) - 10^(seeds.data.int$log.interacts[j])
  }
}
ULmeanDiff.IAR.S <- mean(apply(Unlog.DiffIAR.S, 2, mean)) #-3.67

Unlog.DiffIAR.G <- mat.or.vec(nr=600, nc=NG)
for(i in 1:600){
  for(j in 1:NG){
    Unlog.DiffIAR.G[i,j] <- 10^(pred.IAR.G.real[i,j]) - 10^(grass.data.int$log.interacts[j])
  }
}
ULmeanDiff.IAR.G <- mean(apply(Unlog.DiffIAR.G, 2, mean)) #-0.40

Unlog.DiffIAR.A <- mat.or.vec(nr=600, nc=NAP)
for(i in 1:600){
  for(j in 1:NAP){
    Unlog.DiffIAR.A[i,j] <- 10^(pred.IAR.A.real[i,j]) - 10^(ant.data.int$log.interacts[j])
  }
}
ULmeanDiff.IAR.A <- mean(apply(Unlog.DiffIAR.A, 2, mean)) #4.49

#LSSL
Unlog.DiffLSSL.S <- mat.or.vec(nr=600, nc=NS)
for(i in 1:600){
  for(j in 1:NS){
    Unlog.DiffLSSL.S[i,j] <- 10^(pred.LSSL.S.real[i,j]) - 10^(seeds.data.int$log.interacts[j])
  }
}
ULmeanDiff.LSSL.S <- mean(apply(Unlog.DiffLSSL.S, 2, mean)) #7.23

Unlog.DiffLSSL.G <- mat.or.vec(nr=600, nc=NG)
for(i in 1:600){
  for(j in 1:NG){
    Unlog.DiffLSSL.G[i,j] <- 10^(pred.LSSL.G.real[i,j]) - 10^(grass.data.int$log.interacts[j])
  }
}
ULmeanDiff.LSSL.G <- mean(apply(Unlog.DiffLSSL.G, 2, mean)) #34.1

Unlog.DiffLSSL.A <- mat.or.vec(nr=600, nc=NAP)
for(i in 1:600){
  for(j in 1:NAP){
    Unlog.DiffLSSL.A[i,j] <- 10^(pred.LSSL.A.real[i,j]) - 10^(ant.data.int$log.interacts[j])
  }
}
ULmeanDiff.LSSL.A <- mean(apply(Unlog.DiffLSSL.A, 2, mean)) #4.31

#CCH
Unlog.DiffCCH.S <- mat.or.vec(nr=600, nc=NS)
for(i in 1:600){
  for(j in 1:NS){
    Unlog.DiffCCH.S[i,j] <- 10^(pred.CCH.S.real[i,j]) - 10^(seeds.data.int$log.interacts[j])
  }
}
ULmeanDiff.CCH.S <- mean(apply(Unlog.DiffCCH.S, 2, mean)) #46.8

Unlog.DiffCCH.G <- mat.or.vec(nr=600, nc=NG)
for(i in 1:600){
  for(j in 1:NG){
    Unlog.DiffCCH.G[i,j] <- 10^(pred.CCH.G.real[i,j]) - 10^(grass.data.int$log.interacts[j])
  }
}
ULmeanDiff.CCH.G <- mean(apply(Unlog.DiffCCH.G, 2, mean)) #226

Unlog.DiffCCH.A <- mat.or.vec(nr=600, nc=NAP)
for(i in 1:600){
  for(j in 1:NAP){
    Unlog.DiffCCH.A[i,j] <- 10^(pred.CCH.A.real[i,j]) - 10^(ant.data.int$log.interacts[j])
  }
}
ULmeanDiff.CCH.A <- mean(apply(Unlog.DiffCCH.A, 2, mean)) #316


