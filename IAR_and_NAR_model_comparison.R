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
Sugiara <- subset(AllNetworks, AllNetworks$Dataset=="Sugiara")


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
  z1 ~ dnorm(0, 0.01) 
  z2 ~ dnorm(0, 0.01) 
  # u is given as data as assumption of both models
  sigma1 ~ dunif(0, 100)
  tau1 <- 1/(sigma1*sigma1)
  sigma2 ~ dunif(0, 100)
  tau2 <- 1/(sigma2*sigma2)

  for(i in 1:n.points){
    mu1[i] <- log.c1 + z1 * log.Area[i] #SAR for spp group 1
    log.Spp1[i] ~ dnorm(mu1[i], tau1)
    mu2[i] <- log.c2 + z2 * log.Area[i] #SAR for spp group 2
    log.Spp2[i] ~ dnorm(mu2[i], tau2)
  }
}

# A third model that calculates link-scaling and constant connectance b parameter
LSCCs <- function(){
  log.b1 ~ dnorm(0, 0.01) 
  log.b2 ~ dnorm(0, 0.01) 
  sigma1 ~ dunif(0, 100)
  tau1 <- 1/(sigma1*sigma1)
  sigma2 ~ dunif(0, 100)
  tau2 <- 1/(sigma2*sigma2)
  
  for(i in 1:n.points){    
    mu1[i] <- log.b1 + 1 * log.Spp2[i] #Link-species scaling relationship
    log.interacts.11.1[i] ~ dnorm(mu1[i], tau1)
    mu2[i] <- log.b2 + 2 * log.Spp2[i] #Link-species scaling relationship
    log.interacts.11.2[i] ~ dnorm(mu2[i], tau2)
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
  n.points = dim(Sugiara)[1],
  log.connectance = log10(Sugiara$Connectance/100),
  log.Area = log10(Sugiara$Patch_Area_m2)
)

CAR.S <- jags(data=seeds.data.car, inits = list(list(alpha=1, delta=0), list(alpha=-1, delta=1), list(alpha=3, delta=-1)), parameters.to.save= parameters.CAR, CARs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
CAR.G <- jags(data=grass.data.car, inits = list(list(alpha=1, delta=0), list(alpha=-1, delta=1), list(alpha=3, delta=-1)), parameters.to.save= parameters.CAR, CARs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
CAR.A <- jags(data=ant.data.car, inits = list(list(alpha=1, delta=0), list(alpha=1, delta=1), list(alpha=3, delta=-1)), parameters.to.save= parameters.CAR, CARs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
# Remember to check fits, Rhat, N.eff

#################### Fit SAR models to all 3 datasets #################### 
parameters.SAR <- c("log.c1", "log.c2", "z1", "z2") 

seeds.data.sar <- list(
             n.points = dim(Seeds)[1],
             log.Spp1 = log10(Seeds$Lower_Spp_Rich),
             log.Spp2 = log10(Seeds$Upper_Spp_Rich),
             log.Area = log10(Seeds$Patch_Area_m2)
             )
grass.data.sar <- list(
                   n.points = dim(Grass1)[1],
                   log.Spp1 = log10(Grass1$Lower_Spp_Rich),
                   log.Spp2 = log10(Grass1$Upper_Spp_Rich),
                   log.Area = log10(Grass1$Patch_Area_m2)
)
ant.data.sar <- list(
                 n.points = dim(Sugiara)[1],
                 log.Spp1 = log10(Sugiara$Lower_Spp_Rich),
                 log.Spp2 = log10(Sugiara$Upper_Spp_Rich),
                 log.Area = log10(Sugiara$Patch_Area_m2)
)

SAR.S <- jags(data=seeds.data.sar, inits = list(list(log.c1=3, log.c2=3, z1=0.3, z2=0.35), list(log.c1=1, log.c2=1, z1=0.2, z2=0.25), list(log.c1=-1, log.c2=-1, z1=0.1, z2=0.15)), parameters.to.save= parameters.SAR, SARs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
SAR.G <- jags(data=grass.data.sar, inits = list(list(log.c1=3, log.c2=3, z1=0.3, z2=0.35), list(log.c1=1, log.c2=1, z1=0.2, z2=0.25), list(log.c1=-1, log.c2=-1, z1=0.1, z2=0.15)), parameters.to.save= parameters.SAR, SARs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
SAR.A <- jags(data=ant.data.sar, inits = list(list(log.c1=3, log.c2=3, z1=0.3, z2=0.35), list(log.c1=1, log.c2=1, z1=0.2, z2=0.25), list(log.c1=-1, log.c2=-1, z1=0.1, z2=0.15)), parameters.to.save= parameters.SAR, SARs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
# Remember to check fits, Rhat, N.eff

#################### Fit LSS & CCH models to all 3 datasets #################### 
parameters.LSCC <- c("log.b1", "log.b2") 

seeds.data.lscc <- list(
  n.points = dim(Seeds)[1],
  log.interacts.11.1 = log10(Seeds$No_Int),
  log.interacts.11.2 = log10(Seeds$No_Int),
  log.Spp2 = log10(Seeds$Upper_Spp_Rich)
)
grass.data.lscc <- list(
  n.points = dim(Grass1)[1],
  log.interacts.11.1 = log10(Grass1$No_Int),
  log.interacts.11.2 = log10(Grass1$No_Int),  
  log.Spp2 = log10(Grass1$Upper_Spp_Rich)
)
ant.data.lscc <- list(
  n.points = dim(Sugiara)[1],
  log.interacts.11.1 = log10(Sugiara$No_Int),
  log.interacts.11.2 = log10(Sugiara$No_Int),  
  log.Spp2 = log10(Sugiara$Upper_Spp_Rich)
)

LSCC.S <- jags(data=seeds.data.lscc, inits = list(list(log.b1=0.3, log.b2=-0.5), list(log.b1=0.1, log.b2=-0.01), list(log.b1=0.5, log.b2=-1.1)), parameters.to.save= parameters.LSCC, LSCCs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
LSCC.G <- jags(data=grass.data.lscc, inits = list(list(log.b1=0.3, log.b2=-0.5), list(log.b1=0.1, log.b2=-0.01), list(log.b1=0.5, log.b2=-1.1)), parameters.to.save= parameters.LSCC, LSCCs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
LSCC.A <- jags(data=ant.data.lscc, inits = list(list(log.b1=0.3, log.b2=-0.5), list(log.b1=0.1, log.b2=-0.01), list(log.b1=0.5, log.b2=-1.1)), parameters.to.save= parameters.LSCC, LSCCs, n.burnin=230000, n.iter=250000, n.thin=100, n.chains=3)
# Remember to check fits, Rhat, N.eff


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

# extract parameters
param.S.LSCC <- MCMCpstr(LSCC.S, params = parameters.LSCC, type = "chains")
param.G.LSCC <- MCMCpstr(LSCC.G, params = parameters.LSCC, type = "chains")
param.A.LSCC <- MCMCpstr(LSCC.A, params = parameters.LSCC, type = "chains")

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
#Eq 11, u=1, LSSL
#LSSL = log(b) + u*z2*log(c2) + u*z2*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.LSSL.S[i, ] <- param.S.LSCC$log.b1[i] + 1 * param.S.SAR$z2[i] * param.S.SAR$log.c2[i]  + 1 * param.S.SAR$z2[i] * x.area.S
  pred.LSSL.G[i, ] <- param.G.LSCC$log.b1[i] + 1 * param.G.SAR$z2[i] * param.G.SAR$log.c2[i]  + 1 * param.G.SAR$z2[i] * x.area.G
  pred.LSSL.A[i, ] <- param.A.LSCC$log.b1[i] + 1 * param.A.SAR$z2[i] * param.A.SAR$log.c2[i]  + 1 * param.A.SAR$z2[i] * x.area.A
}

# make predicted line for each posterior draw
#Eq 11, u=2, CCH
#CCH = log(b) + u*z2*log(c2) + u*z2*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.CCH.S[i, ] <- param.S.LSCC$log.b2[i] + 2 * param.S.SAR$z2[i] * param.S.SAR$log.c2[i]  + 2 * param.S.SAR$z2[i] * x.area.S
  pred.CCH.G[i, ] <- param.G.LSCC$log.b2[i] + 2 * param.G.SAR$z2[i] * param.G.SAR$log.c2[i]  + 2 * param.G.SAR$z2[i] * x.area.G
  pred.CCH.A[i, ] <- param.A.LSCC$log.b2[i] + 2 * param.A.SAR$z2[i] * param.A.SAR$log.c2[i]  + 2 * param.A.SAR$z2[i] * x.area.A
}


#################### Plot Figure 2 #################### 

pdf(file = paste0("Figure2_",Sys.Date(),".pdf"), width = 3.25, height = 6)
par(las = 1, mar = c(2, 3, 0.5, 0.5), oma=c(1,0,0,0),
    tcl = -0.2, mgp = c(1.8,0.25,0), mfrow = c(3, 1))

plot(seeds.data.sar$log.Area, seeds.data.lscc$log.interacts.11.1, ylim = c(0.3, 3.2),
     xlab = "", ylab = "log(Interactions)", 
     pch=21, col="black", bg="#00000080", cex = 1.5)
lines(x.area.S, apply(pred.LSSL.S, 2, mean), col = "#929084", lwd = 2, lty = 4)
lines(x.area.S, apply(pred.CCH.S, 2, mean), col = "#2e4052", lwd = 2, lty = 2)
lines(x.area.S, apply(pred.IAR.S, 2, mean), col = "#e5323b", lwd = 2, lty = 1)
title(main = "A. Seed dispersal", adj = 0.01, line = -1.2)

plot(grass.data.sar$log.Area, grass.data.lscc$log.interacts.11.1, ylim = c(0, 2.5),
     xlab = "", ylab = "log(Interactions)", 
     pch=22, col="black", bg="#00000080", cex = 1.5)
lines(x.area.G, apply(pred.LSSL.G, 2, mean), col = "#929084",  lwd = 2, lty = 4)
lines(x.area.G, apply(pred.CCH.G, 2, mean), col = "#2e4052", lwd = 2, lty = 2)
lines(x.area.G, apply(pred.IAR.G, 2, mean), col = "#e5323b", lwd = 2, lty = 1)
title(main = "B. Pollination", adj = 0.01, line = -1.2)

plot(ant.data.sar$log.Area, ant.data.lscc$log.interacts.11.1, ylim = c(0.9, 1.9),
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



######################### CALCULATE SUM OF SQUARES ###############################
NS <- length(seeds.data.lscc$log.interacts.11.1)
NG <- length(grass.data.lscc$log.interacts.11.1)
NAP <- length(ant.data.lscc$log.interacts.11.1)

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
#Eq 11, u=1, LSSL
#LSSL = log(b) + u*z2*log(c2) + u*z2*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.LSSL.S.real[i, ] <- param.S.LSCC$log.b1[i] + 1 * param.S.SAR$z2[i] * param.S.SAR$log.c2[i]  + 1 * param.S.SAR$z2[i] * x.area.S.real
  pred.LSSL.G.real[i, ] <- param.G.LSCC$log.b1[i] + 1 * param.G.SAR$z2[i] * param.G.SAR$log.c2[i]  + 1 * param.G.SAR$z2[i] * x.area.G.real
  pred.LSSL.A.real[i, ] <- param.A.LSCC$log.b1[i] + 1 * param.A.SAR$z2[i] * param.A.SAR$log.c2[i]  + 1 * param.A.SAR$z2[i] * x.area.A.real
}

# make predicted line for each posterior draw
#Eq 11, u=2, CCH
#CCH = log(b) + u*z2*log(c2) + u*z2*log(A)
for(i in 1:length(param.S.SAR$log.c2)) {
  pred.CCH.S.real[i, ] <- param.S.LSCC$log.b2[i] + 2 * param.S.SAR$z2[i] * param.S.SAR$log.c2[i]  + 2 * param.S.SAR$z2[i] * x.area.S.real
  pred.CCH.G.real[i, ] <- param.G.LSCC$log.b2[i] + 2 * param.G.SAR$z2[i] * param.G.SAR$log.c2[i]  + 2 * param.G.SAR$z2[i] * x.area.G.real
  pred.CCH.A.real[i, ] <- param.A.LSCC$log.b2[i] + 2 * param.A.SAR$z2[i] * param.A.SAR$log.c2[i]  + 2 * param.A.SAR$z2[i] * x.area.A.real
}


#SSE for IAR
SquaresIAR.S <- mat.or.vec(600,NS)
SumSquaresIAR.S <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NS){
    SquaresIAR.S[i,j] <- (pred.IAR.S.real[i,j] - seeds.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresIAR.S[i] <- sum(SquaresIAR.S[i,])
}
SquaresIAR.G <- mat.or.vec(600,NG)
SumSquaresIAR.G <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NG){
    SquaresIAR.G[i,j] <- (pred.IAR.G.real[i,j] - grass.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresIAR.G[i] <- sum(SquaresIAR.G[i,])
}
SquaresIAR.A <- mat.or.vec(600,NAP)
SumSquaresIAR.A <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NAP){
    SquaresIAR.A[i,j] <- (pred.IAR.A.real[i,j] - ant.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresIAR.A[i] <- sum(SquaresIAR.A[i,])
}
#SSE for LSSL
SquaresLSSL.S <- mat.or.vec(600,NS)
SumSquaresLSSL.S <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NS){
    SquaresLSSL.S[i,j] <- (pred.LSSL.S.real[i,j] - seeds.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresLSSL.S[i] <- sum(SquaresLSSL.S[i,])
}
SquaresLSSL.G <- mat.or.vec(600,NG)
SumSquaresLSSL.G <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NG){
    SquaresLSSL.G[i,j] <- (pred.LSSL.G.real[i,j] - grass.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresLSSL.G[i] <- sum(SquaresLSSL.G[i,])
}
SquaresLSSL.A <- mat.or.vec(600,NAP)
SumSquaresLSSL.A <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NAP){
    SquaresLSSL.A[i,j] <- (pred.LSSL.A.real[i,j] - ant.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresLSSL.A[i] <- sum(SquaresLSSL.A[i,])
}
#SSE for CCH
SquaresCCH.S <- mat.or.vec(600,NS)
SumSquaresCCH.S <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NS){
    SquaresCCH.S[i,j] <- (pred.CCH.S.real[i,j] - seeds.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresCCH.S[i] <- sum(SquaresCCH.S[i,])
}
SquaresCCH.G <- mat.or.vec(600,NG)
SumSquaresCCH.G <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NG){
    SquaresCCH.G[i,j] <- (pred.CCH.G.real[i,j] - grass.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresCCH.G[i] <- sum(SquaresCCH.G[i,])
}
SquaresCCH.A <- mat.or.vec(600,NAP)
SumSquaresCCH.A <- mat.or.vec(600,1)
for(i in 1:600){
  for(j in 1:NAP){
    SquaresCCH.A[i,j] <- (pred.CCH.A.real[i,j] - ant.data.lscc$log.interacts.11.1[j])^2
  }
  SumSquaresCCH.A[i] <- sum(SquaresCCH.A[i,])
}

sse.table <- data.frame(cbind(model = c("IAR", "LSSL", "CCH")), 
                        rbind(mean(SumSquaresIAR.S), mean(SumSquaresLSSL.S), mean(SumSquaresCCH.S)), 
                        rbind(quantile(SumSquaresIAR.S, 0.025), quantile(SumSquaresLSSL.S, 0.025), quantile(SumSquaresCCH.S, 0.025)),
                        rbind(quantile(SumSquaresIAR.S, 0.975), quantile(SumSquaresLSSL.S, 0.975), quantile(SumSquaresCCH.S, 0.975)),
                        rbind(mean(SumSquaresIAR.G), mean(SumSquaresLSSL.G), mean(SumSquaresCCH.G)), 
                        rbind(quantile(SumSquaresIAR.G, 0.025), quantile(SumSquaresLSSL.G, 0.025), quantile(SumSquaresCCH.G, 0.025)),
                        rbind(quantile(SumSquaresIAR.G, 0.975), quantile(SumSquaresLSSL.G, 0.975), quantile(SumSquaresCCH.G, 0.975)),
                        rbind(mean(SumSquaresIAR.A), mean(SumSquaresLSSL.A), mean(SumSquaresCCH.A)), 
                        rbind(quantile(SumSquaresIAR.A, 0.025), quantile(SumSquaresLSSL.A, 0.025), quantile(SumSquaresCCH.A, 0.025)),
                        rbind(quantile(SumSquaresIAR.A, 0.975), quantile(SumSquaresLSSL.A, 0.975), quantile(SumSquaresCCH.A, 0.975)))
colnames(sse.table) <- c("model", "seeds.mean", "seeds.l95", "seeds.u95", "grass.mean", "grass.l95", "grass.u95", "ant.mean", "ant.l95", "ant.u95")
write.csv(sse.table, file = paste0("Table1_",Sys.Date(),".csv"))

############################ NOT CORRECT FOR ABOVE CODE AND MAY NOT NEED THIS ##############################################
### Calculate R^2
pred.IA.S.R <- array(dim = c(length(param.S$log.c1), length(Seeds$Patch_Area_m2)))
pred.IA.G.R <- array(dim = c(length(param.S$log.c1), length(Grass1$Patch_Area_m2)))
pred.IA.A.R <- array(dim = c(length(param.S$log.c1), length(Sugiara$Patch_Area_m2)))
r2.S <- mat.or.vec(nr=length(param.S$log.c2), nc=1)
r2.G <- mat.or.vec(nr=length(param.S$log.c2), nc=1)
r2.A <- mat.or.vec(nr=length(param.S$log.c2), nc=1)

for(i in 1:length(param.S$log.c2)) {
  pred.IA.S.R[i, ] <- param.S$alpha[i] + param.S$log.c1[i] + param.S$log.c2[i] + (param.S$delta[i] + param.S$z1[i] + param.S$z2[i]) * log10(Seeds$Patch_Area_m2)
  r2.S[i] <- (cor(log10(Seeds$No_Int), pred.IA.S.R[i,]))^2
  
  pred.IA.G.R[i, ] <- param.G$alpha[i] + param.G$log.c1[i] + param.G$log.c2[i] + (param.G$delta[i] + param.G$z1[i] + param.G$z2[i]) * log10(Grass1$Patch_Area_m2)
  r2.G[i] <- (cor(log10(Grass1$No_Int), pred.IA.G.R[i,]))^2  
  
  pred.IA.A.R[i, ] <- param.A$alpha[i] + param.A$log.c1[i] + param.A$log.c2[i] + (param.A$delta[i] + param.A$z1[i] + param.A$z2[i]) * log10(Sugiara$Patch_Area_m2)
  r2.A[i] <- (cor(log10(Sugiara$No_Int), pred.IA.A.R[i,]))^2  
}
