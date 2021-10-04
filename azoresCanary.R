##
# THIS IS A DEPCRECATED VERSION - PLEASE USE VERSION 1.02
#
# Bayesian regression analyses of the number of species in Azores and Canary Islands
# using STAN through the package rethinking, by Mc Elreath
#
# Author: Diogo Barros
# Version: 1.01
##

source("/home/splinter/Documents/BEAG/R scripts/Task1/myFunctions2.R")
setwd("/home/splinter/Documents/R data/Paper1")
library(rethinking)
library(readxl)
library(dplyr)
options(mC.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


###
#Preparing the data
###

#AZORES
azoresFull <- as.data.frame(read_excel("azores_arthopoda.xlsx", sheet=1))
azoresVars <- as.data.frame(read_excel("MAC_EnvData.xlsx", sheet=1, col_names = FALSE))
azoresEnd <- filter(azoresFull, COL.=='END')
azores <- azoresVars[9:17,]
azoresNames = azores[,1]
azores <- df2matrix(azores)
azores <- as.data.frame(azores)
azores$V1 <- azoresNames
colnames(azores) <- azoresVars[7,]
azores$index <- 1:9 
#Getting species richness
azoresData<-azoresFull[20:length(azoresFull)]
azoresEndData<-azoresEnd[20:length(azoresEnd)]
azoresSIEData<- filter(azoresEndData, rowSums(azoresEndData, na.rm = TRUE) == 1)
azores$richness = colSums(azoresData, na.rm = TRUE)/2 #/2 removes total from excel
azores$endmisms = colSums(azoresEndData, na.rm = TRUE)
azores$SIE = colSums(azoresSIEData, na.rm = TRUE)
#Distance matrix
dmAzores = distMatrix(df2matrix(azores[2:3]))
dmAzores=dmAzores/1000  #distance in km
dmAzores=dmAzores/100  #distance in hundreds of Km
colnames(dmAzores) <- colnames(azoresData)
round(dmAzores,3)

#CANARY ISLANDS
canaryFull <- as.data.frame(read_excel("ArthropodaCanaryIslands.xls", sheet=1))
canaryVars <- as.data.frame(read_excel("MAC_EnvData.xlsx", sheet=1, col_names = FALSE))
canaryEnd=filter(canaryFull, END=='End')
canary <- canaryVars[36:46,]
canaryNames = canary[,1]
canary <- df2matrix(canary)
canary <- as.data.frame(canary)
canary$V1 <- canaryNames
colnames(canary) <- canaryVars[7,]
canary=canary[c(1:7),]
canary$index <- 1:7
#Getting species richness
canaryData<-canaryFull[18:length(canaryFull)]
canaryData=data.frame(canaryData$fuerteventura, canaryData$lanzarote,canaryData$canaria,canaryData$tenerife,
                      canaryData$gomera,canaryData$hierro,canaryData$palma)
canaryEndData <-canaryEnd[18:length(canaryEnd)] 
canaryEndData=data.frame(canaryEnd$fuerteventura, canaryEnd$lanzarote,canaryEnd$canaria,canaryEnd$tenerife,
                      canaryEnd$gomera,canaryEnd$hierro,canaryEnd$palma)
canarySIEData<- filter(canaryEndData, rowSums(canaryEndData, na.rm = TRUE) == 1)
canary$richness = colSums(canaryData, na.rm = TRUE)/2 #/2 removes total from excel
canary$endemisms = colSums(canaryEndData, na.rm = TRUE)
canary$SIE = colSums(canarySIEData, na.rm = TRUE)
#Distance matrix
dmCanary = distMatrix(df2matrix(canary[2:3]))
dmCanary=dmCanary/1000  #distance in km
dmCanary=dmCanary/100  #distance in hundreds of Km
colnames(dmCanary) <- colnames(canaryData)
round(dmCanary,3)

###########
#Modeling Species Richness under a Bayesian approach
###########
island=azores
dmisland=dmAzores
#island=canary
#dmisland=dmCanary
richness=island$endmisms

###
#Model 1: spp~area+ distance
###


### log variables
mC.1 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + bp*area,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=log(island$`A (km2)`),
    index=island$index,
    Dmat=dmisland),
  warmup=4000 ,cores = 4, iter=1e4 , chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

### log area + non fixed sigma
mC.1.3 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + bp*area,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=log(island$`A (km2)`),
    index=island$index,
    Dmat=dmisland),
  warmup=4000 ,cores=4, iter=1e4 , chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

### sdt area variables
mC.1.4 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + bp*area,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=strIsland$`A (km2)`,
    index=island$index,
    Dmat=dmisland),
  warmup=2000 , iter=1e4 ,cores = 4, chains=4)


###
#Model 2: spp ~ area + distance + age
####
mC.2 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*age,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=logIsland$`A (km2)`,
    age=island$`Age (Ma)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 , chains=4,cores=4,control=list(max_treedepth =15,adapt_delta=0.9))

##log age: CRAZY SLOW DOES NOT SCALE
#mC.2.2 <- map2stan(
#  alist(
#    richness ~ dpois(lambda),
#    log(lambda) <- a + g[index] + b1*area + b2*age,
#    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01),
#    a ~ dnorm(0,10),
#    b1 ~ dnorm(0,1),
#    b2 ~ dnorm(0,1),
#    etasq ~ dcauchy(0,1),
#    rhosq ~ dcauchy(0,1),
#    sigma ~ dcauchy(0,1)
#  ),
#  data=list(
#    richness=richness,
#    area=island$`A (km2)`,
#    age =logIsland$`Age (Ma)`,
#    index=island$index,
#    Dmat=dmisland),
#  warmup=4000 , iter=1e4 , chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

##standerized everything
mC.2.3 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*age,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=strIsland$`A (km2)`,
    age =strIsland$`Age (Ma)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))


### 
#Model 3: spp ~ area + distance + altitude
###
init_fun <- function() { list(a=sqrt(rnorm(1,0,10)^2),
                              b1=sqrt(rnorm(1,0,1)^2),
                              b2=sqrt(rnorm(1,0,1)^2),
                              etasq=sqrt(rcauchy(1,0,1)^2),
                              rhosq=sqrt(rcauchy(1,0,1)^2),
                              sigma=sqrt(rcauchy(1,0,1)^2))} 
mC.3 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*altitude,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=log((island$`A (km2)`)),
    altitude=(island$`ELEV (m)`),
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 , chains=4,cores=4, init = init_fun, control=list(max_treedepth =15,adapt_delta=0.9))

### standerized Altitude

mA.3.2 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*altitude,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=logIsland$`A (km2)`,
    altitude = ((island$`ELEV (m)`)-mean((island$`ELEV (m)`)))/sd((island$`ELEV (m)`)),
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 , cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

###
#Model 4: spp ~ area + distance + age + altitude 
###

mC.4 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*altitude + b3*age,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=logIsland$`A (km2)`,
    altitude=island$`ELEV (m)`,
    age=island$`Age (Ma)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))


## standerized elevation and age
mC.4.2 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*altitude + b3*age,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=logIsland$`A (km2)`,
    altitude = strIsland$`ELEV (m)`,
    age = strIsland$`Age (Ma)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

###
#Model 5: spp ~ area + distance + distance to mainland 
#A complete biogeographic model
###

mC.5 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*d2MLand,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=logIsland$`A (km2)`,
    d2MLand = island$`DM (km)`/100,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))


###Standerized Distance
mC.5.2 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*d2MLand,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area= island$`A (km2)`,
    d2MLand = strIsland$`DM (km)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

###
#Model 6: spp ~ area + distance + distance to mainland + age + altitude 
#A complete biogeographic model + 2 extra variables
###

mC.6 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area + b2*d2MLand + b3*age + b4*altitude,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~ dnorm(0,1),
    b4 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=logIsland$`A (km2)`,
    d2MLand = island$`DM (km)`/100,
    age=island$index,
    altitude=island$index,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))


###
#Model 7: simple Spp~area regression
###

mC.SCAZ <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a + bp*area,
    a ~ dnorm(0,100),
    bp ~ dnorm(0,1)
  ),
  data=list(
    richness=richness,
    area=logIsland$`A (km2)`),
  warmup=2000 , iter=1e4 ,cores=4, chains=4)

###
#Model 8: null model
###

mC.Null <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a, 
    a ~ dnorm(0,100)
  ),
  data=list(
    richness=richness),
  warmup=2000 , iter=1e4 ,cores=4, chains=4)


###
#Model 9: Non gaussian processess, do not include spatial autocorrelation
###
mC.NGP1 <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a + b1*area+b2*altitude+b3*age+b4*d2MLand,
    a ~ dnorm(0,100),
    c(b1,b2,b3,b4) ~ dnorm(0,1)
  ),
  data=list(
    richness=richness,
    area=log(island$`A (km2)`),
    age=island$`Age (Ma)`,
    altitude = island$`ELEV (m)`,
    d2MLand = island$`DM (km)`/100),
  warmup=2000 , iter=1e4 ,cores=4, chains=4)

##
#Model 10: spatial auto correlation only
##
mC.10 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index],
    a ~ dnorm(0,100),
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 , cores=4,chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

##
#Model 11: Multiplicative
##

dat_list <- list(
  T = as.numeric(island$richness),
  P = as.numeric(island$`A (km2)`),
  index = island$index,
  Dmat=dmisland )
mC.11 <- ulam(
  alist(
    T ~ dpois(lambda),
    lambda <- (a*P^b/g)*exp(k[index]),
    vector[10]:k ~ multi_normal( 0 , SIGMA ),
    matrix[10,10]:SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
    c(a,b,g) ~ dexp( 1 ),
    etasq ~ dexp( 2 ),
    rhosq ~ dexp( 0.5 )
  ), data=dat_list , chains=4 , cores=4 , iter=2000 )




##
#OUTPUT

#write.table(compare(mC.1,mC.1.2,mC.1.2.teste,mC.1.3)@output,"waics2.csv", col.names = FALSE, append= TRUE)
a=compare(mC.1,mC.1.3,mC.1.4,mC.10,mC.2,mC.2.3,mC.3.2,mC.4,mC.5,mC.4.2,mC.5.2,mC.6,mC.NGP1,mC.Null,mC.SCAZ)
a
write.table(a@output,"waicsCanaryEnd.csv")


###
#SAVING DAS MODELENS

#saveRDS()


#######
# Further down the rabbit hole: covariance analyses
#######
# compute posterior median covariance among islands
post = extract.samples(mA.1)
#post = extract.samples(mA.3.2)
dm=dmisland
K <- matrix(0,nrow=length(dmisland[1,]),ncol=length(dmisland[1,]))
for ( i in 1:length(dmisland[1,]) )
  for ( j in 1:length(dmisland[1,]) )
    K[i,j] <- median(post$etasq) *
  exp( -median(post$rhosq) * dm[i,j]^2)

diag(K) <- median(post$etasq) + 0.01

# convert to correlation matrix
Rho <- round( cov2cor(K) , 2 )
# add row/col names for convenience
colnames(Rho) <- colnames(dm)
rownames(Rho) <- colnames(Rho)
Rho

# scale point size to logpop
psize <- island$`A (km2)` / max(island$`A (km2)`)*4
#psize <- exp(psize*1.5)-2
# plot raw data and labels
plot( island$Longitude , island$Latitude , xlab="longitude" , ylab="latitude",
      col=rangi2 , cex=psize, pch=16 , xlim=c(min(island$Longitude)-1,max(island$Longitude)+1))
labels <- as.character(island$`Archipelago/Island`)
text( island$Longitude , island$Latitude , labels=labels , cex=0.7 , pos=c(1,1,1,3,1,3,1,3,3) )
# overlay lines shaded by Rho
for( i in 1:length(dmisland[1,]))
  for ( j in 1:length(dmisland[1,]))
    if ( i < j )
      lines( c( island$Longitude[i],island$Longitude[j] ) ,
             c( island$Latitude[i],island$Latitude[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2 ))


# compute posterior median relationship, ignoring distance
area.seq <- seq( from=4, to=9 , length.out=30 )
lambda <- sapply( area.seq , function(lp) exp( post$a + post$b1*lp))# + post$b2*lp +post$b3*lp))
lambda.median <- apply( lambda , 2 , median )
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.9 )
# plot raw data and labels
  plot( log(island$`A (km2)`) , log(island$richness) , col=rangi2 , cex=psize , pch=16 ,
      xlab="log area" , ylab="species richness",
      xlim=c(min(log(island$`A (km2)`))-1,max(log(island$`A (km2)`))+1),
      ylim=c(min(log(island$richness))-1,max(log(island$richness))+1)
      )
text( log(island$`A (km2)`) , log(island$richness) , labels=labels , cex=0.7 ,
      pos=c(3,3,4,1,2,1,1,3,3,2) )

# display posterior predictions
lines( area.seq , log(lambda.median), lty=2 )
lines( area.seq , log(lambda.PI80[1,]), lty=2 )
lines( area.seq , log(lambda.PI80[2,]), lty=2 )
# overlay correlations
for( i in 1:7 )
  for ( j in 1:7 )
    if ( i < j )
      lines( c( log(island$`A (km2)`[i]),log(island$`A (km2)`[j]) ) ,
             c( log(island$richness[i]),log(island$richness[j]) ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )



# plot the posterior median covariance function
curve( median(post$etasq)*exp(-median(post$rhosq)*x^2) , from=0 , to=10 ,
       xlab="distance (thousand km)" , ylab="covariance" , ylim=c(0,1) ,
       yaxp=c(0,1,4) , lwd=2 )
# plot 100 functions sampled from posterior
for ( i in 1:100 )
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE ,
         col=col.alpha("black",0.2) )


##
#Posterior plots
##

post2=as.matrix(mC.4.2@stanfit)

canary_title1 <- ggtitle("Posterior distributions for Azores model",
                      "with medians and 95% intervals")
mcmc_areas(post2,
           pars = c("a", "b1", "b2", "b3",
                    "g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]", "g[7]")
                    #"etasq","rhosq"),
           ,prob = 0.95) + canary_title1


color_scheme_set("red")
ppc_dens_overlay(y = post$a,
                 yrep = posterior_predict(mC.4.2@stanfit, draws = 50))
azores
plot(log(azores$richness)~log(azores$`A (km2)`),pch=19)
  abline(log(azores$richness)~log(azores$`A (km2)`), col="red", lty=1)

plot(log(canary$richness)~log(canary$`A (km2)`),pch=19)
abline(log(canary$richness), log(canary$`A (km2)`), col="red", lty=1)

#Other plots
plot(island$richness~island$`A (km2)`)


#ggpair
#ggally
#BAYESPLOT


library("bayesplot")
library("rstanarm")
library("ggplot2")

post <- extract.samples(m14.7)
post2=as.matrix(m14.7@stanfit)
title <- ggtitle("Posterior distributions for Azores - 1st Edition")
mcmc_areas(post2,
          pars = colnames(post2)[-length(post2[1,])],
                  #"a", "b", #"b2", "b3",
                  #"k[1]", "k[2]", "k[3]", "k[4]", "k[5]", "k[6]", "k[7]",
                  #etasq","rhosq"),
           prob = 0.95) + title

-c(colnames(length(post2[1,])))

cDf=data.frame(canary$`A (km2)`,canary$`DM (km)`,canary$`Age (Ma)`,canary$`ELEV (m)`)
pairs(cDf)




