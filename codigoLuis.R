###
# Modelling species richness using Bayesina Mehtods for Azores and Canary Islands
#
# NOTE: PLEASE INSERT WORKING DIRECTORY BEFORE RUNING THE SCRIPT
##
source("/home/splinter/Documents/BEAG/R scripts/Task1/myFunctions2.R")
setwd("/home/splinter/Documents/R data/Paper1")
library(rethinking)
library(readxl)
library(dplyr)
options(mC.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Preparing the data:

########
#AZORES

azoresFull <- as.data.frame(read_excel("azores_arthopoda.xlsx", sheet=1))
azoresVars <- as.data.frame(read_excel("MAC_EnvData.xlsx", sheet=1, col_names = FALSE))
azoresEnd=filter(azoresFull, COL.=='END')
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
azores$richness = colSums(azoresData, na.rm = TRUE)/2 #/2 removes total from excel
azores$endmisms = colSums(azoresEndData, na.rm = TRUE)
#Distance matrix
dmAzores = distMatrix(df2matrix(azores[2:3]))
dmAzores=dmAzores/1000  #distancia em km
dmAzores=dmAzores/100  #dist em centenas de Km
colnames(dmAzores) <- colnames(azoresData)
round(dmAzores,3)

##################
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
canary$richness = colSums(canaryData, na.rm = TRUE)/2 #/2 removes total from excel
canary$endemisms = colSums(canaryEndData, na.rm = TRUE)
#Distance matrix
dmCanary = distMatrix(df2matrix(canary[2:3]))
dmCanary=dmCanary/1000  #distancia em km
dmCanary=dmCanary/100  #dist em centenas de Km
colnames(dmCanary) <- colnames(canaryData)
round(dmCanary,3)

###########
#Modeling Species Richness under a Bayesian approach
###########
strAzores = azores[4:length(azores[1,])]
strAzores$`A (km2)`=log(strAzores$`A (km2)`)
for(i in 1:length(strAzores[1,])){
  strAzores[i] = (strAzores[,i]-mean(strAzores[,i]))/sd(strAzores[,i])
}
strCanary = canary[4:length(canary[1,])]
strCanary$`A (km2)`=log(strCanary$`A (km2)`)
for(i in 1:length(strCanary[1,])){
  strCanary[i] = (strCanary[,i]-mean(strCanary[,i]))/sd(strCanary[,i])
}

###
# Azores model: log(S)~ log(c) + z*log(A) + b1*elevation
# Area and elevation are standerized
##
mABoth <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a +  g[index] + b1*log(area) +  b2*altitude + b3*age,
    #log(lambda) <- a +  b1*log(area), 
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=azores$richness,
    area = azores$`A (km2)`,
    altitude = azores$`ELEV (m)`,
    age = azores$`Age (Ma)`,
    index=azores$index,
    Dmat=dmAzores),
  warmup=2000 , iter=6000 , cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))


###
# Canary model: log(S)~ log(c) + z*log(A) + b1*elevation + b2*age
# Area and elevation are standerized
##
mCBoth <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*log(area),  b2*altitude + b3*age,
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
    richness=canary$richness,
    area=canary$`A (km2)`,
    altitude = canary$`ELEV (m)`,
    age = canary$`Age (Ma)`,
    index=canary$index,
    Dmat=dmCanary),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

##
#Visualização dos dados
##
library("bayesplot")
library("rstanarm")
library("ggplot2")

postA <- extract.samples(mA.1)
#postA2 <- as.matrix(mA@stanfit)
#title <- ggtitle("Posterior distributions for Azores - 1st Edition")
#mcmc_areas(postA2,
#           pars = colnames(postA2)[-length(postA2[1,])],
#           prob = 0.95) + title

postC <- extract.samples(mC.1)
#postC2 <- as.matrix(mC@stanfit)
#title <- ggtitle("Posterior distributions for Canary Islands - 1st Edition")
#mcmc_areas(postC2,
#           pars = colnames(postC2)[-length(postC2[1,])],
#           prob = 0.95) + title

#Plot species area relationship
y1 = precis(mA.1)[1,1] + precis(mA.1)[2,1]*log(azores$`A (km2)`)
gg1 = precis(mA.1, depth = 2)[1:9,1]
#alt1 = precis(mA.4)[3,1] 
#age1 = precis(mA.1)[3,1]
yy1 = y1 + gg1 
  plot(log(azores$richness) ~ log(azores$`A (km2)`), xlab= "Ln (Area in km)", ylab ="Ln (Richness)")
  lines(precis(mA.1)[1,1]+precis(mA.1)[2,1]*log(azores$`A (km2)`)~log(azores$`A (km2)`), col="red")
  points(log(azores$`A (km2)`),yy1, pch=20, cex=0.5)
  title("Azores: Age and  Altitude model")

y2 = precis(mC.1)[1,1]  + precis(mC.1)[2,1]*log(canary$`A (km2)`)#+ precis(mC.1)[2,1]*log(canary$`A (km2)`)
gg2 = precis(mC.1, depth = 2)[2:8,1]
#alt2 = precis(mC.1)[3,1]
#age2 = precis(mC.1)[2,1]
#d2ML2= precis(mC.1)[4,1]
#yy2 = y2 #+ alt2 + age2 + d2ML2
plot(log(canary$richness) ~ log(canary$`A (km2)`), xlab= "Ln (Area in km)", ylab ="Ln (Richness)" )
  lines(precis(mC.1)[1,1]+precis(mC.1)[2,1]*log(canary$`A (km2)`)~log(canary$`A (km2)`), col="red")
  points(log(canary$`A (km2)`),yy2, pch=20, cex=0.5)
  title("Canary: Species-Area relationship")
  
  
#################################################################
#######
# Further down the rabbit hole: covariance analyses
#######
#AZORES
# compute posterior median covariance among islands
dm=dmAzores
K <- matrix(0,nrow=length(dm[1,]),ncol=length(dm[1,]))
for ( i in 1:length(dm[1,]) )
  for ( j in 1:length(dm[1,]) )
    K[i,j] <- median(postA$etasq) *
  exp( -median(postA$rhosq) * dm[i,j]^2)

diag(K) <- median(postA$etasq) + 0.01

# convert to correlation matrix
Rho <- round( cov2cor(K) , 2 )
# add row/col names for convenience
colnames(Rho) <- colnames(dm)
rownames(Rho) <- colnames(Rho)

# plot raw data and labels
plot( azores$Longitude , azores$Latitude , xlab="Longitude" , ylab="Latitude",
      col=rangi2 , pch=16 , xlim=c(min(azores$Longitude)-1,max(azores$Longitude)+1))
labels <- as.character(azores$`Archipelago/Island`)
text( azores$Longitude , azores$Latitude , labels=labels , cex=0.7 , pos=c(1,1,1,3,1,3,1,3,3) )
title("Azores correlations by geographic location")
# overlay lines shaded by Rho
for( i in 1:length(dm[1,]))
  for ( j in 1:length(dm[1,]))
    if ( i < j )
      lines( c( azores$Longitude[i],azores$Longitude[j] ) ,
             c( azores$Latitude[i],azores$Latitude[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2 ))

# compute posterior median relationship, ignoring distance
area.seq <- seq( from=1, to=9 , length.out=30 )           
lambda <- sapply( area.seq , function(lp) exp( postA$a + postA$b1*lp + postA$b3*lp))
lambda.median <- apply( lambda , 2 , median)
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.8 )
# plot raw data and labels
plot( log(azores$`A (km2)`) , log(azores$richness) , col=rangi2, pch=16 ,
      xlab="Ln(Area)" , ylab="Ln(Species Richness)",
      xlim=c(min(log(azores$`A (km2)`))-1,max(log(azores$`A (km2)`))+1),
      ylim=c(min(log(azores$richness))-1,max(log(azores$richness))+1)
)
text( log(azores$`A (km2)`) , log(azores$richness) , labels=labels , cex=0.7 ,
      pos=c(rep(1,7)))
title("Azores model fit")
# display posterior predictions
lines( area.seq , log(lambda.median), lty=2 )
lines( area.seq , log(lambda.PI80[1,]), lty=2 )
lines( area.seq , log(lambda.PI80[2,]), lty=2 )
# overlay correlations
for( i in 1:7 )
  for ( j in 1:7 )
    if ( i < j )
      lines( c( log(azores$`A (km2)`[i]),log(azores$`A (km2)`[j]) ) ,
             c( log(azores$richness[i]),log(azores$richness[j]) ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

#CANARY
# compute posterior median covariance among islands
dm=dmCanary
K <- matrix(0,nrow=length(dm[1,]),ncol=length(dm[1,]))
for ( i in 1:length(dm[1,]) )
  for ( j in 1:length(dm[1,]) )
    K[i,j] <- median(postC$etasq) *
  exp( -median(postC$rhosq) * dm[i,j]^2)

diag(K) <- median(postC$etasq) + 0.01

# convert to correlation matrix
Rho <- round( cov2cor(K) , 2 )
# add row/col names for convenience
colnames(Rho) <- colnames(dm)
rownames(Rho) <- colnames(Rho)

# plot raw data and labels
plot( canary$Longitude , canary$Latitude , xlab="Longitude" , ylab="Latitude",
      col=rangi2 , pch=16 , xlim=c(min(canary$Longitude)-1,max(canary$Longitude)+1))
labels <- as.character(canary$`Archipelago/Island`)
text( canary$Longitude , canary$Latitude , labels=labels , cex=0.7 , pos=c(1,1,1,3,1,3,1,3,3) )
title("Canary correlations by geographic location")
# overlay lines shaded by Rho
for( i in 1:length(dm[1,]))
  for ( j in 1:length(dm[1,]))
    if ( i < j )
      lines( c( canary$Longitude[i],canary$Longitude[j] ) ,
             c( canary$Latitude[i],canary$Latitude[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2 ))

# compute posterior median relationship, ignoring distance
area.seq <- seq( from=4, to=9 , length.out=30 )           
lambda <- sapply( area.seq , function(lp) exp( postC$a + postC$b1*lp + postC$b2*lp + postC$b3*lp)) #+ postC$b4*lp))
lambda.median <- apply( lambda , 2 , median)
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.8 )
# plot raw data and labels
plot( log(canary$`A (km2)`) , log(canary$richness) , col=rangi2, pch=16 ,
      xlab="ln(Area)" , ylab="ln(Species Richness)",
      xlim=c(min(log(canary$`A (km2)`))-1,max(log(canary$`A (km2)`))+1),
      ylim=c(min(log(canary$richness))-2,max(log(canary$richness))+1)
)
text( log(canary$`A (km2)`) , log(canary$richness) , labels=labels , cex=0.7 ,
      pos=c(rep(1,7)))
title("Canary Islands model fit")
# display posterior predictions
lines( area.seq , log(lambda.median), lty=2 )
lines( area.seq , log(lambda.PI80[1,]), lty=2 )
lines( area.seq , log(lambda.PI80[2,]), lty=2 )
# overlay correlations
for( i in 1:7 )
  for ( j in 1:7 )
    if ( i < j )
      lines( c( log(canary$`A (km2)`[i]),log(canary$`A (km2)`[j]) ) ,
             c( log(canary$richness[i]),log(canary$richness[j])) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]) )

