##
# Bayesian regression analyses of the number of species in Azores and Canary Islands
# using STAN through the package rethinking, by Mc Elreath
#
#
# Author: Diogo Barros
# Version: 1.02.b*
#
#* ADAPTED TO MAC OS
##

source("/Users/diogobarros/Documents/BEAG/R scripts/Task1/myFunctions2.R")
setwd("/Users/diogobarros/Documents/BEAG/Dados/Paper1")
library(rethinking)
library(readxl)
library(dplyr)
options(mC.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#########
#Preparing the data
#########

#AZORES
azoresFull <- as.data.frame(read_excel("azores_arthopoda.xlsx", sheet=1))
azoresVars <- as.data.frame(read_excel("MAC_EnvData.xlsx", sheet=1, col_names = FALSE))
azoresEnd <- filter(azoresFull, COL.=='END')
azores <- azoresVars[9:17,]
azoresNames = azores[,1]  
azores <- as.data.frame(df2matrix(azores))
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
azores$endNoSIE = azores$endmisms-azores$SIE
#Distance matrix
dmAzores = distMatrix(df2matrix(azores[2:3]))
dmAzores=dmAzores/100000  #distance in hundreds of km
colnames(dmAzores) <- azores$`Archipelago/Island`
rownames(dmAzores) <- azores$`Archipelago/Island`
round(dmAzores,3)


#CANARY ISLANDS
canaryFull <- as.data.frame(read_excel("ArthropodaCanaryIslands.xls", sheet=1))
canaryVars <- as.data.frame(read_excel("MAC_EnvData.xlsx", sheet=1, col_names = FALSE))
canaryEnd=filter(canaryFull, END=='End')
canary <- canaryVars[36:42,]
canaryNames = canary[,1]
canary <- as.data.frame(df2matrix(canary))
canary$V1 <- canaryNames
colnames(canary) <- canaryVars[7,]
canary$index <- 1:7
#Getting species richness
canaryData<-canaryFull[18:length(canaryFull)]
canaryData=data.frame(canaryData$fuerteventura, canaryData$lanzarote,canaryData$canaria,canaryData$tenerife,
                      canaryData$gomera,canaryData$hierro,canaryData$palma)
canaryEndData <-canaryEnd[18:length(canaryEnd)] 
canaryEndData=data.frame(canaryEnd$fuerteventura, canaryEnd$lanzarote,canaryEnd$canaria,canaryEnd$tenerife,
                         canaryEnd$gomera,canaryEnd$hierro,canaryEnd$palma)
canarySIEData <- filter(canaryEndData, rowSums(canaryEndData, na.rm = TRUE) == 1)
canary$richness = colSums(canaryData, na.rm = TRUE)/2 #/2 removes total from excel
canary$endemisms = colSums(canaryEndData, na.rm = TRUE)
canary$SIE = colSums(canarySIEData, na.rm = TRUE)
canary$endNoSIE = canary$endemisms - canary$SIE
#Distance matrix
dmCanary = distMatrix(df2matrix(canary[2:3]))
dmCanary=dmCanary/100000  #distance in hundreds of km
colnames(dmCanary) <- canary$`Archipelago/Island`
rownames(dmCanary) <- canary$`Archipelago/Island`
round(dmCanary,3)


###########
#Modeling Species Richness under a Bayesian approach
###########


island=azores
dmisland=dmAzores
#island=canary
#dmisland=dmCanary
richness=island$endNoSIE

###
#Model 1: spp~area+ distance
###
mA.end2 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*log(area),
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=island$`A (km2)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

pairs(mA.All)
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
    area=log(island$`A (km2)`),
    age=island$`Age (Ma)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 , chains=4,cores=4,control=list(max_treedepth =15,adapt_delta=0.9))

  ### 
#Model 3: spp ~ area + distance + altitude
###

mC.3 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*log(area) + b3*altitude,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b3 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=island$`A (km2)`,
    altitude=island$`ELEV (m)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 , chains=4,cores=4, control=list(max_treedepth =15,adapt_delta=0.9))


###
#Model 4: spp ~ area + distance + age + altitude 
###

mC.4 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*area +  b2*age + b3*altitude,
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
    area=log(island$`A (km2)`),
    altitude=island$`ELEV (m)`,
    age=island$`Age (Ma)`,
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
    log(lambda) <- a + g[index] + b1*area + b4*d2MLand,
    g[index] ~ GPL2(Dmat , etasq , rhosq , sigma ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b4 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=log(island$`A (km2)`),
    d2MLand = island$`DM (km)`/100,
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
    log(lambda) <- a + g[index] + b1*area  + b2*age + b3*altitude + b4*d2MLand,
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
    area=log(island$`A (km2)`),
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
    log(lambda) <- a + b1*area,
    a ~ dnorm(0,100),
    b1 ~ dnorm(0,1)
  ),
  data=list(
    richness=richness,
    area=log(island$`A (km2)`)),
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
    log(lambda) <- a + b1*area, #+b2*age+b3*altitude+b4*d2MLand,
    a ~ dnorm(0,100),
    b1 ~ dnorm(0,1)
    #c(b1,b2,b3,b4) ~ dnorm(0,1)
  ),
  data=list(
    richness=richness,
    area=log(island$`A (km2)`)),
    #age=island$`Age (Ma)`,
    #altitude = island$`ELEV (m)`,
    #d2MLand = island$`DM (km)`/100),
  warmup=2000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

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

###
#Model 11: GDM MODEL (BORREGARD ET AL, 2008)  diversity = [log]Area + Time + Time2
###
mC.GDM<- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + b1*log(area) + b2*age + b3*age^2,
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~ dnorm(0,1)
  ),
  data=list(
    richness=richness,
    area=island$`A (km2)`,
    age=island$`Age (Ma)`),
  warmup=4000 , iter=1e4 ,cores=4, chains=4,control=list(max_treedepth =15,adapt_delta=0.9))

###
#Model 12: GDM MODEL + SA 
###

mC.GDMSA <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + b1*log(area) + b2*age + b3*age^2,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    b1 ~ dnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=richness,
    area=island$`A (km2)`,
    age=island$`Age (Ma)`,
    index=island$index,
    Dmat=dmisland),
  warmup=4000 , iter=1e4 , chains=4,cores=4,control=list(max_treedepth =15,adapt_delta=0.9))




####
#OUTPUT
####

pairs(azores[4:8][-3], pch=19, main="Correlation pairs between variables")

#Create list with all models for 1 island (specified by the dataframe "island")
if (island[1,1] == "Corvo") myPat = "^mA." else myPat = "^mC."
modelNames = ls()[grep(myPat, ls())]
models = c(get(modelNames[1]))

if (length(modelNames)>1)
  for (i in 2:length(modelNames)) {
    modelsIsland = append(models, get(modelNames[i]))
  }
    
#Compare WAICS and wite to .csv. PLEASE SPECIFY NAME
waics=compare(mC.1,mC.2,mC.3,mC.4,mC.5,mC.6,mC.NGP1,mC.Null,mC.SCAZ,mC.10)
write.table(waics,"/home/splinter/Documents/BEAG/Results/Paper1/waicsAzoresTotal2.csv")

#######
# Covariance analyses
#######
#Prep
if (length(modelNames==1)) bestModel = models[[1]] else bestModel = get(rownames(waics[1,]))
nVars = length(precis(bestModel)[,1])-2
post <- extract.samples(bestModel)

#Richness ~ area graphs
y = precis(bestModel)[1,1] + precis(bestModel)[2,1]*log(island$`A (km2)`)
gg <- alt <- age <- d2ML <- 0
if(length(precis(bestModel, depth = 2)[,1])>length(island$`Archipelago/Island`)){
  gg = precis(bestModel, depth = 2)[1:length(island$`Archipelago/Island`),1]
}
if("age" %in% names(bestModel@data)){
  age = mean(post$b2)
}
if("altitude" %in% names(bestModel@data)){
  alt = mean(post$b3)
}
if("d2MLand" %in% names(bestModel@data)){
  d2ML = mean(post$b4)
}
#Final equation
yy = y + gg + age + alt + d2ML
plot(log(richness) ~ log(island$`A (km2)`), xlab= "Ln (Area in km)", ylab ="Ln (Richness)")   
if("area" %in% names(bestModel@data))
  lines(mean(post$a)+mean(post$b1)*log(island$`A (km2)`) ~ log(island$`A (km2)`), col="red")
  a = lm(log(richness)~ log(island$`A (km2)`))
  lines(a, col="blue")
  abline(a=precis(bestModel)[1,1], b=precis(bestModel)[2,1], col="blue")
#else if ("index" %in% names(bestModel@data)) {
# lines(mean(post$a), col="red")
#}
points(log(island$`A (km2)`),yy, pch=20, cex=0.5)
title(sprintf("Azores species-area relationship\n for endemic species (without SIE) "))

#####
# Geographic correlations
#####
dm=dmisland
K <- matrix(0,nrow=length(dm[1,]),ncol=length(dm[1,]))
for ( i in 1:length(dm[1,]) )
  for ( j in 1:length(dm[1,]) )
    K[i,j] <- median(post$etasq) *
  exp( -median(post$rhosq) * dm[i,j]^2)
diag(K) <- median(post$etasq) + 0.01

# convert to correlation matrix
Rho <- round( cov2cor(K) , 2 )
# add row/col names for convenience
colnames(Rho) <- colnames(dm)
rownames(Rho) <- colnames(Rho)
write.csv(Rho, "/home/splinter/Documents/BEAG/Results/Paper1/CanarySIE.csv")

# plot raw data and labels
plot( island$Longitude , island$Latitude , xlab="Longitude" , ylab="Latitude",
      col=rangi2 , pch=16 , xlim=c(min(island$Longitude)-1,max(island$Longitude)+1))
labels <- as.character(island$`Archipelago/Island`)
text( island$Longitude , island$Latitude , labels=labels , cex=0.7 , pos=c(1,1,1,3,1,3,1,3,3) )
title("Azores correlations\n by geographic location")
# overlay lines shaded by Rho
for( i in 1:length(dm[1,]))
  for ( j in 1:length(dm[1,]))
    if ( i < j )
      lines( c( island$Longitude[i],island$Longitude[j] ) ,
             c( island$Latitude[i],island$Latitude[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2))

##
#  Species - Area Correlations
##  

# compute posterior median relationship, ignoring distance
area.seq <- seq( from=1, to=9 , length.out=30 )           
lambda <- sapply( area.seq , function(lp) exp(post$a + post$b1*lp))
lambda.median <- apply( lambda , 2 , median)
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.9 )
# plot raw data and labels
plot( log(azores$`A (km2)`) , log(richness) , col=rangi2, pch=16 ,
      xlab="Ln(Area)" , ylab="Ln(Species Richness)",
      xlim=c(min(log(azores$`A (km2)`))-1,max(log(azores$`A (km2)`))+1),
      ylim=c(min(log(richness))-1,max(log(richness))+1)
)
text( log(azores$`A (km2)`) , log(richness) , labels=labels , cex=0.7 ,
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
      lines( c( log(island$`A (km2)`[i]),log(island$`A (km2)`[j]) ) ,
             c( log(richness[i]),log(richness[j]) ) ,
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
lambda <- sapply( area.seq , function(lp) exp( postC$a + postC$b1*lp + postC$b3*lp))# + postC$b3*lp))
lambda.median <- apply( lambda , 2 , median)
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.8 )
# plot raw data and labels
plot( log(canary$`A (km2)`) , log(canary$richness) , col=rangi2, pch=16 ,
      xlab="ln(Area)" , ylab="ln(Species Richness)",
      xlim=c(min(log(canary$`A (km2)`))-1,max(log(canary$`A (km2)`))+1),
        ylim=c(min(log(canary$richness))-1,max(log(canary$richness))+1)
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
               c( log(canary$richness[i]),log(canary$richness[j]) ) ,
               lwd=2 , col=col.alpha("black",Rho[i,j]^2) )




####
# Creating common species tables
####

#SERECT UR CUR
common = islandSIEData
commonSpp = matrix(nrow=length(island$`Archipelago/Island`),ncol=length(island$`Archipelago/Island`),0)
rownames(commonSpp)=island$`Archipelago/Island`
colnames(commonSpp)=island$`Archipelago/Island`

for (k in 1:length(common[1,])){
  for (j in 1:length(common[1,])){
    for (i in 1:length(common[,1])){
      if(common[i,k] == 1 && common[i,j] == 1) {
        commonSpp[k,j]=commonSpp[k,j]+1
      }
    }
  }
}
#REMEMBER TO CHANGE THE FILE NAME ACORDING TO THE TYPE OF RICHNESS 
write.table(commonSpp, "/home/splinter/Documents/BEAG/Results/Paper1/commonSppSIE.csv")


plot(island$age~island$)


##############################################################################################################
#Garbage/junk code

titles = c("Spatial Autocorrelation", "Age", "Altitude", "Age and Alitude", "Distance to Mainland",
           "Age, Altitude and Distance to Mainland", "No gaussian process", "Null model", "simple species-area",
           "only Spatial Autocorrelation")
#Richness ~area graphs
##
for(i in 1:length(modelsAzores)){
  post <- extract.samples(modelsAzores[[i]])
  y = precis(modelsAzores[[i]])[1,1] + precis(modelsAzores[[i]])[2,1]*log(canary$`A (km2)`)
  gg <- alt <- age <- d2ML <- 0
  if(length(precis(modelsAzores[[i]], depth = 2)[,1])>7){
    gg = precis(modelsAzores[[i]], depth = 2)[1:7,1]
  }
  if("age" %in% names(modelsAzores[[i]]@data)){
    age = mean(post$b2)
  }
  if("altitude" %in% names(modelsAzores[[i]]@data)){
    alt = mean(post$b3)
  }
  if("d2MLand" %in% names(modelsAzores[[i]]@data)){
    d2ML = mean(post$b4)
  }
  yy = y + gg + age + alt + d2ML
  plot(log(canary$richness) ~ log(canary$`A (km2)`), xlab= "Ln (Area in km)", ylab ="Ln (Richness)")
  if("area" %in% names(modelsAzores[[i]]@data)){
    lines(mean(post$a)+mean(post$b1)*log(canary$`A (km2)`)~log(canary$`A (km2)`), col="red")
  }else if ("index" %in% names(modelsAzores[[i]]@data)) {
    lines(mean(post$a), col="red")
  }
  points(log(canary$`A (km2)`),yy, pch=20, cex=0.5)
  title(sprintf("Canary: %s model", titles[i]))
}

lines(mean(post$a)+mean(post$b1)*log(canary$`A (km2)`)~log(canary$`A (km2)`), col="red")
##
# Geographic correlations
##
for(i in 1:length(modelsAzores)){
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
}

##
#  Species - Area Correlations
##  
for(i in 1:length(modelsAzores)){
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
  lambda <- sapply( area.seq , function(lp) exp( postC$a + postC$b1*lp + postC$b3*lp))# + postC$b3*lp))
  lambda.median <- apply( lambda , 2 , median)
  lambda.PI80 <- apply( lambda , 2 , PI , prob=0.8 )
  # plot raw data and labels
  plot( log(canary$`A (km2)`) , log(canary$richness) , col=rangi2, pch=16 ,
        xlab="ln(Area)" , ylab="ln(Species Richness)",
        xlim=c(min(log(canary$`A (km2)`))-1,max(log(canary$`A (km2)`))+1),
        ylim=c(min(log(canary$richness))-1,max(log(canary$richness))+1)
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
               c( log(canary$richness[i]),log(canary$richness[j]) ) ,
               lwd=2 , col=col.alpha("black",Rho[i,j]^2) )
}
