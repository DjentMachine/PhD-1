##
#Gaussian Processeses
##

source("/home/splinter/Documents/BEGC/R scripts/Task1/myFunctions2.R")
setwd("/home/splinter/Documents/R data/Azores")
library(rethinking)
library(readxl)


#Preparing the data
fullData <- as.data.frame(read_excel("BAD_LBA_Final_Calculos.xlsx", sheet=1))
islandVars <- as.data.frame(read_excel("SITES_BALA_100_UTMs.xlsx", sheet=3))
iData <- as.data.frame(read_excel("islandCoords.xlsx"))
iData <- iData[-c(1,6),]
iData$index <- 1:7
coords <- iData[6:7]
iNames <-c("FLO","FAI","PIC","SJG","TER","SMG","SMR")
iData$richness <- c(rep(0,length(iNames)))
iData$maxAlt <- round(as.numeric(islandVars$..4[-c(1,2,7)]),2)
iData$meanAlt <- round(as.numeric(islandVars$..3[-c(1,2,7)]),2)
iData

#Getting species richness


#Distance matrix
dm = distMatrix(coords)
dm=dm/1000  #distancia em km
dm=dm/100  #dist em centenas de Km
colnames(dm) <- c("FLO","FAI","PIC","SJG","TER","SMG","SMR")
round(dm,3)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
###
# Model 1: spp~area+ distance
# Varying intercepts for each society
###

### Normal area
m1.1 <- map2stan(
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
    richness=iData$richness,
    area=iData$Area,
    index=iData$index,
    Dmat=dm),
  warmup=2000 , iter=1e4 , chains=4)

    
### log area
m1.2 <- map2stan(
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
    richness=iData$richness,
    area=log(iData$Area),
    index=iData$index,
    Dmat=dm),
  warmup=4000 , iter=1e4 , chains=4)


###
#Model 2: spp ~ area + distance + age
###
m2 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + bp*area + bc*age,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    bc ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area),
    index=iData$index,
    age = (iData$`Age (MY)`-mean(iData$`Age (MY)`))/sd(iData$`Age (MY)`),#iData$`Age (MY)`,
    Dmat=dm),
    warmup=2000 , iter=1e4 , chains=4)

### 
#Model 3: spp ~ area + distance + altitude
###

m3.2 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + bp*area + bc*altitude,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    bc ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area),
    index=iData$index,
    #altitude = (iData$maxAlt-mean(iData$maxAlt))/sd(iData$maxAlt),
    altitude = (iData$meanAlt-mean(iData$meanAlt))/sd(iData$meanAlt),
    Dmat=dm),
    warmup=2000 , iter=1e4 , chains=4)


###
#Model 4: spp ~ area + distance + age + altitude 
###

m4 <- map2stan(
  alist(
    richness ~ dpois(lambda),
    log(lambda) <- a + g[index] + bp*area + bc*altitude + bd*age,
    g[index] ~ GPL2(Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    bc ~ dnorm(0,1),
    bd ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area),
    index=iData$index,
    age = (iData$`Age (MY)`-mean(iData$`Age (MY)`))/sd(iData$`Age (MY)`),
    altitude = (iData$maxAlt-mean(iData$maxAlt))/sd(iData$maxAlt),
    Dmat=dm),
  warmup=2000 , iter=1e4 , chains=4)




#######
# Further down the rabbit hole: covariance analyses
#######
# compute posterior median covariance among islands
post = extract.samples(m1.2)

K <- matrix(0,nrow=7,ncol=7)
for ( i in 1:7 )
  for ( j in 1:7 )
    K[i,j] <- median(post$etasq) *
              exp( -median(post$rhosq) * dm[i,j])

diag(K) <- median(post$etasq) + 0.01

# convert to correlation matrix
Rho <- round( cov2cor(K) , 2 )
# add row/col names for convenience
colnames(Rho) <- colnames(dm)
rownames(Rho) <- colnames(Rho)
Rho

# scale point size to logpop
psize <- iData$Area / max(iData$Area)*4
#psize <- exp(psize*1.5)-2
# plot raw data and labels
plot( iData$`Longitude  (WGS84)` , iData$`Latitude (WGS84)` , xlab="longitude" , ylab="latitude",
      col=rangi2 , cex=psize, pch=16 , xlim=c(-32,-24) )
labels <- as.character(iData$Island)
text( iData$`Longitude  (WGS84)` , iData$`Latitude (WGS84)` , labels=labels , cex=0.7 , pos=c(1,1,1,3,3,3,3) )
# overlay lines shaded by Rho
for( i in 1:7)
  for ( j in 1:7)
    if ( i < j )
      lines( c( iData$`Longitude  (WGS84)`[i],iData$`Longitude  (WGS84)`[j] ) ,
             c( iData$`Latitude (WGS84)`[i],iData$`Latitude (WGS84)`[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )


# compute posterior median relationship, ignoring distance
area.seq <- seq( from=3, to=7 , length.out=30 )
lambda <- sapply( area.seq , function(lp) exp( post$a + post$bp*lp ) )
lambda.median <- apply( lambda , 2 , median )
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.8 )
# plot raw data and labels
plot( log(iData$Area) , log(iData$richness) , col=rangi2 , cex=psize , pch=16 ,
      xlab="log area" , ylab="species richness" )
text( log(iData$Area) , log(iData$richness) , labels=labels , cex=0.7 ,
      pos=c(4,3,4,2,2,1,4,4,4,2) )
# display posterior predictions
lines( area.seq , log(lambda.median) , lty=2 )
lines( area.seq , log(lambda.PI80[1,]) , lty=2 )
lines( area.seq , log(lambda.PI80[2,]) , lty=2 )
# overlay correlations
for( i in 1:7 )
  for ( j in 1:7 )
    if ( i < j )
      lines( c( log(iData$Area[i]),log(iData$Area[j]) ) ,
             c( log(iData$richness[i]),log(iData$richness[j]) ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )










##################################################
plot(precis(m1,depth=2))
post <- extract.samples(m1)



###
#No Gaussian Process, bayesian inference
###

###
#Model 1: simple Spp~area regression
###

mNGP1 <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a + bp*area,
    a ~ dnorm(0,100),
    bp ~ dnorm(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area)),
    #age=iData$`Age (MY)`),
  warmup=2000 , iter=1e4 , chains=4)

###
#Model 2: spp~area+age
###

mNGP2 <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a + bp*area+bc*age,
    a ~ dnorm(0,100),
    c(bp,bc) ~ dnorm(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area),
    age=iData$`Age (MY)`),
  warmup=2000 , iter=1e4 , chains=4)

###
#Model 3: spp~area+altitude
###

mNGP3 <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a + bp*area+bc*altitude,
    a ~ dnorm(0,100),
    c(bp,bc) ~ dnorm(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area),
    altitude = (iData$maxAlt-mean(iData$maxAlt))/sd(iData$maxAlt)),
  warmup=2000 , iter=1e4 , chains=4)

###
#Model 4: spp~area+altitude+age
###

mNGP4 <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a + bp*area+bc*altitude+bd*age,
    a ~ dnorm(0,100),
    c(bp,bc,bd) ~ dnorm(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area),
    age=iData$`Age (MY)`,
    altitude = (iData$maxAlt-mean(iData$maxAlt))/sd(iData$maxAlt)),
  warmup=2000 , iter=1e4 , chains=4)



###
#Model 5: spp~area+altitude+age
###

mNGP5 <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a + bp*area+bc*altitude+bd*age,
    a ~ dnorm(0,100),
    c(bp,bc,bd) ~ dnorm(0,1)
  ),
  data=list(
    richness=iData$richness,
    area=log(iData$Area),
    age=iData$`Age (MY)`,
    altitude = (iData$maxAlt-mean(iData$maxAlt))/sd(iData$maxAlt)),
  warmup=2000 , iter=1e4 , chains=4)

###
#Model 6: null model
###

mNull <- map2stan(
  alist(
    richness ~ dpois( lambda ),
    log(lambda) <- a, 
    a ~ dnorm(0,100)
  ),
  data=list(
    richness=iData$richness),
  warmup=2000 , iter=1e4 , chains=4)


a=compare(m1.2,m2,m3,m4,mNGP1,mNGP2,mNGP3,mNGP4,mNGP5,mNull)
b=as.data.frame(a@output)
write.csv(b, "teste.csv")






# plot the posterior median covariance function
curve( median(post$etasq)*exp(-median(post$rhosq)*x^2) , from=0 , to=3 ,
       xlab="distance (hundreds km)" , ylab="covariance" , ylim=c(0,1) ,
       yaxp=c(0,1,4) , lwd=2 )
# plot 100 functions sampled from posterior
for ( i in 1:100 )
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE ,
         col=col.alpha("black",0.2))
