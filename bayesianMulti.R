###
#BAYSIAN MODELING WITH STAN:
#ABUNDANCE OF DIFFERENT SPECIES
###

setwd("C:/Users/Diogo Barros/Documents/Diogo/R data")
data <- read.csv("landscape.csv", sep=";")
attach(data) 
library("rstan") # observe startup messages
library("beepr") #BEEP
library("mailR") #EMAIL
rstan_options(auto_write = TRUE) # for multiple cores
options(mc.cores = parallel::detectCores())

###
#Running STAN with mutiple variables - Park metrics
###

#Variables - 
y <- CR #Response

buffer<-cbind(data$F100.,F500.,F1000.) #Buffer matrix

mVar <- cbind(rep(1,11),
              GreenUsable,
              Habitats,
              ShapeInd,
              rep(1,11))

colnames(mVar)[1]<-"Intercept"

waics<-matrix(1,5,3)


#ADDING EXTRA VARIABLES 
for(i in 1:length(buffer[1,])){
  mVar[,5] <- buffer[i] 
  
  file <-"pois_ln_regression_WAIC.stan"
  source("C:/Users/Splinterq/Documents/Diogo/The future/BEGC/R scripts/Task1/stanCaller.R")
    
  waics[,i]<- myStan(y,mVar,file)
  print(i)
  if(i == length(buffer[1,]))
    beep(sound=3, expr = NULL)
    send.mail(from = "dingobarros@gmail.com",
            to = "Eu <dingobarros@gmail.com>",
            subject = "Modulation complete",
            body = "IT IS DONE!!!11",
            smtp = list(host.name = "aspmx.l.google.com", port = 25),
            authenticate = FALSE, send = TRUE, debug=TRUE)
}

waics

###
#Running STAN with mutiple variables - Park metrics
###
y <- CR #Response

vars<-cbind(People,Animals,soundCentr) #Buffer matrix
mVar <- cbind(rep(1,11),
              GreenUsable)
colnames(mVar)[1]<-"Intercept"
waics<-matrix(1,5,3)

for(i in 1:length(vars[1,])){
  mVar <- cbind(mVar,vars[,i])
  colnames(mVar)[2+i]<- names(vars[1,i])
  
  file <-"pois_ln_regression_WAIC.stan"
  source("C:/Users/Splinterq/Documents/Diogo/The future/BEGC/R scripts/Task1/stanCaller.R")
  
  waics[,i]<- myStan(y,mVar,file)
  print(i)
  if(i == length(buffer[1,]))
    beep(sound=3, expr = NULL)
    send.mail(from = "dingobarros@gmail.com",
            to = "Eu <dingobarros@gmail.com>",
            subject = "Modulation complete",
            body = "IT IS DONE!!!11",
            smtp = list(host.name = "aspmx.l.google.com", port = 25),
            authenticate = FALSE, send = TRUE, debug=TRUE)
}

waics





###
#Running STAN with a singe variable
###

#Variables


y2 <- MM #Response
mVar2 <- cbind(rep(1,11),
               Area,
               Perimeter,
               NrFrags,
               ShrubArea,
               HerbsArea,
               NakedArea,
               MonsantoDist,
               ClosestDist,
               LightMax,
               soundCentr,
               soundMarg,
               People,
               GreenUsable,
               LogGreen,
               Habitats,
               ShapeInd,
               F100,
               F500,
               F1000)



waics2<-c(rep(1,length(mVar2[1,])))
file <-"pois_ln_regression_WAIC.stan"
source("C:/Users/Splinterq/Documents/Diogo/The future/BEGC/R scripts/Task1/stanCaller.R")
waics2<- myStan2(y2,mVar2,file)

write.csv(waics2, file="waicsTemp.csv")  

