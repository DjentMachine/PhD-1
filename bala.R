###
# Dados BALA - Paulo Borges
#
# 5 ilhas => Faial (8T), Flores(12T), Pico(16T),
###

setwd("C:/Users/Diogo/Documents/Diogo/R data")
fullData <- read.csv("bala.csv", sep=",")
vars <- read.csv("balaVars.csv",sep=",")


#MODEL FUNCTIONS

##
#Frequentist Speies-area regression: log(S)=log(c)+zlog(A)
##
freqRegression <- function(S,A){
  lm1 = lm(log(S)~log(A))
  #cat(sprintf("C is %s \nZ is %f", exp(1)^lm1$coefficients[1], lm1$coefficients[2]))
  return(lm1)
}

##  
#Bayesian Species-area regression:
##
bayesianRegression <- function(S,A,file){
  stan_data <- list(N=length(S),
                    y=log(S),
                    x=log(A))
  fit <- stan(file=file,data=stan_data,iter=5000,chain=4,
              cores = 4, 
              control=list(max_treedepth =20,adapt_delta=0.95))
  return(fit)
}


##
#Getting data by Island
##

#Acumultion curves
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE)  
options(mc.cores = parallel::detectCores()) # for multiple cores
library(vegan)
library(beepr)
#Parameters
a�oresI = c("FAI","FLO","PIC","SJG","SMG","SMR","TER")
resultsM = matrix(nrow=7, ncol=4) 
colnames(resultsM) = c("c.F","z.F","c.B","z.B")
rownames(resultsM) = a�oresI
aCurves = vector("list",length(a�oresI))
fModels = vector("list",length(a�oresI))
bModels = vector("list",length(a�oresI))
file = "regressionSIMPLE.stan"
for(i in 1:length(aCurves)){
  #Selecting data
  selection = grep(a�oresI[i],colnames(fullData))
  selection = selection[-length(selection)]
  iData = fullData[selection]
  iData = iData[1:(length(iData[,1])-1),] 
  #Saving Acurve in vector and estimates of C and Z in matrix
  aCurves[[i]]  = specaccum(t(iData), method="random")
  fModels[[i]]  = freqRegression(aCurves[[i]]$richness,c(1:length(iData[1,])))
  resultsM[i,1] = exp(fModels[[i]]$coefficients[1])
  resultsM[i,2] = fModels[[i]]$coefficients[2] 
  bModels[[i]]  = bayesianRegression(aCurves[[i]]$richness,c(1:length(iData[1,])),file)
  resultsM[i,3] = exp(mean(extract(bModels[[i]])[[2]]))
  resultsM[i,4] = mean(extract(bModels[[i]])[[1]])
  sprintf("Set %s finished, commencing set %s",i,i+1)
  if(i==length(aCurves))beep(sound=3)  
}

##
#Hierarchical models
##

#Exploratory: getting posterior histograms, fitting normal curve
M=1
posterior = extract(bModels[[M]])
for (i in 1:(length(names(posterior))-2)){ 
  h <- hist(posterior[[i]], main = sprintf("%s : Posterior of %s",a�oresI[M], names(posterior)[i]),
            xlab = sprintf("Value of %s", names(posterior)[i]))
  xfit <- seq(min(posterior[[i]]), max(posterior[[i]]), length = 40) 
  yfit <- dnorm(xfit, mean = mean(posterior[[i]]), sd = sd(posterior[[i]])) 
  yfit <- yfit * diff(h$mids[1:2]) * length(posterior[[i]]) 
  lines(xfit, yfit, col = "black", lwd = 2)
}
bModels[[M]]

#All islands, priors = z~normal(0.466,0.002)  ; logC~(3.651,0.037)
# prior parameters got from mean estimations
islands <- read.csv("balaISLANDS.csv", sep=",")
len = length(islands[,1])
#Getting island area in one vector
AI = c(rep(-1,length(islands[1,])))
for(i in 1:length(AI)){
  AI[i] = islands[length(islands[,1]),i]   
}
#Getting island Richness in one vector
SI = AI
for(i in 1:length(SI)){
  SI[i] = islands[len-1,i]   
}
#Accumulation curve for ALL species on ALL islands
allACurve = specaccum(t(islands), method="random")

#Getting BETA diversity in one vector
bDiv = AI
for(i in 1:length(SI)){
  SI[i] = islands[len-1,i]   
}

#Finding out prior information
c = resultsM[,3]
z = resultsM[,4]
hist(c, breaks=7)
hist(z, breaks=7)


#Go Hierarchical model go!\
S=allACurve$richness
A=c(1,2,3,4,5,6,7)
stan_data <- list(N=length(S),
                  y=log(S),
                  x=log(A))
fit <- stan(file="regressionH2.stan",data=stan_data,iter=10000,chain=4,
            cores = 4, 
            control=list(max_treedepth =20,adapt_delta=0.9999))

fit
hist(extract(fit)$z)
hist(extract(fit)$z[,1])
fit

freqRegression(S,A)




#Getting SA relatioship variablles
ST= rep(0,length(TER[1,]))
for(i in 1:length(ST)){
  ST[i] = as.numeric(TER[length(TER[,1]),i])  
}

AT = vars$Area.Fragment..ha.[61:99]
AT= 0.01*AT

##
#Extracting results
#TARGET: z= 0.22, logC=15.4
#saveRDS(fit, "normals.rds")  #GRAVA MODELO
#READING PREVIOUS MODELS 


fit = readRDS("uniforms.rds")
y = S
-x = A
posterior = extract(fit) #Full posterior of parameters

#plotting the posterior distribution for the parameters
plot(y~x, pch=20, xlab="Log")
abline( mean(posterior$logC), mean(posterior$z))
abline(lm1, col = 2, lty = 2, lw = 3)
for (i in 1:500) {
-  abline(posterior$logC[i], posterior$z[i], col = "gray", lty = 1)
}
plot(exp(y)~exp(x), pch=20, xlab="Log",
     xlim=c(0,max(exp(x[length(x)-1]))))


#Fit S=CA^Z
lines(exp(x),exp(mean(posterior$logC))*exp(x)^mean(posterior$z),col="red")




