###
#BAYSIAN MODELING WITH STAN
###


setwd("C:/Users/Diogo/Documents/Diogo/R data")
data <- read.csv("landscape.csv", sep=",")
attach(data) #Dependent variable
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE) # for multiple cores
options(mc.cores = parallel::detectCores())


#Parameters
y= RMM
n <- 11
p <- 2
X <- cbind(rep(1,11),
      GreenUsable/10000)
      #LogGreen)
      #F100.,
      #F500)
      #F1000)
      #Habitats
      #ShapeInd)
      #People)
      
log_offset <- 1

stan_data <- list(y=y,
                  n=n,
                  p=p,
                  X=X)


fit <- stan(file="beta4.stan",       
            #file="pois_ln_regression_WAIC.stan",
            data=stan_data,iter=500,chains=1,
            control=list(max_treedepth =10),
            init=0)
    
##
#Extracting results
##
  
print(fit)
results <-summary(fit)

#Selecting the Lambda estimations 
a<-extract(fit)
lambdas <- a[["lambda"]]

#Sampling Lambdas
source("C:/Users/Diogo/Documents/Diogo/The future/BEGC/R scripts/Task1/myFunctions.R")
sampleSize<-1000
sample <- sampler(lambdas,sampleSize,n)  #amostra 1 vez do output do STAN
s.sample  <- s.gen(sample)               #Generate a matrix of S using Poiss(lambda)

#Calculating Likelihood and WAIC  
waic1<-waic(fit)
waic1
saveRDS(fit, "teste.rds")  #GRAVA MODELO
teste <- readRDS("teste.rds")      #LÊ MODELO

#####
#Plots
#####

#Plots for general info
#plot(fit, pars =names(fit))
#traceplot(fit, pars =names(fit))

#Histograms for 11 lambda distributions 
par(mfrow=c(3,4))
for(i in 1:dim(a$lambda)[2]){
  hist(a$lambda[,i],main=data$Name[i],xlim=c(0,12),xlab="Lambdas")
  abline(v=Richness[i],col="red")
  m <- paste("mean=",round(mean(a$lambda[,i]),2))
  s <- paste("sd=",round(sd(a$lambda[,i]),2))
  legend("topright",c(m,s),bty="n")
}

#Histograms for the 11 Richness estimates
par(mfrow=c(3,4))
for(i in 1:11){
  hist(s.sample[,i],main=data$Name[i],xlim=c(0,12),xlab="Predicted Richness")
  abline(v=Richness[i],col="red")
  m <- paste("mean=",round(mean(s.sample[,i]),2))
  s <- paste("sd=",round(sd(s.sample[,i]),2))
  legend("topleft",c(m,s),bty="n")
}





######################################################################################################
# Extract estimates
mu.est <- mean(extract(m8_4_stanfit, pars = 'mu', inc_warmup = FALSE))
sigma.est <- mean(extract(m8_4_stanfit, pars = 'sigma', inc_warmup = FALSE))
residual.est <- (extract(m8_4_stanfit, pars = 'residual', inc_warmup = FALSE))

#mu_est
for(i in 1:n){
  mu_est[i] <- mean(extract(m8_4_stanfit, pars = 'mu[i]', inc_warmup = FALSE)) # doesn't work - can't index mu
}
mu.1 <- (extract(m8_4_stanfit, pars = 'mu[1]', inc_warmup = FALSE))


