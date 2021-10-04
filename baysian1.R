###
#BAYSIAN MODELING WITH STAN
###

setwd("C:/Users/Diogo Barros/Documents/Diogo/R data")
data <- read.csv("landscape.csv")
attach(data) #Richeness
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE) # for multiple cores
options(mc.cores = parallel::detectCores())

#Response Variable
y <- Richness

#Parameters
n <- 11
log_offset <- 1
X <- rep(1,11)

teste <- cbind(X,
          GreenUsable,
          F100,
          Habitats2,
          ShapeInd)

waics <- c(rep(0,5))

for(k in 1:4){
    p <- k+1 
    X <- cbind(X,teste[,k+1])
    
    stan_data <- list(y=y,
                  n=n,
                  p=p,
                  X=X)
  
    fit <- stan(file="pois_ln_regression_WAIC.stan",
              data=stan_data,iter=1000,
              control=list(adapt_delta=0.9), max_treedepth=15)
  


  
  ##
  #Extracting results
  ##
  
  results <-summary(fit)

  #Selecting the Lambda estimations 
  a<-extract(fit)
  lambdas <- a[["lambda"]]

  #Sampling Lambdas
  source("C:/Users/Diogo Barros/Documents/Diogo/The future/BEGC/R scripts/Task1/myFunctions.R")
  sampleSize<-1000
  sample <- sampler(lambdas,sampleSize,n)  #amostra 1 vez do output do STAN
  s.sample  <- s.gen(sample)               #Generate a matrix of S using Poiss(lambda)

  #Sampling WAIC
  waics[k] <- waic(fit)[[1]]

  #####
  #Plots
  #####

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
  for(i in 1:dim(a$lambda)[2]){
    hist(s.sample[,i],main=data$Name[i],xlim=c(0,12),xlab="Predicted Richness")
    abline(v=Richness[i],col="red")
    legend("topleft",c(m,s),bty="n")
  }

} #FIM DO LOOP ENORME


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

6
