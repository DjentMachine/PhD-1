###
#BAYSIAN MODELING WITH STAN
###

setwd("C:/Users/Diogo Barros/Documents/Diogo/R data")
data <- read.csv("landscape.csv")
#attach(data) #Richeness
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE) # for multiple cores
options(mc.cores = parallel::detectCores())

#Parameters
n <- 11
p <- 1
X <- cbind(rep(1,11))
#          LogGreen,
#           F500)
y <- Richness

stan_data <- list(y=y,
                  n=n,
                  p=p,
                  X=X,
                  log_offset=log(1))

fit2 <- stan(file="pois_ln_regression.stan",data=stan_data,iter=3000)

#####
#Extracting results
#####

print(fit)

results <-summary(fit)
write.csv(results[[1]], "C:/Users/Diogo Barros/Desktop/NewTask1/STAN/area+F100/summaryF1000.csv", col.names = TRUE)
plot(fit, pars =names(fit))
#traceplot(fit, pars =names(fit))
#names <- names(fit)

mean <- c(rep(0,11))
sd <- c(rep(0,11))
a<-extract(fit)
par(mfrow=c(3,4))
for(i in 1:dim(a$lambda)[2]){
  hist(a$lambda[,i],main=i,xlim=c(0,12))
  mean[i]=round(mean(a$lambda[,i]),2)
  m <- paste("mean=",mean[i])
  sd[i]= round(sd(a$lambda[,i]),2)
  s <- paste("sd=",sd[i])
  legend(x=2,y=500,c(m,s),bty="n")
}
write.csv(cbind(mean,sd), "C:/Users/Diogo Barros/Desktop/NewTask1/STAN/area+F1000/meanF1000.csv", col.names = TRUE)

all<-extract(fit)
boxplot(mean)

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

