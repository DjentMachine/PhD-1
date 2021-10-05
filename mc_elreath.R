# load the distance matrix
library(rethinking)
library(bayesplot)

data(islandsDistMatrix)
# display (measured in thousands of km)
Dmat <- islandsDistMatrix
colnames(Dmat) <- c("Ml","Ti","SC","Ya","Fi","Tr","Ch","Mn","To","Ha")
round(Dmat,1)

data(Kline2) # load the ordinary data, now with coordinates
d <- Kline2
d$society <- 1:10 # index observations

14.38
dat_list <- list(
  T = d$total_tools,
  P = d$population,
  society = d$society,
  Dmat=islandsDistMatrix )
m14.7 <- ulam(
  alist(
    T ~ dpois(lambda),
    lambda <- (a*P^b/g)*exp(k[society]),
    vector[10]:k ~ multi_normal( 0 , SIGMA ),
    matrix[10,10]:SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
    c(a,b,g) ~ dexp( 1 ),
    etasq ~ dexp( 2 ),
    rhosq ~ dexp( 0.5 )
  ), data=dat_list , chains=4 , cores=4 , iter=2000 )

post <- extract.samples(m14.7)
mcmc_areas(post2[,1:11])

post2=as.matrix(m14.7@stanfit)

canary_title1 <- ggtitle("Posterior distributions for Mc Elreaths' Island models",
                         "with medians and 95% intervals")
mcmc_areas(post2,
           pars = c("a", "b", #"b2", "b3",
                     "k[1]", "k[2]", "k[3]", "k[4]", "k[5]", "k[6]", "k[7]",
           "etasq","rhosq"),
           ,prob = 0.95) + canary_title1

