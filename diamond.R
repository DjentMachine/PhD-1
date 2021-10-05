##
#Part - I Fitting Bayesian Regression to Diamond data
#
#Simple species Area relationship:
#s= c*A^Z <=> log(s) = log(c) + zlog(A) ==> REGRESSION
#S = number of ssp
#A = area
#c = Param1
#Z = Param2
#####



setwd("C:/Users/Diogo/Documents/Diogo/R data")
data <- read.csv("diamond.csv", sep=",")
attach(data) #Dependent variable
S = data[,1]
A = data[,2]


##
#Frequentist methods
##

# Simple SA relationship

lm1 = lm(log(S)~log(A))
cat(sprintf("C is %s \nZ is %f", exp(1)^lm1$coefficients[1], lm1$coefficients[2]))
summary(lm1)

##
# Adding more variables to the relationship: S = c * (1+ 0.027* L/1000) * (e^(-0.693d/1620) * A^z 
#log(S)=lob(C)+log(1+0.027*L/1000)+(-0.69d/1620)+zlog(A)
#
# 1 => Total relative number of species at sea level 
# 0.027 => Assumes that per each 1000 feet in altitude, 2.7% of speices are associated with altiture.
#         Its an increment to total species. The higher mountains, the more species
# L/1000 => altitude incremente per 1000 feet  
# d => distance to New Guinea (continent=?) in MILES
# e => nepper's number,  2.71^-0.693 =~ 0,5 . It takes 1620 miles for the S to halve (or decrese by a 2 fold factor)
##


gam


##
#Bayesian Methods
##
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE)  
options(mc.cores = parallel::detectCores()) # for multiple cores

#Parameters
N = length(S)
#K = length(colnames(data))-1
y = log(data[,1])
x <- log(data[,2])     
stan_data <- list(N=N,
                  y=y,
                  x=x)

fit <- stan(file="regressionSIMPLE.stan",data=stan_data,iter=4000,chain=4,
            cores = 4, 
            #adapt_delta=0.9,
            control=list(max_treedepth =20))

#saveRDS(fit, "normals.rds")  #GRAVA MODELO
##
#Extracting results
#TARGET: z= 0.22, logC=15.4
##

fit          #summary 
posterior = extract(fit) #Full posterior of parameters

#plotting the posterior distribution for the parameters
plot(y~x, pch=20, xlab="Log")
abline( mean(posterior$logC), mean(posterior$z))
abline(lm1, col = 2, lty = 2, lw = 3)
for (i in 1:500) {
  abline(posterior$logC[i], posterior$z[i], col = "gray", lty = 1)
}
plot(exp(y)~exp(x), pch=20, xlab="Log",
     xlim=c(0,max(exp(x[length(x)-1]))))


#Fit S=CA^Z
lines(exp(x),exp(mean(posterior$logC))*exp(x)^mean(posterior$z),col="red")

#Getting posterior histograms
for (i in 1:(length(names(posterior))-1)){ 
  hist(posterior[[i]], main = sprintf("Posterior of %s", names(posterior)[i]),
     xlab = sprintf("Value of %s", names(posterior)[i]))
}


layout(matrix(1:2,ncol=2,byrow=T))
plot(wing ~ age, pch = 20, xlab="Age (in days)", ylab="Wing length (in cm)", main="a")
for (i in 1:1000) {
  abline(posterior$alpha[i], posterior$beta[i], col = "gray", lty = 1)
}
abline(mean(posterior$alpha), mean(posterior$beta), col = 6, lw = 2)
points(age,wing,pch=20)
plot(wing ~ age, pch = 20, xlab="Age (in days)", ylab="Wing length (in cm)", main="b",yl
     for (i in 1:1000) {
       points(age,rnorm(length(age),posterior$alpha[i]+age*posterior$beta[i],posterior$sigma[]
     }
     points(age,wing,pch=20)


##
#READING PREVIOUS MODELS
##
fit = readRDS("uniforms.rds")

