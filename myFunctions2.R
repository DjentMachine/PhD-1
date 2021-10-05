###
# Helper functions for Macaronesian Species Richness
###

##
# radians to degrees and vice versa converters
##
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

##
#Getting a distance matrix out of 2 collomns of WGS86 coordinates:
#Latitude first col, longitude 2nd col
##
distMatrix <- function(coords){
  dist <- matrix(ncol=length(coords[,1]),nrow=length(coords[,1]))
  coords=deg2rad(coords)
  for(i in 1:length(dist[,1])){
    for(j in 1:length(dist[,1])){
      dlon = coords[j,2] - coords[i,2]
      dlat = coords[j,1] - coords[i,1] 
      a = (sin(dlat/2))^2 + cos(coords[i,1]) * cos(coords[j,1]) * (sin(dlon/2))^2
      c = 2 * atan2( sqrt(a), sqrt(1-a) )
      dist[i,j] = 6378100*c #(where R is the radius of the Earth)
    }
  }
  return(dist)
}

##
#Function to turn a dataframe containing numbers to a matrix
##
df2matrix <- function(df){
  m = matrix(0,nrow=length(df[,1]),ncol=length(df))
  for(i in 1:length(m[1,])){
    for(j in 1:length(m[,1])){
      m[j,i]= as.numeric(as.character(df[j,i]))
    }
  }
  return(m)  
}

##
#Extracting posterior from model
##
getPost <- function(model){
  post = extract(model)
  plot(density(post$z), xlab="alpha (the intercept)", main = "a")
  plot(density(post$logC), xlab="beta (the slope)", main = "b")
  plot(density(post$sigma), xlab="sigma", main = "c")
  return(post)
}

##
#Getting graphs from model
##
modelGraphs <- function(modelFit, y, x){
  posterior = extract(modelFit)
  layout(matrix(1:6,ncol=3,byrow=T))
  plot(y ~ x, pch = 20, ylim=c(min(y)*0.90,max(y)*1.1), ylab="Log S (richness)",
       xlab="Log Area (no units)", main="a")
  for (i in 1:1000) {
    abline(posterior$logC[i], posterior$z[i],  col = "gray", lty = 1)
  }
  abline(mean(posterior$logC),mean(posterior$z),  col = 6, lw = 2)
  points(A,S,pch=20)
  plot(S ~ A, pch = 20, ylim=c(min(S)*0.90,max(S)*1.1), ylab="Log S (richness)",
       xlab="Log Area (no units)", main="b")
  for (i in 1:1000) {
    points(A,rnorm(length(S),posterior$logC[i]+A*posterior$z[i],posterior$sigma[i]),
           col = "gray", lty = 1)
  }
  points(A,S,pch=20)
}
