###
#Exploratoory stuff
###

setwd("C:/Users/Diogo/Documents/Diogo/R data")
data <- read.csv("landscape.csv", sep=";")
attach(data) #Dependent variable

rel <- data[52:55]

for(i in 1:length(rel)){
  v <- rel[,i]
  print(names(rel[i]))
  print(mean(v))
  print(var(v))
  print("")
}



