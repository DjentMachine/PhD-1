########
#Landscape Models - 2
########

setwd("C:/Users/Diogo Barros/Documents/Diogo/R data")
data <- read.csv("landscape.csv")
attach(data)
detach(data)
library(lme4)
#library(nlme)
names(data)

###Richness
fit1<-glm( Richness~ GreenArea, family=poisson)
fit2<-glm( Richness~ GreenArea + F100, family = poisson)
fit3<-glm( Richness~ GreenArea + F100 + Habitats)
fit4<-glm( Richness~ GreenArea + F100 + ShapeInd + Habitats2)
fit5<-glm( Richness~ GreenArea + F500, family = poisson)
fit6<-glm( Richness~ GreenArea + F500 + Habitats)
fit7<-glm( Richness~ GreenArea + F500 + ShapeInd + Habitats2)
fit8<-glm( Richness~ GreenArea + F1000, family = poisson)
fit9<-glm( Richness~ GreenArea + F1000 + Habitats)
fit10<-glm( Richness~ GreenArea + F1000 + ShapeInd + Habitats2)

aics<-data.frame(paste("fit",1:10,sep=""),c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,
                                            fit5$aic,fit6$aic,fit7$aic,fit8$aic,
                                            fit9$aic,fit10$aic),row.names=NULL)
colnames(aics)<-c("model","AIC")
aics<-aics[order(-aics$AIC),]
for(i in 1:dim(aics)[1]){
  aics$diff[i]<-aics$AIC[1]-aics$AIC[i]}
aics$wi<-2.71828182845904523536^(-0.5*aics$diff)
aics$aic.weights<-aics$wi/sum(aics$wi)
write.csv(aics, "AICS1.csv")


### TRUE Richness
fit1<-glm( RichMax~ GreenArea, family=poisson)
fit2<-glm( RichMax~ GreenArea + F100, family = poisson)
fit3<-glm( RichMax~ GreenArea + F100 + Habitats)
fit4<-glm( RichMax~ GreenArea + F100 + ShapeInd + Habitats2)
fit5<-glm( RichMax~ GreenArea + F500, family = poisson)
fit6<-glm( RichMax~ GreenArea + F500 + Habitats)
fit7<-glm( RichMax~ GreenArea + F500 + ShapeInd + Habitats2)
fit8<-glm( RichMax~ GreenArea + F1000, family = poisson)
fit9<-glm( RichMax~ GreenArea + F1000 + Habitats)
fit10<-glm( RichMax~ GreenArea + F1000 + ShapeInd + Habitats2)

aics<-data.frame(paste("fit",1:10,sep=""),c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,
                                            fit5$aic,fit6$aic,fit7$aic,fit8$aic,
                                            fit9$aic,fit10$aic),row.names=NULL)
colnames(aics)<-c("model","AIC")
aics<-aics[order(-aics$AIC),]
for(i in 1:dim(aics)[1]){
  aics$diff[i]<-aics$AIC[1]-aics$AIC[i]}
aics$wi<-2.71828182845904523536^(-0.5*aics$diff)
aics$aic.weights<-aics$wi/sum(aics$wi)
write.csv(aics, "AICS15.csv")



###Saturation
fit1<-glm( Saturation~ GreenArea)
fit2<-glm( Saturation~ GreenArea + F100)
fit3<-glm( Saturation~ GreenArea + F100 + Habitats)
fit4<-glm( Saturation~ GreenArea + F100 + ShapeInd + Habitats2)
fit5<-glm( Saturation~ GreenArea + F500)
fit6<-glm( Saturation~ GreenArea + F500 + Habitats)
fit7<-glm( Saturation~ GreenArea + F500 + ShapeInd + Habitats2)
fit8<-glm( Saturation~ GreenArea + F1000)
fit9<-glm( Saturation~ GreenArea + F1000 + Habitats)
fit10<-glm( Saturation~ GreenArea + F1000 + ShapeInd + Habitats2)

aics<-data.frame(paste("fit",1:10,sep=""),c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,
                                            fit5$aic,fit6$aic,fit7$aic,fit8$aic,
                                            fit9$aic,fit10$aic),row.names=NULL)

colnames(aics)<-c("model","AIC")
aics<-aics[order(-aics$AIC),]
for(i in 1:dim(aics)[1]){
  aics$diff[i]<-aics$AIC[1]-aics$AIC[i]}
aics$wi<-2.71828182845904523536^(-0.5*aics$diff)
aics$aic.weights<-aics$wi/sum(aics$wi)
write.csv(aics, "AICS2.csv")

###Full Saturation
fit1<-glm( SaturationFull~ GreenArea)
fit2<-glm( SaturationFull~ GreenArea + F100)
fit3<-glm( SaturationFull~ GreenArea + F100 + Habitats)
fit4<-glm( SaturationFull~ GreenArea + F100 + ShapeInd + Habitats2)
fit5<-glm( SaturationFull~ GreenArea + F500)
fit6<-glm( SaturationFull~ GreenArea + F500 + Habitats)
fit7<-glm( SaturationFull~ GreenArea + F500 + ShapeInd + Habitats2)
fit8<-glm( SaturationFull~ GreenArea + F1000)
fit9<-glm( SaturationFull~ GreenArea + F1000 + Habitats)
fit10<-glm( SaturationFull~ GreenArea + F1000 + ShapeInd + Habitats2)

aics<-data.frame(paste("fit",1:10,sep=""),c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,
                                            fit5$aic,fit6$aic,fit7$aic,fit8$aic,
                                            fit9$aic,fit10$aic),row.names=NULL)

colnames(aics)<-c("model","AIC")
aics<-aics[order(-aics$AIC),]
for(i in 1:dim(aics)[1]){
  aics$diff[i]<-aics$AIC[1]-aics$AIC[i]}
aics$wi<-2.71828182845904523536^(-0.5*aics$diff)
aics$aic.weights<-aics$wi/sum(aics$wi)
write.csv(aics, "AICS25.csv")

###Shannon

fit1<-glm( Shannon ~ GreenArea)
fit2<-glm( Shannon~ GreenArea + F100)
fit3<-glm( Shannon~ GreenArea + F100 + Habitats)
fit4<-glm( Shannon~ GreenArea + F100 + ShapeInd + Habitats2)
fit5<-glm( Shannon~ GreenArea + F500)
fit6<-glm( Shannon~ GreenArea + F500 + Habitats)
fit7<-glm( Shannon~ GreenArea + F500 + ShapeInd + Habitats2)
fit8<-glm( Shannon~ GreenArea + F1000)
fit9<-glm( Shannon~ GreenArea + F1000 + Habitats)
fit10<-glm( Shannon~ GreenArea + F1000 + ShapeInd + Habitats2)

aics<-data.frame(paste("fit",1:10,sep=""),c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,
                                            fit5$aic,fit6$aic,fit7$aic,fit8$aic,
                                            fit9$aic,fit10$aic),row.names=NULL)

colnames(aics)<-c("model","AIC")
aics<-aics[order(-aics$AIC),]
for(i in 1:dim(aics)[1]){
  aics$diff[i]<-aics$AIC[1]-aics$AIC[i]}
aics$wi<-2.71828182845904523536^(-0.5*aics$diff)
aics$aic.weights<-aics$wi/sum(aics$wi)
write.csv(aics, "AICS3.csv")

###Simpson

fit1<-glm( Simpson ~ GreenArea)
fit2<-glm( Simpson~ GreenArea + F100)
fit3<-glm( Simpson~ GreenArea + F100 + Habitats)
fit4<-glm( Simpson~ GreenArea + F100 + ShapeInd + Habitats2)
fit5<-glm( Simpson~ GreenArea + F500)
fit6<-glm( Simpson~ GreenArea + F500 + Habitats)
fit7<-glm( Simpson~ GreenArea + F500 + ShapeInd + Habitats2)
fit8<-glm( Simpson~ GreenArea + F1000)
fit9<-glm( Simpson~ GreenArea + F1000 + Habitats)
fit10<-glm( Simpson~ GreenArea + F1000 + ShapeInd + Habitats2)

aics<-data.frame(paste("fit",1:10,sep=""),c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,
                                            fit5$aic,fit6$aic,fit7$aic,fit8$aic,
                                            fit9$aic,fit10$aic),row.names=NULL)

colnames(aics)<-c("model","AIC")
aics<-aics[order(-aics$AIC),]
for(i in 1:dim(aics)[1]){
  aics$diff[i]<-aics$AIC[1]-aics$AIC[i]}
aics$wi<-2.71828182845904523536^(-0.5*aics$diff)
aics$aic.weights<-aics$wi/sum(aics$wi)
write.csv(aics, "AICS4.csv")