########
#Models#
########


setwd("C:/Users/Diogo Barros/Documents/Diogo/R data")
data <- read.csv("landscape.csv")
attach(data)
detach(data)
library(lme4)
#library(nlme)
names(data)

summary(Richness)

m<- mean(Richness)
s <- var(Richness)
print(m)
print(s)
install.packages("AER")

####
#Diversity
###


#0 - No effects, null model
fit0<-lm(Richness ~ 1)
summary(fit0)
AIC(fit0)

#1 - All variables, no random effects
fit<-lm(Richness~ soundCentr + People + Area + GreenUsable + ShapeInd + NrFrags + Habitats + MonsantoDist + ClosestDist)
summary(fit)
AIC(fit)

#2 Luz
fit<-lm(Richness ~ LightMax)
summary(fit)
AIC(fit)

#3 Sound
fit<-lm(Richness ~ soundCentr + soundMarg)
summary(fit)
AIC(fit)

#4 Light and Sound
fit<-lm(Richness ~ soundCentr + soundMarg + LightMax)
summary(fit)
AIC(fit)

#4 Light and Sound
fit<-lm(Richness ~ soundCentr + soundMarg + LightMax)
summary(fit)
AIC(fit)

#5 People
fit<-lm(Richness ~ People)
summary(fit)
AIC(fit)

#6 Light and Sound and people
fit<-lm(Richness ~ soundCentr + soundMarg + LightMax + People)
summary(fit)
AIC(fit)

#7 Area + Perimeter + shape
fit<-lm(Richness ~ Area)
summary(fit)
AIC(fit)

#8 Distance
fit<-lm(Richness ~ ClosestDist + MonsantoDist + ClosestPark)
summary(fit)
AIC(fit)

#8 Area + Perimeter + shape
fit<-lm(Richness ~ Area + Perimeter + ShapeInd)
summary(fit)
AIC(fit)

#9 Best areas
fit<-lm(Richness ~ GreenUsable +  
          Area + Perimeter + ShapeInd)
summary(fit)
AIC(fit)



# - All variables, parks are random effects
fit<-lmer(Richness~ 1+ (1|Name)+ soundCentr + People + Area + GreenUsable + ShapeInd + NrFrags + Habitats + MonsantoDist + ClosestDist)
summary(fit)


AIC(fit0,fit1, fit2)


###
#Abundance
###
abun <- read.csv("hairIndex.csv")

#####MMUSCULUS#######
fit<-lm(abun$MM~ soundCentr + People + Area + GreenUsable + ShapeInd + NrFrags + Habitats + MonsantoDist + ClosestDist)
summary(fit)

#0 - No effects, null model
fit0<-lm(abun$MM ~ 1)
summary(fit0)
AIC(fit0)

#1 - All variables, no random effects
fit<-lm(abun$MM~ soundCentr + People + Area + GreenUsable + ShapeInd + NrFrags + Habitats + MonsantoDist + ClosestDist)
summary(fit)
AIC(fit)

#2 Luz
fit<-lm(abun$MM ~ LightMax)
summary(fit)
AIC(fit)

#3 Sound
fit<-lm(abun$MM ~ soundCentr + soundMarg)
summary(fit)
AIC(fit)

#4 Light and Sound
fit<-lm(abun$MM ~ soundCentr + soundMarg + LightMax)
summary(fit)
AIC(fit)


#5 People
fit<-lm(abun$MM ~ People)
summary(fit)
AIC(fit)

#6 Light and Sound and people
fit<-lm(abun$MM ~ soundCentr + soundMarg + LightMax + People)
summary(fit)
AIC(fit)

#7 Area + Perimeter + shape
fit<-lm(abun$MM ~ Area + Perimeter + ShapeInd)
summary(fit)
AIC(fit)

#8 Distance
fit<-lm(abun$MM ~ ClosestDist + MonsantoDist + ClosestPark)
summary(fit)
AIC(fit)


#9 Best areas
fit<-lm(abun$MM ~ GreenUsable +  
          Area + Perimeter + ShapeInd)
summary(fit)
AIC(fit)


#####CRUSSULA#######
fit<-lm(abun$CR~ soundCentr + People + Area + GreenUsable + ShapeInd + NrFrags + Habitats + MonsantoDist + ClosestDist)
summary(fit)

#0 - No effects, null model
fit0<-lm(abun$CR ~ 1)
summary(fit0)
AIC(fit0)

#1 - All variables, no random effects
fit<-lm(abun$CR~ soundCentr + People + Area + GreenUsable + ShapeInd + NrFrags + Habitats + MonsantoDist + ClosestDist)
summary(fit)
AIC(fit)

#2 Luz
fit<-lm(abun$CR ~ LightMax)
summary(fit)
AIC(fit)

#3 Sound
fit<-lm(abun$CR ~ soundCentr + soundMarg)
summary(fit)
AIC(fit)

#4 Light and Sound
fit<-lm(abun$CR ~ soundCentr + soundMarg + LightMax)
summary(fit)
AIC(fit)


#5 People
fit<-lm(abun$CR ~ People)
summary(fit)
AIC(fit)

#6 Light and Sound and people
fit<-lm(abun$CR ~ soundCentr + soundMarg + LightMax + People + Animals)
summary(fit)
AIC(fit)

#7 Area + Perimeter + shape
fit<-lm(abun$CR ~ Area + Perimeter + ShapeInd)
summary(fit)
AIC(fit)

#8 Distance
fit<-lm(abun$CR ~ ClosestDist + MonsantoDist + ClosestPark)
summary(fit)
AIC(fit)


#9 Best areas
fit<-lm(abun$CR ~ GreenUsable +  
          Area + Perimeter + ShapeInd)
summary(fit)
AIC(fit)








###
#CCA
###

require(ggplot2)
require(GGally)
require(CCA)

###
#CCA 2
###

library(vegan)
data(varespec) # matriz com dados vegetacionais, com  24 linhas (parcelas) e 44 colunas (espécies).
data(varechem) # matriz com dados ambientais, com 24 linhas (parcelas) e 14 colunas (variáveis), referente às características edáficas de cada parcela
varespec # visualizar matriz
varechem # visualizar matriz

vare.cca <- cca(varespec ~ Baresoil+Humdepth+pH+N+P+K+Ca+Mg+S+Al+Fe, data=varechem)
vare.cca
plot(vare.cca)
summary(vare.cca)


###
#Landscape 
###

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



