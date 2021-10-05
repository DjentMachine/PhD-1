###########
#Landscape#
###########


setwd("C:/Users/Splinterq/Documents/Diogo/R data")
data <- read.csv("landscape.csv", sep=";")
attach(data)
detach(data)



##
#Cluster analyses
##

#DendrogramSound

dataSound <- data[c(1,24)]
d <- dist(as.matrix(dataSound))   #dist matrix
hc <- hclust(d)
summary(hc)
plot(hc,
     labels = nomeCurto 
)

#DendrogramArea

dataSound <- data[c(1,5)]
d <- dist(as.matrix(dataSound))   #dist matrix
hc <- hclust(d)
summary(hc)
plot(hc,
     labels = nomeCurto 
)

#DendrogramGreen

dataSound <- data[c(1,8)]
d <- dist(as.matrix(dataSound))   #dist matrix
hc <- hclust(d)
summary(hc)
plot(hc,
     labels = nomeCurto 
)

#DendrogramAnthro

dataSound <- data[c(1,26,27,28)]
d <- dist(as.matrix(dataSound))   #dist matrix
hc <- hclust(d)
summary(hc)
plot(hc,
     labels = nomeCurto 
)


#Dendrogram

d <- dist(as.matrix(data))   #dist matrix
hc <- hclust(d)
summary(hc)
plot(hc,
     labels = nomeCurto 
     )

# k-means

install.packages('rattle')
data(wine, package='rattle')
head(wine)

set.seed(123456789) ## to fix the random starting clusters
greenArea <- kmeans(food[,c("", "")], centers=3, nstart=10)
grpMeat
###


library(reshape2)
pivot <-melt (data)

library(cluster)
D=daisy(pivot, metric='gower')

H.fit <- hclust(D, method="ward.D")

plot(H.fit) # display dendrogram

groups <- cutree(H.fit, k=4) # cut tree into 4 clusters

# draw dendogram with red borders around the 4 clusters
rect.hclust(H.fit, k=4, border="red") 

clusplot(pivot, groups, color=TRUE, shade=TRUE,
         labels=2, lines=0, main= 'Customer segments')

####
#MORE CLUSTER ANALYSES - THE PCA
###

#excepções <- names(data) %in% c("nomeCurto", "ShrubArea", "NakedArea", "HerbsArea", "Light",
#                                "NrFragsShrub", "NrFragsNaked", "NrFragsHerbs") 
#pcaData<- data[!excepções]
#logData <-cbind((log(pcaData[3:8])),(log(pcaData[10:15])))
pcaData <- names(data) %in% c("ShapeInd", "Area", "Perimeter", "GreenUsable",
                              "NrFragsGreenUs", "MonsantoDist", "Habitats2",
                              "soundCentr","soundMarg","People","Vehicles","Animals",
                              "RCR","RMM","RMS","RRN")
pcaData <-data[pcaData]
myPca <- prcomp(pcaData,center = TRUE,scale. = TRUE)
plot(myPca, type ="lines")
plot( myPca[,1], myPca[,2] )
a<-summary(myPca)
print(myPca)
biplot(myPca, xlabs=(data$nomeCurto), xlab= "PC1\nVariation exlpained: 39.87%",
       ylab="PC2\nVariation exlpained: 24.43% ")

cor.test(Richness,Habitats2,pairwise=TRUE)

#GGPLOT STUFF
library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

