#species rarefaction curves for 18s data

CountsPath<-"./Data/18s/all.otutab.csv"
#tbc bactLoadCorrectedCountsPath<-"./Data/all.otutab_Corrected_valuesOnly.csv"
otuTaxonomyPath<-"./Data/18s/taxa_fixed18s.csv"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep=",")
countsDF<-read.csv(CountsPath, sep="\t")
head(otuTaxonomy)
data<-cbind(countsDF[order(countsDF$X.OTU.ID),],otuTaxonomy[order(otuTaxonomy$sequence_identifier),])
head(data)
unique(data$Phylum)
#Filtering
countSum<-apply(countsDF[2:45],1,sum) #sum the rows
xdata<-cbind(data,countSum)
xdata<-xdata[xdata$countSum >=5,]
xdata<-xdata%>%filter(Domain!="Bacteria")
unique(xdata$Phylum)
head(xdata)
#remove otu 1 and 2 because of coamplification
which(xdata$X.OTU.ID=="OTU_1")
which(xdata$X.OTU.ID=="OTU_2")
xdata<-xdata[-c(1,85),]
#take rows for species T. aurantium and transpose for analysis
data<-t(as.data.frame(xdata[,c(2:12,45)]))
head(data)
library(vegan)
S <- specnumber(data)
raremax<-min(rowSums(data))
Srare <- rarefy(data, raremax)
par(mfrow=c(1,1))
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(data, step =1, sample = raremax, col = "blue", cex = 0.4)
#For T meloni
data_me<-t(as.data.frame(xdata[,c(13:23,44,45)]))
S_me <- specnumber(data_me)
raremax_me<-min(rowSums(data_me))
Srare_me <- rarefy(data_me, raremax_me)
par(mfrow=c(1,1))
plot(S_me, Srare_me, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(data_me, step =1, sample = raremax_me, col = "blue", cex = 0.4)
#for T citrina
data_ci<-t(as.data.frame(xdata[,24:42]))
S_ci <- specnumber(data_ci)
raremax_ci<-min(rowSums(data_ci))
Srare_ci <- rarefy(data_ci, raremax_ci)
par(mfrow=c(1,1))
plot(S_ci, Srare_ci, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(data_ci, step =1, sample = raremax_ci, col = "blue", cex = 0.4)


png("Figures/18splots/Rarefaction_curve_18s.png",width=20,height=30,units="cm",res=300)
par(mfrow=c(3,2))
#tau
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",main=substitute(paste(italic("Tethya aurantium"))))
abline(0, 1)
rarecurve(data, step =1, sample = raremax, col = "blue", cex = 0.4)
#tme
plot(S_me, Srare_me, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",main=substitute(paste(italic("Tethya meloni"))))
abline(0, 1)
rarecurve(data_me, step =1, sample = raremax_me, col = "blue", cex = 0.4)
#tci
plot(S_ci, Srare_ci, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",main=substitute(paste(italic("Tethya citrina"))))
abline(0, 1)
rarecurve(data_ci, step =1, sample = raremax_ci, col = "blue", cex = 0.4)

dev.off()
