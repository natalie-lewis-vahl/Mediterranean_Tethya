library(dplyr)
library(ggplot2)
taxa<-read.table("./Data/18s/taxa_fixed18s.csv",header=TRUE,sep=",")
otus<-read.table("./Data/18s/all.otutab.csv",sep="\t")
names(taxa)
names(otus)
reorder_id<-match(taxa$sequence_identifier,otus$V1)
otus<-otus[reorder_id,]

merged<-merge(taxa,otus,by.x="sequence_identifier",by.y="V1")
str(merged)
dim(merged)
xmerged<-merged[!(merged$Domain=="Bacteria"),]
count<-apply(xmerged[-c(1:11)],1,sum)
xmerged<-cbind(xmerged,count)
#Uncomment to also take out OTU 1 and 2
xfiltered<-xmerged[!(xmerged$sequence_identifier=="OTU_1"|xmerged$sequence_identifier=="OTU_2"),]

#Remove OTUS with less than 5 reads
filtered<-xfiltered[xfiltered$count>5,]
filtered[filtered==""]<-"Unclassified"

plot1data<-filtered%>%
  count(Phylum)

data<-plot1data$n
names(data)<-plot1data$Phylum
barplot(data)
#The large amount of otus belonging to unclassified makes other phyla harder
#to differentiate - could represent same phyla or lots of different ones
#remove 9215 unknown OTUs 
plot2data<-plot1data[!(plot1data$Phylum=="Unclassified"),]
data2<-plot2data$n
names(data2)<-plot2data$Phylum
barplot(data2)
#delete row with unclassified phyla
filtered<-filtered[!(filtered$Phylum=="Unclassified"),]
counts<-table(filtered$Class,filtered$Phylum)
jpeg("./Figures/18splots/overall18s_taxa_distribution.jpeg",units="in", width=15, height=5, res=300)
barplot(counts,xlab="Phyla",ylab="Number of OTUs",cex.names= 0.5,ylim=c(0,25))
dev.off()
??barplot

unique(filtered$Class[filtered$Phylum=="Annelida"])
length(filtered$Class=="Polychaeta"~filtered$Phylum=="Annelida")


head(filtered)
#Tau V2 till v12 col 12 to 22
#Tci v13 till v34 col 23 to 44
#Tme v35 till v45 col 45 till 55

unique(filtered$Class[filtered$Phylum=="Ciliophora"])

length(which(filtered$Class=="Heterotrichea"))

unique(filtered$Order[filtered$Class=="Spirotrichea"])       

length(which(filtered$Order=="Sporadotrichida"))

length(which(!filtered$Phylum=="Unclassified"))

