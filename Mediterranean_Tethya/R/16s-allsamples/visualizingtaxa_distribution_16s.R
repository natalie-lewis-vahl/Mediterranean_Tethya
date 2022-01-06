library(dplyr)
library(ggplot2)
taxa<-read.table("./Data/16s_allsamples_taxa.csv",header=TRUE,sep=";")
otus<-read.table("./Data/16s_allsamples_otu.csv",sep="\t")
names(taxa)
names(otus)
reorder_id<-match(taxa$sequence_identifier,otus$V1)
otus<-otus[reorder_id,]

merged<-merge(taxa,otus,by.x="sequence_identifier",by.y="V1")
str(merged)
dim(merged)
count<-apply(merged[-c(1:8)],1,sum)
merged<-cbind(merged,count)

#Remove OTUS with less than 50 reads
filtered<-merged[merged$count>50,]
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
jpeg("./Figures/AS_16splots/overall16s_taxa_distribution.jpeg",units="in", width=15, height=5, res=300)
barplot(counts,xlab="Phyla",ylab="Number of OTUs",cex.names= 0.5,ylim=c(0,160))
dev.off()
??barplot

head(filtered)
#Tau V2 till v12 col 12 to 22
#Tci v13 till v34 col 23 to 44
#Tme v35 till v45 col 45 till 55

unique(filtered$Class[filtered$Phylum=="Proteobacteria"])

length(which(filtered$Phylum=="Proteobacteria"))

