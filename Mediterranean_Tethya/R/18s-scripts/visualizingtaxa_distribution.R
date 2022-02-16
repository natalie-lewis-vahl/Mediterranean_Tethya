library(dplyr)
library(ggplot2)
library(grDevices)
taxa<-read.table("./Data/18s/taxa_fixed18s.csv",header=TRUE,sep=";")
otus<-read.table("./Data/18s/all.otutab.csv",sep=";")
names(taxa)
names(otus)
reorder_id<-match(taxa$sequence_identifier,otus$V1)
otus<-otus[reorder_id,]

merged<-merge(taxa,otus,by.x="sequence_identifier",by.y="V1")
str(merged)
dim(merged)
amerged<-merged[!(merged$Domain=="Bacteria"),]
xmerged<-amerged[-c(1:11)]
#Convert columns being read as characcters into numeric
df2<-lapply(xmerged,as.numeric)
count<-apply(as.data.frame(df2),1,sum)
xmerged<-cbind(amerged,count)
#Uncomment to also take out OTU 1 and 2
xfiltered<-xmerged[!(xmerged$sequence_identifier=="OTU_4076"|xmerged$sequence_identifier=="OTU_4077"),]

#Remove OTUS with less than 5 reads
filtered<-xfiltered[xfiltered$count>=5,]
filtered[filtered==""]<-"Unclassified"

plot1data<-filtered%>%
  dplyr::count(Phylum)
View(plot1data2)
plot1data2<-filtered%>%
  dplyr::count(Phylum,Class)
#The large amount of otus belonging to unclassified makes other phyla harder
#to differentiate - could represent same phyla or lots of different ones
#remove 9215 unknown OTUs 

plot1data<-plot1data[!(plot1data$Phylum=="Unclassified"),]

#delete row with unclassified phyla
filtered<-filtered[!(filtered$Phylum=="Unclassified"),]
data<-plot1data$n
names(data)<-plot1data$Phylum
barplot(data)
View(data)
counts<-table(filtered$Class,filtered$Phylum)
View(counts)
png("./Figures/18splots/overall18s_taxa_distribution.png",units="in", width=15, height=5, res=300)
barplot(counts,xlab="Phyla",ylab="Number of OTUs",cex.names= 0.5,ylim=c(0,25))
dev.off()
##for class
counts2<-table(filtered$Class)
png("./Figures/18splots/overall18s_class_distribution.png",units="in", width=15, height=5, res=300)
barplot(counts2,xlab="Class",ylab="Number of OTUs",cex.names= 0.5)
dev.off()
View(counts2)
dim(counts2)
#using ggplot
#Arrange table for this and order to have nice plot
head(plot1data)
datax<-plot1data
datax<-datax %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))
unique(datax$Phylum)

################
###########
##Set colours
Unclassified<-"slategray"
Ciliophora<-"green4"
Platyhelminthes<-"wheat3"
Rhodophyta<-"navyblue"
Porifera<-"yellow2"
Annelida<-"violetred4"
Nematoda<-"violetred2"
Cnidaria<-"violetred"
Euglenoza<-"violet"
Bacillariophyta<-"turquoise2"
Magnoliophyta<-"darkgreen"
Mollusca<-"tomato3"                 
Bryoza<-"tomato1"   
Echinodermata<-"purple"   
Kinorhyncha<-"thistle4"                    
Chlorophyta<-"tan4"               
Arthropoda<-"tan3"
Dinoflagellata<-"steelblue4"
Brachiopoda<-"steelblue2"
Chordata<-"yellowgreen"
####
head(datax)
datax[datax$phylum==""]<-"Unclassified"
png("./Figures/18splots/overall18s_taxa_distributionunclass.png",units="in", width=18, height=7, res=300)
ggplot(datax,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+scale_y_continuous(limits=c(0,250),expand=c(0,0),minor_breaks = seq(0 , 250, 10), breaks = seq(0, 250, 50))+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_manual(values = c(Unclassified,Ciliophora,Porifera,Bacillariophyta,Annelida,Cnidaria,Arthropoda,
                               Dinoflagellata,Rhodophyta,Chlorophyta,Chordata,Mollusca,Brachiopoda,Bryoza,
                               Echinodermata,Euglenoza,Kinorhyncha,Magnoliophyta,Nematoda,Platyhelminthes)) +
  labs(x="Phylum",y="Number of OTUs")+ggtitle("All 3 tethya species combined")+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9),legend.position="none",plot.title=element_text(hjust=0.5))+
  theme(text=element_text(size=12),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ) )
dev.off()

#species seperately
#Tau V2 till v12+v45 (col 9-19,52)
dim(filtered)
filteredxx<-select(filtered,-size,-sequence_score,-identity)
yfiltered<-filteredxx
yfiltered[,9:53]<-lapply(yfiltered[,9:53],as.numeric)
filteredtau<-yfiltered[,-c(20:51)]
dim(filteredtau)
taux<-filteredtau[!apply(filteredtau[,9:20]==0,1,all),]
dim(taux)
str(filteredtau)
rowSums(filteredtau[,9:20])
#Tme v13 till v23 +v43tillv44 (col 20-30,50,51)
filteredtme<-filtered[,-c(9:19,31:49,52)]
head(filteredtme)
tmex<-filteredtme[!apply(filteredtme[,9:21]==0,1,all),]

#Tci v35 till v42 (col 42-49)
filteredtci<-filtered[-c(9:30,50:52)]
tcix<-filteredtci[!apply(filteredtci[,9:27]==0,1,all),]
#scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
head(taux)
plot2data<-taux%>%
  dplyr::count(Phylum)
datatau<-plot2data %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))
##########
head(datatau)
tiff("./Figures/18splots/tau18s_taxa_distribution2.tiff",units="in", width=6, height=5, res=300)
ggplot(datatau,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c(Ciliophora,Bacillariophyta,Porifera,Annelida,Cnidaria,Arthropoda,
                               Rhodophyta,Chlorophyta,Dinoflagellata,Chordata,Mollusca,Brachiopoda,Bryoza,
                               Echinodermata,Euglenoza,Kinorhyncha,Nematoda,Platyhelminthes)) +
  labs(x="Phylum",y="Number of OTUs")+
  ggtitle(expression(italic("Tethya aurantium")))+scale_y_continuous(limits = c(0,25),expand=c(0,0),minor_breaks = seq(0 , 25, 1), breaks = seq(0, 100, 5))+
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=0.8,size=6),legend.position="none")+
  theme(text=element_text(size=9),axis.title.x=element_text(vjust=3),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ))
dev.off()

####Tme
plot3data<-tmex%>%
  dplyr::count(Phylum)
datatme<-plot3data %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))
##########
tiff("./Figures/18splots/tme18s_taxa_distribution2.tiff",units="in", width=6, height=5, res=300)
ggplot(datatme,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c(Ciliophora,Porifera,Bacillariophyta,Annelida,Arthropoda,Dinoflagellata,
                               Rhodophyta,Chlorophyta,Chordata,Cnidaria,Mollusca,Brachiopoda,Bryoza,
                               Echinodermata,Euglenoza,Kinorhyncha,Magnoliophyta,Nematoda)) +
  labs(x="Phylum",y="Number of OTUs")+
  ggtitle(expression(italic("Tethya meloni")))+scale_y_continuous(limits = c(0,25),expand=c(0,0),minor_breaks = seq(0 , 25, 1), breaks = seq(0, 100, 5))+
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=0.8,size=6),legend.position="none")+
  theme(text=element_text(size=9),axis.title.x=element_text(vjust=3),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ))

dev.off()

####Tci
plot4data<-tcix%>%
  dplyr::count(Phylum)
datatci<-plot4data %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))

##########
tiff("./Figures/18splots/tci18s_taxa_distribution2.tiff",units="in", width=6, height=5, res=300)
ggplot(datatci,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c(Ciliophora,Porifera,Annelida,Bacillariophyta,Cnidaria,Arthropoda,
                               Chordata,Dinoflagellata,Mollusca,Bryoza,Chlorophyta,
                               Echinodermata,Euglenoza,Kinorhyncha,Rhodophyta)) +
  labs(x="Phylum",y="Number of OTUs")+
  ggtitle(expression(italic("Tethya citrina")))+scale_y_continuous(limits = c(0,25),expand=c(0,0),minor_breaks = seq(0 , 25, 1), breaks = seq(0, 100, 5))+
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=0.8,size=6),legend.position="none")+
  theme(text=element_text(size=9),axis.title.x=element_text(vjust=3),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ))


dev.off()


