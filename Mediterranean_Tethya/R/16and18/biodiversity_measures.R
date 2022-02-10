#Alpha, beta and gamma
library(matrixStats)
library(dplyr)
library(vegan)
library(ggplot2)
library(cowplot)

taotus<-read.csv("./Data/16and18s_otu.csv", sep=";")
otus<-arrange(taotus, X.OTU.ID)
#BInd taxonomy to remove bacterial OTUs from 18s dataset 
taxa<-read.csv("./Data/16ans18s_taxa.csv",sep=";")
head(otus)
taxa<-arrange(taxa,sequence_identifier)
otu_and_taxa<-bind_cols(otus,taxa)
otu_and_taxa<-otu_and_taxa%>%filter(!(Domain=="Bacteria" & dataset=="18s"))

#and remove two otus fromcoamplified sponges for 18s
#Take out OTU 1 for sponge amplification
head(otu_and_taxa)
otu_and_taxa<-otu_and_taxa[-c(47)]
countSum<-apply(otu_and_taxa[3:46],1,sum)
otus_and_taxa<-cbind(otu_and_taxa,countSum)
coamplif_sponge<-otus_and_taxa[otus_and_taxa$dataset=="18s",]%>%
  arrange(desc(countSum))
#
coamplif_sponge$sequence_identifier[1]#OTU_4076
otus_and_taxa<-otus_and_taxa[!(otus_and_taxa$X.OTU.ID=="OTU_4076"),]
coamplif_sponge$sequence_identifier[2]#OTU_4077
otus_and_taxa<-otus_and_taxa[!(otus_and_taxa$X.OTU.ID=="OTU_4077"),]

#filtering: should be done seperate for both 16 and 18s data set
#Keep rows that have more than 50 counds and are from 16s data set OR have more than 5 counts and belong to 18s ddataset
otus_and_taxab<-otus_and_taxa[otus_and_taxa$countSum >= 50 & otus_and_taxa$dataset=="16s"|otus_and_taxa$countSum >= 5 & otus_and_taxa$dataset=="18s",]
head(otus_and_taxab)
data_all<-otus_and_taxab
library(writexl)

write_xlsx(data_all,"Figures/data_full_filtered.xlsx")

#For absence or presence data
data_all<-as.data.frame(data_all[,-c(1:2,47:55)]>0)
data_all<-cbind(data_all,otus_and_taxab[1])
#turns into logical value, true or false if present or not

#just 16s
data_16<-data_all[data_all$dataset=="16s",]
head(data_16)
dim(data_16)
data_16<-data_16[,-45]
#just 18s
data_18<-data_all[data_all$dataset=="18s",]
dim(data_18)
data_18<-data_18[,-45]

##
data_16
nrow(data_16)
#404 is gamma diversity = in all sponge samples
g_persample<-apply(data_16,2,sum)
g_persample<-as.matrix(g_persample)
meang_persample<-apply(g_persample,2,mean)

beta<-nrow(data_16)/meang_persample
beta
#2.435 higher differentiation between samples
data_18
nrow(data_18)
#404 is gamma diversity = in all sponge samples
g_persample18<-apply(data_18,2,sum)
g_persample18<-as.matrix(g_persample18)
meang_persample18<-apply(g_persample18,2,mean)

beta18<-nrow(data_18)/meang_persample18
beta18
#18s data much more differentiated

###Make graphs for Alpha diversity
#Turn tables back to counts before they were made into logical
data_all<-otus_and_taxab
#just 16s
data_16<-data_all[data_all$dataset=="16s",]
data_16<-data_16[,-c(1,2,47:55)]
#just 18s
data_18<-data_all[data_all$dataset=="18s",]
data_18<-data_18[,-c(1,2,47:55)]
#Code taken from RPubs Ecological Diversity Dr. Ro Allen M.Res. Marine Biology 2019
data_16<-as.matrix(t(data_16))
H<-diversity(data_16)
richness<-specnumber(data_16)
species<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
spcols<-c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green")
evenness<-H/log(richness)#Peilous richness
alpha<-cbind(shannon=H, richness=richness, pielou=evenness,species)
alpha<-as.data.frame(alpha)
str(alpha)
alpha$shannon<-as.numeric(alpha$shannon)
alpha$richness<-as.numeric(alpha$richness)
alpha$pielou<-as.numeric(alpha$pielou)

shanplot16<-ggplot(alpha,aes(species,shannon,colour=spcols))+
  geom_point()+
  ylab("Shannon's H'") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
richplot16<-ggplot(alpha,aes(species,richness,colour=spcols))+
  geom_point()+
  ylab("Species richness") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
pielouplot16<-ggplot(alpha,aes(species,pielou,colour=spcols))+
  geom_point()+
  ylab("Pielou's evenness") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
png("Figures/AS_16splots/alphadivplots16s.png",units="cm",width=30,height=10,res=300)
plot_grid(shanplot16,richplot16,pielouplot16,ncol=3)
dev.off()
#For 18s data
data_18<-as.matrix(t(data_18))

H18<-diversity(data_18)
richness18<-specnumber(data_18)
species<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
evenness18<-H18/log(richness18)#Peilous richness
alpha18<-cbind(shannon=H18, richness=richness18, pielou=evenness18,species)
alpha18<-as.data.frame(alpha18)
alpha18$shannon<-as.numeric(alpha18$shannon)
alpha18$richness<-as.numeric(alpha18$richness)
alpha18$pielou<-as.numeric(alpha18$pielou)

shanplot18<-ggplot(alpha18,aes(species,shannon,colour=spcols))+
  geom_point()+
  ylab("Shannon's H'") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
richplot18<-ggplot(alpha18,aes(species,richness,colour=spcols))+
  geom_point()+
  ylab("Species richness") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
pielouplot18<-ggplot(alpha18,aes(species,pielou,colour=spcols))+
  geom_point()+
  ylab("Pielou's evenness") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
png("Figures/18splots/alphadivplots18s.png",units="cm",width=30,height=10,res=300)
plot_grid(shanplot18,richplot18,pielouplot18,ncol=3)
dev.off()
####For 18 and 16 merged
data_all<-t(otus_and_taxab[,-c(1:2,47:55)])

Hb<-diversity(data_all)
richnessb<-specnumber(data_all)
species<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
evennessb<-Hb/log(richnessb)#Peilous richness
alphab<-cbind(shannon=Hb, richness=richnessb, pielou=evennessb,species)
alphab<-as.data.frame(alphab)
alphab$shannon<-as.numeric(alphab$shannon)
alphab$richness<-as.numeric(alphab$richness)
alphab$pielou<-as.numeric(alphab$pielou)

shanplotb<-ggplot(alphab,aes(species,shannon,colour=spcols))+
  geom_point()+
  ylab("Shannon's H'") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
richplotb<-ggplot(alphab,aes(species,richness,colour=spcols))+
  geom_point()+
  ylab("Species richness") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
pielouplotb<-ggplot(alphab,aes(species,pielou,colour=spcols))+
  geom_point()+
  ylab("Pielou's evenness") + 
  xlab("" )+
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Tau"="T.aurantium","Tme"="T.meloni","Tci"="T.citrina"))
png("Figures/16and18splots/alphadivplots16and18s.png",units="cm",width=30,height=10,res=300)
plot_grid(shanplotb,richplotb,pielouplotb,ncol=3)
dev.off()
