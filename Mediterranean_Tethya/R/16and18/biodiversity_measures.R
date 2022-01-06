#Alpha, beta and gamma
library(matrixStats)
library(dplyr)
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
coamplif_sponge$sequence_identifier[2]#OTU_4206
otus_and_taxa<-otus_and_taxa[!(otus_and_taxa$X.OTU.ID=="OTU_4206"),]

#filtering: should be done seperate for both 16 and 18s data set
#Keep rows that have more than 50 counds and are from 16s data set OR have more than 5 counts and belong to 18s ddataset
otus_and_taxab<-otus_and_taxa[otus_and_taxa$countSum > 50 & otus_and_taxa$dataset=="16s"|otus_and_taxa$countSum > 5 & otus_and_taxa$dataset=="18s",]
head(otus_and_taxab)
data_all<-otus_and_taxab
View(data_all)
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
sapply(data_16)