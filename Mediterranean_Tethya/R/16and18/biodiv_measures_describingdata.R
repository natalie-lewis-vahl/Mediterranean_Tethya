#just 16s
data_16<-data_all[data_all$dataset=="16s",]
head(data_16)
dim(data_16)
data_16<-data_16[,-45]


otus16<-otus_and_taxab[otus_and_taxab$dataset=="16s",2]
head(otus16)
data_all<-as.data.frame(otus_and_taxab[otus_and_taxab$dataset=="16s",-c(1:2,47:55)]>0)
yorn16<-cbind(otus16,data_all)
head(yorn18)
which(yorn16$otus16=="OTU_1"|yorn16$otus16=="OTU_2") #row 1 and 102

yorn16[c(1,102),]

#####

otus18<-otus_and_taxab[otus_and_taxab$dataset=="18s",2]
head(otus18)
data_all<-as.data.frame(otus_and_taxab[otus_and_taxab$dataset=="18s",-c(1:2,47:55)]>0)
yorn18<-cbind(otus18,data_all)
head(yorn18)
which(yorn16$otus16=="OTU_1"|yorn16$otus16=="OTU_2") #row 1 and 102

data_all<-as.data.frame(otus_and_taxa[otus_and_taxa$dataset=="18s",-c(1:2,47:55)]>0)
head(otus_and_taxab)
