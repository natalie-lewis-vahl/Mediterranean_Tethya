library(vegan)
library(matrixStats)
library(dplyr)

otus<-read.csv("./Data/16and18s_otu.csv", sep=";")
head(otus)
otus<-arrange(otus,X.OTU.ID)
#BInd taxonomy to remove bacterial OTUs from 18s dataset 
taxa<-read.csv("./Data/16ans18s_taxa.csv",sep=";")
head(taxa)
taxa<-arrange(taxa,sequence_identifier)
otu_and_taxa<-bind_cols(otus,taxa)
#remove bacteria from 18s data set
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

#fILTER TO KEEP TOP 95%
#pct<-otus_and_taxa$countSum/sum(otus_and_taxa$countSum)
#otus_and_taxab<-cbind(otus_and_taxa, pct)
#otus_and_taxac<-otus_and_taxab%>%
#  arrange(desc(pct))
#otus_and_taxac<-otus_and_taxac%>%
#  mutate(cumpct=cumsum(pct))

#filtering: should be done seperate for both 16 and 18s data set
#Keep rows that have more than 50 counds and are from 16s data set OR have more than 5 counts and belong to 18s ddataset
otus_and_taxad<-otus_and_taxa[otus_and_taxa$countSum >= 50 & otus_and_taxa$dataset=="16s"|otus_and_taxa$countSum >= 5 & otus_and_taxa$dataset=="18s",]
head(otus_and_taxad)
dim(otus_and_taxad)
#otus_and_taxac<-otus_and_taxac[otus_and_taxac$cumpct < 0.95,]

oat<-otus_and_taxad[,-c(1,47:55)]
#rEMOVE SEQUENCE id AS A COLUMN AND MAKE THEM THE ROW NAMES
row.names(oat)<-oat$X.OTU.ID
#remove column with id names and taxa columns again
oat<-oat[,-1]

#groups 
#T aurantium: GW1941 till GW1951 inclusive +GW1984
#T meloni: GW1952 till GW1962 +GW1982 + GW1983
# T citroni: GW1963 till GW1981 inclusive

#group vectors for non bactloadCorrected datasets
num_groups<-c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,1)
label_groups<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")

head(oat)
dim(oat)
nmds<-metaMDS(t(oat), try=50)

jpeg("./Figures/16and18splots/ordination_nmds_plot.jpeg",height=8,width=8,units="in",res=300)

plot(nmds, type = "t", display="sites")

points(nmds, col=c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"), cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

dev.off()

x<-anosim(t(oat), num_groups)
summary(x)
x$statistic
hist(x$perm)
#vector 
species_vector<-data.frame(Species=c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau"))

#wFreqs_more50Counts % more1PctCounts
#temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","T","C","C","C","T")
#                           , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P"))

#adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)
meta_cca<-cca(t(oat)~Species, data=species_vector)
#bact_cca<-cca(t(bact_matrix[-c(17,18)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

jpeg("./Figures/16and18splots/ordination_cca_plot.jpeg",height=8,width=8,units="in",res=300)

plot(meta_cca, type="t", display="sites")
points(meta_cca, col=num_groups, cex=1.5, pch=16)
ordispider(meta_cca, label_groups, label=T)
dev.off()
anova(meta_cca)
anova(meta_cca, by="term")

anova(meta_cca, by="margin")

#highly variable otus

is.matrix(oat)
matrixoat<-as.matrix(oat)
matrix_rMeans<-rowMeans(matrixoat)
matrix_rSds<-rowSds(matrixoat)
matrix_zscores<-((matrixoat-matrix_rMeans)/matrix_rSds)

variability<-rowVars(matrixoat)
OTUVars<-cbind(oat,variability)
OTUVars<-arrange(OTUVars,by="variability")
variableOTUs<-rownames(head(OTUVars, n=25L))


jpeg("./Figures/16and18splots/heatmap16and18.jpeg",height=8,width=8,units="in",res=300)
heatmap(as.matrix(matrix_zscores),ColSideColors = c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green")) 
dev.off()
#here

jpeg("./Figures/16and18splots/mostvariable25_heatmap16and18.jpeg",height=8,width=8,units="in",res=300)

heatmap(matrix_zscores[rownames(matrix_zscores) %in% variableOTUs,],ColSideColors = c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"))

dev.off()

####
####Vegemite
data(dune)
sponges<-d18s_matrix[,-c(45,46,47)]

dis<-vegdist(t(sponges))
clus<-hclust(dis,"ward.D2")
plot(clus)
#single sites to be joined into large clusters
#or use complete for compact clusters or average -neutral grouping
den<-as.dendrogram(clus)
#ord<-cca(sponges)
wa<-scores(nmds,display="sites",choices=1)
oden<-reorder(den,wa,mean)
vegemite(t(sponges),use=oden,zero="-",scale="log")
tabasco(t(sponges),use=oden,zero="-",scale="Hill") #HIll is the scale publication standard

tabasco(t(sponges),use=oden, zero="-",col = heat.colors(5), scale="Hill")
