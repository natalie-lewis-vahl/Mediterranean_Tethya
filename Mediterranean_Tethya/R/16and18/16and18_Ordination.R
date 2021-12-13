library(vegan)
library(matrixStats)
library(dplyr)
matrix_file<-"./Data/16and18s_otu.csv"

otus<-read.csv("./Data/16and18s_otu.csv", sep=";")
#BInd taxonomy to remove bacterial OTUs from 18s dataset 
taxa<-read.csv("./Data/16ans18s_taxa.csv",sep=";")
head(otus)
otus<-arrange(otus, X.OTU.ID)
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

#fILTER TO KEEP TOP 95%
pct<-otus_and_taxa$countSum/sum(otus_and_taxa$countSum)
otus_and_taxab<-cbind(otus_and_taxa, pct)
otus_and_taxac<-otus_and_taxab%>%
  arrange(desc(pct))
otus_and_taxac<-otus_and_taxac%>%
  mutate(cumpct=cumsum(pct))
View(otus_and_taxac)
#filtering: should it be done seperate for both 16 and 18s data set or
#together with counts or with percentage- taken <97% 
#correspnding to exactly 19 reads or more
#bact_matrix<-bact_matrix[bact_matrix$countSum > 50,] 
head(bact_matrix)

otus_and_taxac<-otus_and_taxac[otus_and_taxac$cumpct < 0.97,]

oat<-otus_and_taxac[,-c(1,47:57)]
#rEMOVE SEQUENCE id AS A COLUMN AND MAKE THEM THE ROW NAMES
row.names(oat)<-oat$X.OTU.ID
#remove column with id names and taxa columns again
oat<-oat[,-1]

#groups 
#T aurantium: GW1941 till GW1951 inclusive 
#T meloni: GW1952 till GW1962
# T citroni: GW1963 till GW1984 inclusive

#group vectors for non bactloadCorrected datasets
num_groups<-c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
label_groups<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci")

head(oat)
dim(oat)
nmds<-metaMDS(t(oat), try=50)


plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

x<-anosim(t(oat), num_groups)
hist(x)
x
x$statistic
hist(x$perm)
#vector 
species_vector<-data.frame(Species=c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci"))

#wFreqs_more50Counts % more1PctCounts
#temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","T","C","C","C","T")
#                           , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P"))

#adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)
meta_cca<-cca(t(oat)~Species, data=species_vector)
#bact_cca<-cca(t(bact_matrix[-c(17,18)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

plot(meta_cca, type="t", display="sites")
points(meta_cca, col=num_groups, cex=1.5, pch=16)
ordispider(meta_cca, label_groups, label=T)
summary(meta_cca)
anova(meta_cca)
anova(meta_cca, by="term")

anova(meta_cca, by="margin")

#highly variable otus

is.matrix(oat)
matrixoat<-as.matrix(oat)
matrix_rMeans<-rowMeans(matrixoat)
matrix_rSds<-rowSds(matrixoat)
matrix_zscores<-((matrixoat-matrix_rMeans)/matrix_rSds)


OTUVars<-sort(rowVars(matrixoat,useNames = T),decreasing = T)
variableOTUs<-names(head(OTUVars, n=25L))
heatmap(as.matrix(matrix_zscores)) #ColSideColors = c("blue","red","red","blue","blue","blue","blue","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","blue","red","red","blue","blue","blue","blue","blue","blue","green","green","red","green","green","red","red","red","blue","blue","green","green","green","green","green","green","green"))
#here
??heatmap
heatmap(as.matrix(matrix_zscores[rownames(matrix_zscores) %in% variableOTUs,]))

matrix_zscores[rownames(matrix_zscores) %in% variableOTUs,]
View(otus_and_taxa)

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
