library(vegan)
library(matrixStats)
library(dplyr)
matrix_file<-"./Data/18s/all.otutab.csv"

otus<-read.csv("./Data/18s/all.otutab.csv", sep="\t")
#BInd taxonomy to remove bacteria
taxa<-read.csv("./Data/18s/taxa_fixed18s.csv",sep=",")
otus<-arrange(otus, X.OTU.ID)
taxa<-arrange(taxa,sequence_identifier)
otu_and_taxa<-bind_cols(otus,taxa)
otu_and_taxa<-otu_and_taxa%>%filter(Domain!="Bacteria")
#Take out OTU 1 for coral amplification
otus_and_taxa<-otu_and_taxa[!(otu_and_taxa$X.OTU.ID=="OTU_1"),]
#Uncomment to also take out OTU 2
otus_and_taxa<-otu_and_taxa[!(otu_and_taxa$X.OTU.ID=="OTU_2"),]
head(otus_and_taxa)
#Remove OTUS with less than 5 reads
countSum<-apply(otus_and_taxa[2:45],1,sum)
otus_and_taxa<-cbind(otus_and_taxa,countSum)
#same thing as using size column from taxa dataset
otus_and_taxa<-otus_and_taxa[otus_and_taxa$countSum>18,]

d18s_matrix<-otus_and_taxa[,-c(46:57)]
#rEMOVE SEQUENCE id AS A COLUMN AND MAKE THEM THE ROW NAMES
row.names(d18s_matrix)<-d18s_matrix$X.OTU.ID
#remove column with id names and taxa columns again
d18s_matrix<-d18s_matrix[,-1]

#if filtering is wanted at a more thorough level of 95% 
#countSum<-apply(d18s_matrix,1,sum)#sum the rows
#pct<-countSum/sum(countSum)

#d18s_matrix<-cbind(d18s_matrix, countSum, pct)
#d18s_matrix[desc(d18s_matrix$pct),]
#d18s_matrix<-d18s_matrix%>%
#  mutate(cumpct=cumsum(pct))

#d18s_matrix<-d18s_matrix%>%
#  arrange(desc(pct))
#d18s_matrix<-d18s_matrix[d18s_matrix$cumpct < 0.95,]

#358 rows with only OTU1 removed and 499 rows with OTU2 also removed

#groups 
#T aurantium: GW1941 till GW1951 inclusive 
#T meloni: GW1952 till GW1962
# T citroni: GW1963 till GW1984 inclusive

#group vectors for non bactloadCorrected datasets
num_groups<-c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
label_groups<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci")

head(d18s_matrix)
dim(d18s_matrix)
nmds<-metaMDS(t(d18s_matrix[,-c(45,46,47)]), try=50)


plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

x<-anosim(t(d18s_matrix[,-c(45,46,47)]), num_groups)
hist(x)
x
x$statistic
hist(x$perm)
abline(v=x$statistic)
#vector 
species_vector<-data.frame(Species=c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci"))

#wFreqs_more50Counts % more1PctCounts
#temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","T","C","C","C","T")
#                           , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P"))

#adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)
bact_cca<-cca(t(d18s_matrix[,-c(45,46,47)])~Species, data=species_vector)
#bact_cca<-cca(t(bact_matrix[-c(17,18)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

plot(bact_cca, type="t", display="sites")
points(bact_cca, col=num_groups, cex=1.5, pch=16)
ordispider(bact_cca, label_groups, label=T)
summary(bact_cca)
anova(bact_cca)
anova(bact_cca, by="term")

anova(bact_cca, by="margin")

#highly variable otus

is.matrix(d18s_matrix)
xd18s_matrix<-as.matrix(d18s_matrix)
head(xd18s_matrix)
matrix_rMeans<-rowMeans(xd18s_matrix[,-c(45,46,47)])
matrix_rSds<-rowSds(xd18s_matrix[,-c(45,46,47)])
matrix_zscores<-((xd18s_matrix[,-c(45,46,47)]-matrix_rMeans)/matrix_rSds)


OTUVars<-sort(rowVars(xd18s_matrix[,-c(45,46,47)],useNames = T),decreasing = T)
variableOTUs<-names(head(OTUVars, n=25L))
heatmap(as.matrix(matrix_zscores),ColSideColors = c("blue","red","red","blue","blue","blue","blue","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","blue","red","red","blue","blue","blue","blue","blue","blue","green","green","red","green","green","red","red","red","blue","blue","green","green","green","green","green","green","green"))

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
