library(vegan)
library(matrixStats)
library(dplyr)
matrix_file<-"./Data/all.otutab_raw.csv"

bact_matrix<-read.csv(matrix_file, head=T, sep="\t", row.names=1)
countSum<-apply(bact_matrix,1,sum)#sum the rows
pct<-countSum/sum(countSum)

bact_matrix<-cbind(bact_matrix, countSum, pct)
bact_matrix<-bact_matrix%>%
  arrange(desc(pct))
bact_matrix<-bact_matrix%>%
  mutate(cumpct=cumsum(pct))
dim(bact_matrix)

View(bact_matrix)
#if filtering is wanted:(doesn't have large enough sample size)
#bact_matrix<-bact_matrix[bact_matrix$countSum > 50,] 
head(bact_matrix)

bact_matrix<-bact_matrix[bact_matrix$cumpct < 0.95,]


dim(bact_matrix)


#groups 
#Tau:"GW1941" "GW1942" "GW1944" "GW1945" "GW1946" "GW1947" "GW1948"
#Tme:"GW1952" "GW1953" "GW1954" "GW1955" "GW1957" "GW1958" "GW1959"
#Tci:"GW1964" "GW1967" "GW1968" "GW1969" "GW1973" "GW1979" "GW1982"

#group vectors for non bactloadCorrected datasets
num_groups<-c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
label_groups<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci")

#head(bact_matrix)

nmds<-metaMDS(t(bact_matrix[,-c(22,23,24)]), try=50)

#With filtered data get Warning message: stress is (nearly) zero: you may have insufficient data

plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

#anosim(t(bact_matrix), num_groups)

#vector 
species_vector<-data.frame(Species=c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci"))

#wFreqs_more50Counts % more1PctCounts
#temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","T","C","C","C","T")
#                           , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P"))

#adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)
??cca
bact_cca<-cca(t(bact_matrix[,-c(22,23,24)])~Species, data=species_vector)
#bact_cca<-cca(t(bact_matrix[-c(17,18)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

plot(bact_cca, type="t", display="sites")
points(bact_cca, col=num_groups, cex=1.5, pch=16)
ordispider(bact_cca, label_groups, label=T)
summary(bact_cca)
anova(bact_cca)
anova(bact_cca, by="term")

anova(bact_cca, by="margin")

#cca with all core species inclusive (core of any of the species).
coreOTUs<-c("OTU_1", "OTU_10","OTU_102","OTU_105","OTU_1068", "OTU_1074", "OTU_1095", "OTU_11","OTU_110","OTU_12", "OTU_128", "OTU_129","OTU_13",
            "OTU_1379","OTU_14","OTU_151","OTU_16","OTU_17","OTU_170","OTU_179","OTU_19","OTU_2","OTU_20","OTU_2075","OTU_2078","OTU_21","OTU_210",
            "OTU_22", "OTU_235","OTU_24", "OTU_241","OTU_25","OTU_26","OTU_261","OTU_27","OTU_295","OTU_30","OTU_31", "OTU_32", "OTU_325","OTU_326",
            "OTU_33", "OTU_34", "OTU_35", "OTU_36", "OTU_374", "OTU_39", "OTU_41", "OTU_42", "OTU_43", "OTU_455",  "OTU_500", 
            "OTU_504","OTU_525","OTU_53", "OTU_532","OTU_54", "OTU_56","OTU_57","OTU_60","OTU_64","OTU_66", "OTU_7","OTU_716","OTU_718",
            "OTU_73", "OTU_8","OTU_809","OTU_81", "OTU_85", "OTU_86", "OTU_87", "OTU_9","OTU_907","OTU_92", "OTU_931")
coreOnlyCounts<-bact_matrix[rownames(bact_matrix[,-c(22,23,24)]) %in% as.character(coreOTUs),-c(22,23,24)]
core_cca<-cca(t(coreOnlyCounts)~Species, data=species_vector)#correspondence analysis,
plot(core_cca, type="t", display="sites")
points(core_cca, col=num_groups, cex=1.5, pch=16)
ordispider(core_cca, label_groups, label=T)

anova(core_cca)



#cca with non-core community OTUs (those not included in the all core species in any species)
nonCoreCounts<-bact_matrix[!(rownames(bact_matrix[,-c(22,23,24)]) %in% as.character(coreOTUs)),-c(22,23,24)]
nonCore_cca<-cca(t(nonCoreCounts)~Species, data=species_vector)

plot(nonCore_cca, type="t", display="sites")
points(nonCore_cca, col=num_groups, cex=1.5, pch=16)
ordispider(nonCore_cca, label_groups, label=T)

anova(nonCore_cca)



#vegemite(t(coreOnlyCounts), core_cca, scale = "log")


#highly variable otus
head(bact_matrix)
is.matrix(bact_matrix)
xbact_matrix<-as.matrix(bact_matrix)
bact_matrix_rMeans<-rowMeans(xbact_matrix[,-c(22,23,24)])
bact_matrix_rSds<-rowSds(xbact_matrix[,-c(22,23,24)])
bact_matrix_zscores<-((xbact_matrix[,-c(22,23,24)]-bact_matrix_rMeans)/bact_matrix_rSds)


OTUVars<-sort(rowVars(xbact_matrix[,-c(22,23,24)],useNames = T),decreasing = T)
variableOTUs<-names(head(OTUVars, n=25L))
heatmap(as.matrix(bact_matrix_zscores),ColSideColors = c("green","green","green","green","green","green","green","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue"))
heatmap(as.matrix(bact_matrix_zscores[rownames(bact_matrix_zscores) %in% variableOTUs,]),ColSideColors = c("green","green","green","green","green","green","green","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue"))
bact_matrix_zscores[rownames(bact_matrix_zscores) %in% variableOTUs,]
