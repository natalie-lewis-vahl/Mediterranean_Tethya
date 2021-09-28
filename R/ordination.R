library(vegan)

#setwd("/naslx/projects/uk213/lu43sur/cbas_temperature_bacts/R/")

matrix_file<-"./Data/cbas_tempVSctrl.otutab.csv"
#matrix_file<-"./cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv"

#matrix_file<-"./cbas_tempVSctrl_NoOTUCountLess50.otutab.csv"
#matrix_file<-"./cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly_wFreqs_more50Counts.csv"
#matrix_file<-"./cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly_wFreqs_more1PctCounts.csv"

bact_matrix<-read.csv(matrix_file, head=T, sep="\t", row.names=1)

countSum<-apply(bact_matrix,1,sum)
pct<-countSum/sum(countSum)

bact_matrix<-cbind(bact_matrix, countSum, pct)

dim(bact_matrix)


#if filtering is wanted:
#bact_matrix<-bact_matrix[bact_matrix$countSum > 50,]
#bact_matrix<-bact_matrix[bact_matrix$pct > 0.001,]


dim(bact_matrix)


#groups = GreenControl = 1
#         GreenTemp = 2
#         PurpleControl = 3
#         PurpleTemp = 4

#group vectors for non bactloadCorrected datasets
num_groups<-c(1,1,1,2,2,1,1,1,2,3,3,3,4,4,4,3,3,3,4,4)
label_groups<-c("GC","GC","GC","GT","GT","GC","GC","GC","GT","PC","PC","PC","PT","PT","PT","PC","PC","PC","PT","PT")

#wFreqs_more50Counts % more1PctCounts
#num_groups<-c(1,1,1,2,2,1,1,1,2,3,3,4,3,3,3,4)
#label_groups<-c("GC","GC","GC","GT","GT","GC","GC","GC","GT","PC","PC","PT","PC","PC","PC","PT")

#head(bact_matrix)

nmds<-metaMDS(t(bact_matrix[,-c(21,22)]), try=50)

#uncomment if using bact load corrected matrices
#nmds<-metaMDS(t(bact_matrix[,-c(17,18)]), try=50)

plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

#anosim(t(bact_matrix), num_groups)

#vector for non bactloadCorrected datasets
temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","C","T","T","T","C","C","C","T","T")
                            , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P","P","P","P","P"))

#wFreqs_more50Counts % more1PctCounts
#temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","T","C","C","C","T")
#                           , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P"))

#adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)

bact_cca<-cca(t(bact_matrix[,-c(21,22)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)
#bact_cca<-cca(t(bact_matrix[-c(17,18)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

plot(bact_cca, type="t", display="sites")
points(bact_cca, col=num_groups, cex=1.5, pch=16)
ordispider(bact_cca, label_groups, label=T)
anova(bact_cca)
anova(bact_cca, by="term")

anova(bact_cca, by="margin")

#cca with core only.
coreOTUs<-c("OTU_1","OTU_11","OTU_13","OTU_18","OTU_2","OTU_21","OTU_3","OTU_388","OTU_4","OTU_47","OTU_48","OTU_5","OTU_50","OTU_59","OTU_6","OTU_8","OTU_89")

coreOnlyCounts<-bact_matrix[rownames(bact_matrix[,-c(21,22)]) %in% as.character(coreOTUs),-c(21,22)]
core_cca<-cca(t(coreOnlyCounts)~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

plot(core_cca, type="t", display="sites")
points(core_cca, col=num_groups, cex=1.5, pch=16)
ordispider(core_cca, label_groups, label=T)

anova(core_cca)
anova(core_cca, by="term")
anova(core_cca, by="margin")


#cca with non-core community OTUs
nonCoreCounts<-bact_matrix[!(rownames(bact_matrix[,-c(21,22)]) %in% as.character(coreOTUs)),-c(21,22)]
nonCore_cca<-cca(t(nonCoreCounts)~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

plot(nonCore_cca, type="t", display="sites")
points(nonCore_cca, col=num_groups, cex=1.5, pch=16)
ordispider(nonCore_cca, label_groups, label=T)

anova(nonCore_cca)
anova(nonCore_cca, by="term")
anova(nonCore_cca, by="margin")


#vegemite(t(coreOnlyCounts), core_cca, scale = "log")


#highly variable otus


dim(bact_matrix)
is.matrix(xbact_matrix)
xbact_matrix<-as.matrix(bact_matrix)
bact_matrix_rMeans<-rowMeans(xbact_matrix[,-c(21,22)])
bact_matrix_rSDs<-rowSds(xbact_matrix[,-c(21,22)])
bact_matrix_zscores<-((xbact_matrix[,-c(21,22)]-bact_matrix_rMeans)/bact_matrix_rSDs)
#Search ## 

OTUVars<-sort(rowVars(xbact_matrix[,-c(21,22)],useNames=T),decreasing = T)
variableOTUs<-names(head(OTUVars, n=25L))
heatmap(as.matrix(bact_matrix_zscores))
heatmap(as.matrix(bact_matrix_zscores[rownames(bact_matrix_zscores) %in% variableOTUs,]))

bact_matrix_zscores[rownames(bact_matrix_zscores) %in% variableOTUs,]

