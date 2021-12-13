library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)

######################
#Private functions
#
#
# with minor modifications from Analytical Methods Using R and Excel
# by: Mark Gardener
# http://www.gardenersown.co.uk/Education/Lectures/Community%20Ecology%20Support%20Files.htm
#
#
####################################################
#significance test of differences in deviance of RAD models
rad_test<-function (x, conf.level = 0.95, log = TRUE) 
{
  if (inherits(x, what = "radfit.frame") == FALSE) 
    stop("\nInput must be the result of a radfit command\n")
  d <- sapply(x, function(y) unlist(lapply(y$models, deviance)))
  a <- sapply(x, function(y) unlist(lapply(y$models, AIC)))
  dev <- t(d)
  dev <- as.data.frame(dev)
  dev <- stack(dev)
  names(dev) <- c("deviance", "model")
  #levels(dev$model) <- c("Log", "Mb", "BS", "Pre", "Zip")
  if (log == TRUE) 
    dev.aov <- aov(log(deviance) ~ model, data = dev)
  else dev.aov <- aov(deviance ~ model, data = dev)
  dev.hsd <- TukeyHSD(dev.aov, ordered = TRUE, conf.level = conf.level)
  result <- list(original = dev, aov = dev.aov, post.hoc = dev.hsd, 
                 aic = a, deviance = d, log = log, conf.level = conf.level)
  class(result) <- "rad.htest"
  invisible(result)
}

print.rad.test<-function (x, ...) 
{
  devmean <- apply(x$deviance, MARGIN = 1, FUN = mean)
  devstdv <- apply(x$deviance, MARGIN = 1, FUN = sd)
  devSE <- devstdv/sqrt(ncol(x$deviance))
  aicmean <- apply(x$aic, MARGIN = 1, FUN = mean)
  result <- rbind(devmean, devSE, aicmean)
  rownames(result) <- c("Mean deviance", "SE of deviance", 
                        "Mean AIC")
  print(result, ...)
}

plot.rad.test<-function (x, which = "original", las = 1, ...) 
{
  WHICH <- c("original", "post.hoc")
  which <- tolower(which)
  which <- match.arg(which, WHICH)
  bpm <- x$original$model
  #bpm <- reorder(x$original$model, -x$original$deviance, FUN = mean)
  if (which == "original") {
    if (x$log == TRUE) 
      boxplot(log(deviance) ~ bpm, data = x$original, 
              las = las, ...)
    else boxplot(deviance ~ bpm, data = x$original, las = las, 
                 ...)
  }
  if (which == "post.hoc") {
    plot(x$post.hoc, las = las, ...)
  }
}

summary.rad.htest<-function (x, ...) 
{
  saov <- summary(x$aov)
  phoc <- x$post.hoc$model
  result <- list(ANOVA = saov, Tukey = phoc, Logarithmic = x$log, 
                 Confidence.Level = x$conf.level)
  print(result, ...)
}

################
#Input/Output paths
###########################

taxa<-read.csv("./Data/18s/taxa_fixed18s.csv",sep=",")
otus<-read.csv("./Data/18s/all.otutab.csv", sep="\t")
otus<-arrange(otus, X.OTU.ID)
taxa<-arrange(taxa,sequence_identifier)

plotsPath<-"./Figures/18splots/"


################
#User defined variables
###########################

plotSVG<-FALSE
plotPNG<-TRUE

################

#get OTU counts and add them to DF
countSum<-apply(otus[-1],1,sum)
countsDF<-cbind(otus, countSum)
head(countsDF)
#bind otu and taxonomy
cat<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], taxa[order(taxa$sequence_identifier),])
head(cat)
#FIlter out unnecessary data
##OR get rid of reads less than 3
cat<-cat[cat$countSum>5,]

#FIlter out bacteria
cat<-cat%>%filter(Domain!="Bacteria")

#remove otu 1 and 2 because of coamplification
#FInd which rows to delete
which(cat$X.OTU.ID=="OTU_1")
which(cat$X.OTU.ID=="OTU_2")
#Delete row 1 and 95
cat<-cat[-c(1,85),]

head(cat)
#Till which column do the samples go on to?
which(colnames(cat)=="GW1984")
# Phylum abundance per sample
cat[cat==""]<-"Unclassified"
#AND Getting rid of "unclassified" phylums
abundanceByPhylumBySample<-cat %>% group_by(Phylum) %>% filter(Phylum != "Unclassified") %>% mutate_at(colnames(cat)[2:45],sum) %>% distinct_at(vars(colnames(cat)[c(2:45,49)])) 
abundanceByClassBySample<-cat%>% group_by(Class)%>%filter(Class != "Unclassified") %>% mutate_at(colnames(cat)[2:45],sum) %>% distinct_at(vars(colnames(cat)[c(2:45,50)])) 
#########################
#Fit rank-abundance dominance models with phylum data
###########################################
library(lattice)
head(abundanceByPhylumBySample)
phylumRAD<-radfit(t(as.data.frame(abundanceByPhylumBySample)[,1:44]))
summary(phylumRAD)
phylumRADPlot<-plot(phylumRAD)
phylumRADPlot
s5rad<-radfit(t(as.data.frame(abundanceByPhylumBySample)[,5]))
View(abundanceByPhylumBySample)

##############3

#For specie T. auratis
s1tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,1]))
s2tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,2]))
s3tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,3]))
s4tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,4]))
s5tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,5]))#sample5 wont converge
s6tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,6]))
s7tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,7]))
s8tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,8]))
s9tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,9]))
s10tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,10]))
s11tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,11]))

phylumRADPlot
#For specie T.meloni
s12tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,12]))
s13tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,13]))
s14tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,14]))
s15tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,15]))
s16tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,16]))
s17tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,17]))
s18tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,18]))
s19tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,19]))
s20tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,20]))
s21tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,21]))
s22tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,22]))

###########
#For species T.citroni
s23tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,23]))
s24tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,24]))
s25tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,25]))
s26tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,26]))
s27tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,27]))
s28tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,28]))
s29tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,29]))
s30tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,30]))
s31tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,31]))
s32tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,32]))
s33tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,33]))
s34tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,34]))
s35tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,35]))
s36tciphylumRAD<-rad.null(t(as.data.frame(abundanceByPhylumBySample)[,36]))
s37tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,37]))
s38tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,38]))
s39tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,39]))
s40tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,40]))
s41tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,41]))
s42tciphylumRAD<-rad.lognormal(t(as.data.frame(abundanceByPhylumBySample)[,42]))
s43tciphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,43]))
s44tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,44]))

###############

#############
#Plot everything together
png("Figures/18splots/18s_phylumrad_grouped.png")
par(mfrow=c(3,4))
#For sp 1
max(s2tauphylumRAD$y)#Maximum for y is 6518

plot(s1tauphylumRAD,pch=1,lty=1,col=6,ylim=c(1,6520),xlim=c(0,12))
mtext(side=3,"Preemption",cex=1)
points(s3tauphylumRAD,col=6)
lines(s3tauphylumRAD,col=6)
points(s4tauphylumRAD,col=6)
lines(s4tauphylumRAD,col=6)
points(s6tauphylumRAD,col=6)
lines(s6tauphylumRAD,col=6)
points(s9tauphylumRAD,col=6)
lines(s9tauphylumRAD,col=6)
points(s10tauphylumRAD,col=6)
lines(s10tauphylumRAD,col=6)
points(s11tauphylumRAD,col=6)
lines(s11tauphylumRAD,col=6)

plot(s2tauphylumRAD,pch=1,lty=1,col=7,ylim=c(1,6520),xlim=c(0,12))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s5tauphylumRAD,col=7)
lines(s5tauphylumRAD,col=7)
points(s7tauphylumRAD,col=7)
lines(s7tauphylumRAD,col=7)
points(s8tauphylumRAD,col=7)
lines(s8tauphylumRAD,col=7)
plot.new()
plot.new()

#For sp 2 #phylumRADPlot
plot(s12tmephylumRAD,pch=1,lty=1,col=2,ylim=c(1,6250),xlim=c(0,12))
mtext(side=3,"Zipf",cex=1)
points(s13tmephylumRAD,col=2)
lines(s13tmephylumRAD,col=2)
points(s14tmephylumRAD,col=2)
lines(s14tmephylumRAD,col=2)
points(s15tmephylumRAD,col=2)
lines(s15tmephylumRAD,col=2)
points(s16tmephylumRAD,col=2)
lines(s16tmephylumRAD,col=2)
points(s17tmephylumRAD,col=2)
lines(s17tmephylumRAD,col=2)
points(s18tmephylumRAD,col=2)
lines(s18tmephylumRAD,col=2)
points(s19tmephylumRAD,col=2)
lines(s19tmephylumRAD,col=2)
points(s20tmephylumRAD,col=2)
lines(s20tmephylumRAD,col=2)
points(s21tmephylumRAD,col=2)
lines(s21tmephylumRAD,col=2)
points(s22tmephylumRAD,col=2)
lines(s22tmephylumRAD,col=2)
plot.new()
plot.new()
plot.new()

#########################
#For species 3
plot(s23tciphylumRAD,col=6,ylim=c(1,6520),xlim=c(0,12))
mtext(side=3,"Preemption",cex=1)
points(s25tciphylumRAD,col=6)
lines(s25tciphylumRAD,col=6)
points(s27tciphylumRAD,col=6)
lines(s27tciphylumRAD,col=6)
points(s28tciphylumRAD,col=6)
lines(s28tciphylumRAD,col=6)
points(s29tciphylumRAD,col=6)
lines(s29tciphylumRAD,col=6)
points(s30tciphylumRAD,col=6)
lines(s30tciphylumRAD,col=6)
points(s34tciphylumRAD,col=6)
lines(s34tciphylumRAD,col=6)
points(s39tciphylumRAD,col=6)
lines(s39tciphylumRAD,col=6)
points(s41tciphylumRAD,col=6)
lines(s41tciphylumRAD,col=6)
points(s44tciphylumRAD,col=6)
lines(s44tciphylumRAD,col=6)

plot(s25tciphylumRAD,col=2,ylim=c(1,6520),xlim=c(0,12))
mtext(side=3,"Zipf",cex=1)
points(s23tciphylumRAD,col=2)
lines(s23tciphylumRAD,col=2)
points(s24tciphylumRAD,col=2)
lines(s24tciphylumRAD,col=2)
points(s26tciphylumRAD,col=2)
lines(s26tciphylumRAD,col=2)
points(s31tciphylumRAD,col=2)
lines(s31tciphylumRAD,col=2)
points(s32tciphylumRAD,col=2)
lines(s32tciphylumRAD,col=2)
points(s33tciphylumRAD,col=2)
lines(s33tciphylumRAD,col=2)
points(s35tciphylumRAD,col=2)
lines(s35tciphylumRAD,col=2)
points(s37tciphylumRAD,col=2)
lines(s37tciphylumRAD,col=2)
points(s38tciphylumRAD,col=2)
lines(s38tciphylumRAD,col=2)
points(s40tciphylumRAD,col=2)
lines(s40tciphylumRAD,col=2)
points(s43tciphylumRAD,col=2)
lines(s43tciphylumRAD,col=2)

plot(s42tciphylumRAD,col=11,ylim=c(1,6520),xlim=c(0,12))
mtext(side=3,"Lognormal",cex=1)

plot(s36tciphylumRAD,col=4,ylim=c(1,6520),xlim=c(0,12))
mtext(side=3,"Null",cex=1)

dev.off()


phylumRAD_test<-rad_test(phylumRAD, conf.level = 0.95, log = T)
summary.rad.htest(phylumRAD_test)
print.rad.test(phylumRAD_test)
par(mfrow=c(1,2),mai=c(1,1,1,1))

plot.rad.test(phylumRAD_test)
plot.rad.test(phylumRAD_test, "post.hoc")
#For class- but misses all the blanks and unclassified excluded because if not it would falsely summarize it all into one group
classRAD<-radfit(t(as.data.frame(abundanceByClassBySample)[,1:44]))
summary(classRAD)
classRADPlot<-plot(classRAD)
classRADPlot
classRAD_test<-rad_test(classRAD, conf.level = 0.95, log = T )
summary.rad.htest(classRAD_test)
#Can't do the rad test for some reason??
##########################
#rank-abundance dominance with OTU data
##########################################################
#Without unclassified phyla
otuRAD<-radfit(t(as.data.frame(cat[(cat$Phylum != "Unclassified"),c(2:45)])))
summary(otuRAD)
otuRADPlot<-plot(otuRAD)

if(plotPNG){
  
  otuRADPlotPath<-paste(plotsPath,"otuRADPlot.png",sep="")
  png(otuRADPlotPath)
  print(otuRADPlot)
  dev.off()
  
}

if(plotSVG){
  
  otuRADPlotPath<-paste(plotsPath,"otuRADPlot.svg",sep="")
  svg(otuRADPlotPath)
  print(otuRADPlot)
  dev.off()
  
}

#To group the same model lines for each sp
##############
#For specie T. auratis
s1otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(2)]))
s2otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(3)]))
s3otuRADRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(4)]))
s4otuRADRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(5)]))
s5otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(6)]))
s6otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(7)]))
s7otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(8)]))
s8otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(9)]))
s9otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(10)]))
s10otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(11)]))
s11otuRAD<-rad.lognormal(t(cat[(cat$Phylum != ""),c(12)]))
#For specie T.meloni
s12otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(13)]))
s13otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(14)]))
s14otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(15)]))
s15otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(16)]))
s16otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(17)]))
s17otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(18)]))
s18otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(19)]))
s19otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(20)]))
s20otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(21)]))
s21otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(22)]))
s22otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(23)]))
#For species T.citroni
s23otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(24)]))
s24otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(25)]))
s25otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(26)]))
s26otuRAD<-rad.preempt(t(cat[(cat$Phylum != ""),c(27)]))
s27otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(28)]))
s28otuRAD<-rad.preempt(t(cat[(cat$Phylum != ""),c(29)]))
s29otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(30)]))
s30otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(31)]))
s31otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(32)]))
s32otuRAD<-rad.preempt(t(cat[(cat$Phylum != ""),c(33)]))
s33otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(34)]))
s34otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(35)]))
s35otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(36)]))
s36otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(37)]))
s37otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(38)]))
s38otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(39)]))
s39otuRAD<-rad.preempt(t(cat[(cat$Phylum != ""),c(40)]))
s40otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(41)]))
s41otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(42)]))
s42otuRAD<-rad.zipfbrot(t(cat[(cat$Phylum != ""),c(43)]))
s43otuRAD<-rad.lognormal(t(cat[(cat$Phylum != ""),c(44)]))
s44otuRAD<-rad.zipf(t(cat[(cat$Phylum != ""),c(45)]))
######
max(s8otuRAD$y)#y max 6397
View(s8otuRAD)
#Plot everything together
png("Figures/18splots/oturad_grouped.png")
par(mfrow=c(3,3))
#For sp 1
plot(s2otuRAD,pch=1,lty=1,col=2,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Zipf",cex=1)
points(s3otuRAD,col=2)
lines(s3otuRAD,col=2)
points(s5otuRAD,col=2)
lines(s5otuRAD,col=2)
points(s6otuRAD,col=2)
lines(s6otuRAD,col=2)
points(s7otuRAD,col=2)
lines(s7otuRAD,col=2)
points(s8otuRAD,col=2)
lines(s8otuRAD,col=2)

plot(s1otuRAD,pch=1,lty=1,col=7,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s4otuRAD,col=7)
lines(s4otuRAD,col=7)
points(s9otuRAD,col=7)
lines(s9otuRAD,col=7)
points(s10otuRAD,col=7)
lines(s10otuRAD,col=7)

plot(s11tciphylumRAD,col=11,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Lognormal",cex=1)

#For sp 2 #phylumRADPlot
plot(s12otuRAD,pch=1,lty=1,col=2,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Zipf",cex=1)
points(s14otuRAD,col=2)
lines(s14otuRAD,col=2)
points(s15otuRAD,col=2)
lines(s15otuRAD,col=2)
points(s16otuRAD,col=2)
lines(s16otuRAD,col=2)
points(s15otuRAD,col=2)
lines(s15otuRAD,col=2)
points(s18otuRAD,col=2)
lines(s18otuRAD,col=2)
points(s22otuRAD,col=2)
lines(s22otuRAD,col=2)


plot(s13otuRAD,pch=1,lty=1,col=7,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s22otuRAD,col=7)
lines(s22otuRAD,col=7)
points(s20otuRAD,col=7)
lines(s20otuRAD,col=7)
points(s21otuRAD,col=7)
lines(s21otuRAD,col=7)

#########################
#For species 3
plot(s26tciphylumRAD,col=6,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Preemption",cex=1)
points(s28otuRAD,col=6)
lines(s28otuRAD,col=6)
points(s32otuRAD,col=6)
lines(s32otuRAD,col=6)

plot(s23otuRAD,col=2,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Zipf",cex=1)
points(s24otuRAD,col=2)
lines(s24otuRAD,col=2)
points(s25otuRAD,col=2)
lines(s25otuRAD,col=2)
points(s27otuRAD,col=2)
lines(s27otuRAD,col=2)
points(s29otuRAD,col=2)
lines(s29otuRAD,col=2)
points(s30otuRAD,col=2)
lines(s30otuRAD,col=2)
points(s31otuRAD,col=2)
lines(s31otuRAD,col=2)
points(s33otuRAD,col=2)
lines(s33otuRAD,col=2)
points(s34otuRAD,col=2)
lines(s34otuRAD,col=2)
points(s35otuRAD,col=2)
lines(s35otuRAD,col=2)
points(s36otuRAD,col=2)
lines(s36otuRAD,col=2)
points(s38otuRAD,col=2)
lines(s38otuRAD,col=2)
points(s40otuRAD,col=2)
lines(s40otuRAD,col=2)
points(s41otuRAD,col=2)
lines(s41otuRAD,col=2)
points(s44otuRAD,col=2)
lines(s44otuRAD,col=2)

plot(s37otuRAD,pch=1,lty=1,col=7,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s42otuRAD,col=7)
lines(s42otuRAD,col=7)


plot(s43otuRAD,col=11,ylim=c(1,6400),xlim=c(0,110))
mtext(side=3,"Lognormal",cex=1)



dev.off()


otuRAD_test<-rad_test(otuRAD, conf.level = 0.95, log = T)
summary.rad.htest(otuRAD_test)
print.rad.test(otuRAD_test)
par(mfrow=c(1,2),mai=c(1,1,1,1))
plot.rad.test(otuRAD_test)
plot.rad.test(otuRAD_test, "post.hoc")
#######

# AND
#########################
#Fit rank-abundance dominance models  for NOT Getting rid of OTUs belonging to "unclassified" phyla
###########################################
unique(cat$Phylum)
##########################
#rank-abundance dominance with OTU data
#With unclassified phyla
head(cat)
xotuRAD<-radfit(t(cat[c(2:45)]))

summary(xotuRAD)
xotuRADPlot<-plot(xotuRAD)
xotuRADPlot

if(plotPNG){
  
  xotuRADPlotPath<-paste(plotsPath,"xotuRADPlot.png",sep="")
  png(xotuRADPlotPath)
  print(xotuRADPlot)
  dev.off()
  
}

if(plotSVG){
  
  xotuRADPlotPath<-paste(plotsPath,"xotuRADPlot.svg",sep="")
  svg(xotuRADPlotPath)
  print(xotuRADPlot)
  dev.off()
  
}

#columns 5 11 and 22 didnt't converge so remove sample 4, 10 and 21

#For specie T. auratis
xs1otuRAD<-rad.zipfbrot(t(cat[c(2)]))
xs2otuRAD<-rad.zipf(t(cat[c(3)]))
xs3otuRAD<-rad.zipf(t(cat[c(4)]))
# xs4otuRAD<-rad.zipfbrot(t(cat[c(5)])) didn't converge
xs5otuRAD<-rad.zipf(t(cat[c(6)]))
xs6otuRAD<-rad.zipf(t(cat[c(7)]))
xs7otuRAD<-rad.zipf(t(cat[c(8)]))
xs8otuRAD<-rad.zipf(t(cat[c(9)]))
xs9otuRAD<-rad.zipfbrot(t(cat[c(10)]))
#xs10otuRAD<-rad.zipfbrot(t(cat[c(11)]))#didn't converge
xs11otuRAD<-rad.lognormal(t(cat[c(12)]))
#For specie T.meloni
xs12otuRAD<-rad.zipf(t(cat[c(13)]))
xs13otuRAD<-rad.zipfbrot(t(cat[c(14)]))
xs14otuRAD<-rad.zipf(t(cat[c(15)]))
xs15otuRAD<-rad.zipf(t(cat[c(16)]))
xs16otuRAD<-rad.zipf(t(cat[c(17)]))
xs17otuRAD<-rad.zipf(t(cat[c(18)]))
xs18otuRAD<-rad.zipf(t(cat[c(19)]))
xs19otuRAD<-rad.zipfbrot(t(cat[c(20)]))
xs20otuRAD<-rad.zipfbrot(t(cat[c(21)]))
#xs21otuRAD<-rad.zipfbrot(t(cat[c(22)])) didn't converge
xs22otuRAD<-rad.zipf(t(cat[c(23)]))
#For species T.citroni
xs23otuRAD<-rad.zipf(t(cat[c(24)]))
xs24otuRAD<-rad.zipf(t(cat[c(25)]))
xs25otuRAD<-rad.zipf(t(cat[c(26)]))
xs26otuRAD<-rad.preempt(t(cat[c(27)]))
xs27otuRAD<-rad.zipf(t(cat[c(28)]))
xs28otuRAD<-rad.preempt(t(cat[c(29)]))
xs29otuRAD<-rad.zipf(t(cat[c(30)]))
xs30otuRAD<-rad.zipf(t(cat[c(31)]))
xs31otuRAD<-rad.zipf(t(cat[c(32)]))
xs32otuRAD<-rad.preempt(t(cat[c(33)]))
xs33otuRAD<-rad.zipf(t(cat[c(34)]))
xs34otuRAD<-rad.zipf(t(cat[c(35)]))
xs35otuRAD<-rad.zipf(t(cat[c(36)]))
xs36otuRAD<-rad.zipf(t(cat[c(37)]))
xs37otuRAD<-rad.zipfbrot(t(cat[c(38)]))
xs38otuRAD<-rad.zipf(t(cat[c(39)]))
xs39otuRAD<-rad.preempt(t(cat[c(40)]))
xs40otuRAD<-rad.zipf(t(cat[c(41)]))
xs41otuRAD<-rad.zipf(t(cat[c(42)]))
xs42otuRAD<-rad.zipfbrot(t(cat[c(43)]))
xs43otuRAD<-rad.lognormal(t(cat[c(44)]))
xs44otuRAD<-rad.zipf(t(cat[c(45)]))
######
#Plot everything together

#Plot everything together
png("./Figures/18splots/xoturad_grouped.png")
par(mfrow=c(3,4))#13650
#For sp 1
plot(xs2otuRAD,pch=1,lty=1,col=2,ylim=c(1,13650),xlim=c(0,110))
mtext(side=3,"Zipf",cex=1)
points(xs3otuRAD,col=2)
lines(xs3otuRAD,col=2)
points(xs5otuRAD,col=2)
lines(xs5otuRAD,col=2)
points(xs6otuRAD,col=2)
lines(xs6otuRAD,col=2)
points(xs7otuRAD,col=2)
lines(xs7otuRAD,col=2)
points(xs8otuRAD,col=2)
lines(xs8otuRAD,col=2)

plot(xs1otuRAD,pch=1,lty=1,col=7,ylim=c(1,13650),xlim=c(0,110))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(xs9otuRAD,col=7)
lines(xs9otuRAD,col=7)

plot(xs11otuRAD,pch=1,lty=1,col=11,ylim=c(1,13650),xlim=c(0,110))
mtext(side=3,"Lognormal",cex=1)

plot.new()

#For sp 2 #phylumRADPlot
plot(xs12otuRAD,pch=1,lty=1,col=2,ylim=c(1,13650),xlim=c(0,110))
mtext(side=3,"Zipf",cex=1)
points(xs14otuRAD,col=2)
lines(xs14otuRAD,col=2)
points(xs15otuRAD,col=2)
lines(xs15otuRAD,col=2)
points(xs16otuRAD,col=2)
lines(xs16otuRAD,col=2)
points(xs17otuRAD,col=2)
lines(xs17otuRAD,col=2)
points(xs18otuRAD,col=2)
lines(xs18otuRAD,col=2)
points(xs22otuRAD,col=2)
lines(xs22otuRAD,col=2)
plot(xs13otuRAD,pch=1,lty=1,col=7,ylim=c(1,13650),xlim=c(0,110))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(xs19otuRAD,col=7)
lines(xs19otuRAD,col=7)
points(xs20otuRAD,col=7)
lines(xs20otuRAD,col=7)
plot.new()
plot.new()
#########################
#For species 3
plot(xs26otuRAD,col=6,ylim=c(1,13650),xlim=c(0,110))
mtext(side=3,"Preemption",cex=1)
points(xs28otuRAD,col=6)
lines(xs28otuRAD,col=6)
points(xs32otuRAD,col=6)
lines(xs32otuRAD,col=6)
points(xs39otuRAD,col=6)
lines(xs39otuRAD,col=6)

plot(xs23otuRAD,pch=1,lty=1,col=2,ylim=c(1,13650),xlim=c(0,110))
mtext(side=3,"Zipf",cex=1)
points(xs24otuRAD,col=2)
lines(xs24otuRAD,col=2)
points(xs25otuRAD,col=2)
lines(xs25otuRAD,col=2)
points(xs27otuRAD,col=2)
lines(xs27otuRAD,col=2)
points(xs29otuRAD,col=2)
lines(xs29otuRAD,col=2)
points(xs30otuRAD,col=2)
lines(xs30otuRAD,col=2)
points(xs31otuRAD,col=2)
lines(xs31otuRAD,col=2)
points(xs33otuRAD,col=2)
lines(xs33otuRAD,col=2)
points(xs34otuRAD,col=2)
lines(xs34otuRAD,col=2)
points(xs35otuRAD,col=2)
lines(xs35otuRAD,col=2)
points(xs36otuRAD,col=2)
lines(xs36otuRAD,col=2)
points(xs38otuRAD,col=2)
lines(xs38otuRAD,col=2)
points(xs40otuRAD,col=2)
lines(xs40otuRAD,col=2)
points(xs44otuRAD,col=2)
lines(xs44otuRAD,col=2)

plot(xs43otuRAD,col=11,ylim=c(1,6520),xlim=c(0,110))
mtext(side=3,"Lognormal",cex=1)

plot(xs37otuRAD,col=7,ylim=c(1,6520),xlim=c(0,110))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(xs42otuRAD,col=7)
lines(xs42otuRAD,col=7)


dev.off()
#5 11 and 22 didnt't converge so remove these
#head(cat)
xotuRAD<-radfit(t(cat[c(2,3,4,6:10,12:21,23:45)]))
summary(xotuRAD)

xotuRAD_test<-rad_test(xotuRAD, conf.level = 0.95, log = T )
??rad_test()
summary.rad.htest(xotuRAD_test)
par(mfrow=c(1,2),mai=c(1,1,1,1))
print.rad.test(xotuRAD_test)
plot.rad.test(xotuRAD_test)
plot.rad.test(xotuRAD_test, "post.hoc")



