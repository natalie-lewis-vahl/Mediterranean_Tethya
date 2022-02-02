library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)

######################
#Private functions
#
#Modified from abundance models script from Sergio Vargas et al
#https://gitlab.lrz.de/cbas/CBAS_16S/-/tree/master/CBAS_HeatStress
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

taxa<-read.csv("./Data/16s_allsamples_taxa.csv",sep=";")
otus<-read.csv("./Data/16s_allsamples_otu.csv", sep=";")


plotsPath<-"./Figures/AS_16splots/"


################

#get OTU counts and add them to DF
countSum<-apply(otus[-1],1,sum)
countsDF<-cbind(otus, countSum)
head(countsDF)
head(taxa)
#bind otu and taxonomy
cat<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], taxa[order(taxa$sequence_identifier),])
head(cat)
#FIlter out unnecessary data
##OR get rid of reads less than 50
cat<-cat[cat$countSum>=50,]

head(cat)
#Till which column do the samples go on to?
which(colnames(cat)=="GW1984")
#column 45
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

##############
#For specie T. auratis
s1tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,1]))
s2tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,2]))
s3tauphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,3]))
s4tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,4]))
s5tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,5]))#sample5 wont converge
s6tauphylumRAD<-rad.lognormal(t(as.data.frame(abundanceByPhylumBySample)[,6]))
s7tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,7]))
s8tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,8]))
s9tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,9]))
s10tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,10]))
s11tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,11]))
s44tauphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,44]))

phylumRADPlot
#For specie T.meloni
s12tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,12]))
s13tmephylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,13]))
s14tmephylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,14]))
s15tmephylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,15]))
s16tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,16]))
s17tmephylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,17]))
s18tmephylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,18]))
s19tmephylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,19]))
s20tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,20]))
s21tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,21]))
s22tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,22]))
s42tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,42]))
s43tmephylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,43]))

###########
#For species T.citroni
s23tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,23]))
s24tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,24]))
s25tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,25]))
s26tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,26]))
s27tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,27]))
s28tciphylumRAD<-rad.lognormal(t(as.data.frame(abundanceByPhylumBySample)[,28]))
s29tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,29]))
s30tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,30]))
s31tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,31]))
s32tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,32]))
s33tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,33]))
s34tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,34]))
s35tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,35]))
s36tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,36]))
s37tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,37]))
s38tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,38]))
s39tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,39]))
s40tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,40]))
s41tciphylumRAD<-rad.lognormal(t(as.data.frame(abundanceByPhylumBySample)[,41]))
###############

#############
#Plot everything together
png("Figures/AS_16splots/AS16s_phylumrad_grouped.png",width=8,height=8,units="in",res=300)
par(mfrow=c(3,4))
#For sp 1
max(s4tauphylumRAD$y)#Maximum for y is 24556

plot(s3tauphylumRAD,pch=1,lty=1,col=6,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3, "Preemption",cex=1)

plot(s1tauphylumRAD,pch=1,lty=1,col=2,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Zipf",cex=1)
points(s2tauphylumRAD,col=2)
lines(s2tauphylumRAD,col=2)
points(s4tauphylumRAD,col=2)
lines(s4tauphylumRAD,col=2)
points(s6tauphylumRAD,col=2)
lines(s6tauphylumRAD,col=2)
points(s7tauphylumRAD,col=2)
lines(s7tauphylumRAD,col=2)
points(s8tauphylumRAD,col=2)
lines(s8tauphylumRAD,col=2)


plot(s9tauphylumRAD,pch=1,lty=1,col=7,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s10tauphylumRAD,col=7)
lines(s10tauphylumRAD,col=7)
points(s11tauphylumRAD,col=7)
lines(s11tauphylumRAD,col=7)
points(s44tauphylumRAD,col=7)
lines(s44tauphylumRAD,col=7)

plot(s5tauphylumRAD,pch=1,lty=1,col=11,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Lognormal",cex=1)

#For sp 2
plot(s13tmephylumRAD,pch=1,lty=1,col=6,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Preemption",cex=1)
points(s19tmephylumRAD,col=6)
lines(s19tmephylumRAD,col=6)
points(s14tmephylumRAD,col=6)
lines(s14tmephylumRAD,col=6)
points(s15tmephylumRAD,col=6)
lines(s15tmephylumRAD,col=6)
points(s17tmephylumRAD,col=6)
lines(s17tmephylumRAD,col=6)
points(s18tmephylumRAD,col=6)
lines(s18tmephylumRAD,col=6)

plot(s16tmephylumRAD,pch=1,lty=1,col=2,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Zipf",cex=1)
points(s20tmephylumRAD,col=2)
lines(s20tmephylumRAD,col=2)
points(s21tmephylumRAD,col=2)
lines(s21tmephylumRAD,col=2)
points(s22tmephylumRAD,col=2)
lines(s22tmephylumRAD,col=2)

plot(s12tmephylumRAD,pch=1,lty=1,col=7,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s43tmephylumRAD,col=7)
lines(s43tmephylumRAD,col=7)

plot(s42tciphylumRAD,col=11,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Lognormal",cex=1)
#########################
#For species 3

plot(s27tciphylumRAD,col=6,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Preemption",cex=1)
points(s38tciphylumRAD,col=6)
lines(s38tciphylumRAD,col=6)

plot(s23tciphylumRAD,col=7,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s24tciphylumRAD,col=7)
lines(s24tciphylumRAD,col=7)
points(s25tciphylumRAD,col=7)
lines(s25tciphylumRAD,col=7)
points(s26tciphylumRAD,col=7)
lines(s26tciphylumRAD,col=7)
points(s29tciphylumRAD,col=7)
lines(s29tciphylumRAD,col=7)
points(s30tciphylumRAD,col=7)
lines(s30tciphylumRAD,col=7)
points(s31tciphylumRAD,col=7)
lines(s31tciphylumRAD,col=7)
points(s32tciphylumRAD,col=7)
lines(s32tciphylumRAD,col=7)
points(s33tciphylumRAD,col=7)
lines(s33tciphylumRAD,col=7)
points(s34tciphylumRAD,col=7)
lines(s34tciphylumRAD,col=7)
points(s35tciphylumRAD,col=7)
lines(s35tciphylumRAD,col=7)
points(s36tciphylumRAD,col=7)
lines(s36tciphylumRAD,col=7)
points(s37tciphylumRAD,col=7)
lines(s37tciphylumRAD,col=7)
points(s39tciphylumRAD,col=7)
lines(s39tciphylumRAD,col=7)
points(s40tciphylumRAD,col=7)
lines(s40tciphylumRAD,col=7)


plot(s28tciphylumRAD,col=11,ylim=c(1,25000),xlim=c(0,30))
mtext(side=3,"Lognormal",cex=1)
points(s41tciphylumRAD,col=11)
lines(s41tciphylumRAD,col=11)

plot.new()
dev.off()
########

phylumRAD_test<-rad_test(phylumRAD, conf.level = 0.95, log = T)
summary.rad.htest(phylumRAD_test)
print.rad.test(phylumRAD_test)
png("Figures/AS_16splots/phylumradtest.png",height=20,width=30,units="cm",res=300)
par(mfrow=c(1,2),mai=c(1,1,1,1))
plot.rad.test(phylumRAD_test)
plot.rad.test(phylumRAD_test, "post.hoc")
dev.off()
#For class- but misses all the blanks and unclassified excluded because if not it would falsely summarize it all into one group
classRAD<-radfit(t(as.data.frame(abundanceByClassBySample)[,1:44]))
summary(classRAD)
classRADPlot<-plot(classRAD)
classRADPlot
##########################
#rank-abundance dominance with OTU data
##########################################################
#Without unclassified phyla
otuRAD<-radfit(t(as.data.frame(cat[(cat$Phylum != "Unclassified"),c(2:45)])))
summary(otuRAD)
otuRADPlot<-plot(otuRAD)
png("Figures/AS_16splots/otuRADplot_no_unclassified.png",height=20,width=10,units="cm",res=300)
classRADPlot
dev.off()

  
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
png("./Figures/AS_16splots/radplototu.png",width=8,height=8,units="in",res=300)
xotuRADPlot
dev.off()
#For specie T. auratis
xs1otuRAD<-rad.zipf(t(cat[c(2)]))
xs2otuRAD<-rad.zipf(t(cat[c(3)]))
xs3otuRAD<-rad.zipfbrot(t(cat[c(4)]))
xs4otuRAD<-rad.zipf(t(cat[c(5)]))
xs5otuRAD<-rad.zipf(t(cat[c(6)]))
xs6otuRAD<-rad.zipf(t(cat[c(7)]))
xs7otuRAD<-rad.zipf(t(cat[c(8)]))
xs8otuRAD<-rad.zipf(t(cat[c(9)]))
xs9otuRAD<-rad.zipf(t(cat[c(10)]))
xs10otuRAD<-rad.zipf(t(cat[c(11)]))
xs11otuRAD<-rad.lognormal(t(cat[c(12)]))
xs44otuRAD<-rad.preempt(t(cat[c(45)]))

#For specie T.meloni
xs12otuRAD<-rad.zipfbrot(t(cat[c(13)]))
xs13otuRAD<-rad.zipfbrot(t(cat[c(14)]))
xs14otuRAD<-rad.zipfbrot(t(cat[c(15)]))#COULD NOT CONVERGE
xs15otuRAD<-rad.lognormal(t(cat[c(16)]))
xs16otuRAD<-rad.zipf(t(cat[c(17)]))
xs17otuRAD<-rad.preempt(t(cat[c(18)]))
xs18otuRAD<-rad.lognormal(t(cat[c(19)]))
xs19otuRAD<-rad.zipfbrot(t(cat[c(20)]))
xs20otuRAD<-rad.zipf(t(cat[c(21)]))
xs21otuRAD<-rad.zipf(t(cat[c(22)]))
xs22otuRAD<-rad.zipf(t(cat[c(23)]))
xs42otuRAD<-rad.lognormal(t(cat[c(43)]))
xs43otuRAD<-rad.zipfbrot(t(cat[c(44)]))
#For species T.citrina
xs23otuRAD<-rad.zipf(t(cat[c(24)]))
xs24otuRAD<-rad.zipfbrot(t(cat[c(25)]))
xs25otuRAD<-rad.zipf(t(cat[c(26)]))
xs26otuRAD<-rad.zipf(t(cat[c(27)]))
xs27otuRAD<-rad.zipfbrot(t(cat[c(28)]))
xs28otuRAD<-rad.zipfbrot(t(cat[c(29)]))
xs29otuRAD<-rad.zipf(t(cat[c(30)]))
xs30otuRAD<-rad.zipf(t(cat[c(31)]))
xs31otuRAD<-rad.zipf(t(cat[c(32)]))
xs32otuRAD<-rad.zipfbrot(t(cat[c(33)]))
xs33otuRAD<-rad.zipfbrot(t(cat[c(34)]))
xs34otuRAD<-rad.zipf(t(cat[c(35)]))
xs35otuRAD<-rad.zipf(t(cat[c(36)]))
xs36otuRAD<-rad.zipf(t(cat[c(37)]))
xs37otuRAD<-rad.zipf(t(cat[c(38)]))
xs38otuRAD<-rad.zipfbrot(t(cat[c(39)]))
xs39otuRAD<-rad.zipf(t(cat[c(40)]))
xs40otuRAD<-rad.zipfbrot(t(cat[c(41)]))
xs41otuRAD<-rad.zipf(t(cat[c(42)]))
######
#Plot everything together
#Order: Preempt, zipf,zipf-mandel,lognormal,null
#Plot everything together
max(xs14otuRAD$y) #16039

png("./Figures/AS_16splots/xoturad_grouped.png",width=8,height=8,units="in",res=300)
par(mfrow=c(3,4))
#For sp 1
plot(xs44otuRAD,col=6,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Preemption",cex=1)
plot(xs1otuRAD,pch=1,lty=1,col=2,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Zipf",cex=1)
points(xs2otuRAD,col=2)
lines(xs2otuRAD,col=2)
points(xs4otuRAD,col=2)
lines(xs4otuRAD,col=2)
points(xs5otuRAD,col=2)
lines(xs5otuRAD,col=2)
points(xs6otuRAD,col=2)
lines(xs6otuRAD,col=2)
points(xs7otuRAD,col=2)
lines(xs7otuRAD,col=2)
points(xs8otuRAD,col=2)
lines(xs8otuRAD,col=2)
points(xs9otuRAD,col=2)
lines(xs9otuRAD,col=2)
points(xs10otuRAD,col=2)
lines(xs10otuRAD,col=2)

plot(xs3otuRAD,pch=1,lty=1,col=7,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Zipf-mandelbrot",cex=1)

plot(xs11otuRAD,pch=1,lty=1,col=11,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Lognormal",cex=1)

#For sp 2 #phylumRADPlot
plot(xs17otuRAD,col=6,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Preemption",cex=1)

plot(xs16otuRAD,pch=1,lty=1,col=2,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Zipf",cex=1)
points(xs20otuRAD,col=2)
lines(xs20otuRAD,col=2)
points(xs21otuRAD,col=2)
lines(xs21otuRAD,col=2)
points(xs22otuRAD,col=2)
lines(xs22otuRAD,col=2)

plot(xs12otuRAD,pch=1,lty=1,col=7,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(xs13otuRAD,col=7)
lines(xs13otuRAD,col=7)
points(xs19otuRAD,col=7)
lines(xs19otuRAD,col=7)
points(xs43otuRAD,col=7)
lines(xs43otuRAD,col=7)

plot(xs15otuRAD,pch=1,lty=1,col=11,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Lognormal",cex=1)
points(xs18otuRAD,col=11)
lines(xs18otuRAD,col=11)
points(xs42otuRAD,col=11)
lines(xs42otuRAD,col=11)
#########################
#For species 3

plot(xs23otuRAD,pch=1,lty=1,col=2,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Zipf",cex=1)
points(xs25otuRAD,col=2)
lines(xs25otuRAD,col=2)
points(xs26otuRAD,col=2)
lines(xs26otuRAD,col=2)
points(xs29otuRAD,col=2)
lines(xs29otuRAD,col=2)
points(xs30otuRAD,col=2)
lines(xs30otuRAD,col=2)
points(xs31otuRAD,col=2)
lines(xs31otuRAD,col=2)
points(xs34otuRAD,col=2)
lines(xs34otuRAD,col=2)
points(xs35otuRAD,col=2)
lines(xs35otuRAD,col=2)
points(xs36otuRAD,col=2)
lines(xs36otuRAD,col=2)
points(xs37otuRAD,col=2)
lines(xs37otuRAD,col=2)
points(xs39otuRAD,col=2)
lines(xs39otuRAD,col=2)
points(xs41otuRAD,col=2)
lines(xs41otuRAD,col=2)


plot(xs24otuRAD,col=7,ylim=c(1,16100),xlim=c(0,310))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(xs27otuRAD,col=7)
lines(xs27otuRAD,col=7)
points(xs28otuRAD,col=7)
lines(xs28otuRAD,col=7)
points(xs32otuRAD,col=7)
lines(xs32otuRAD,col=7)
points(xs33otuRAD,col=7)
lines(xs33otuRAD,col=7)
points(xs38otuRAD,col=7)
lines(xs38otuRAD,col=7)
points(xs40otuRAD,col=7)
lines(xs40otuRAD,col=7)

plot.new()
plot.new()
dev.off()

#head(cat)
xotuRAD<-radfit(t(cat[,c(2:14,16:45)]))
xotuRAD_test<-rad_test(xotuRAD, conf.level = 0.95, log = T )
summary.rad.htest(xotuRAD_test)

png("Figures/AS_16splots/allOTU16radtest.png",height=20,width=30,units="cm",res=300)
par(mfrow=c(1,2),mai=c(1,1,1,1))
plot.rad.test(xotuRAD_test)
plot.rad.test(xotuRAD_test, "post.hoc")
dev.off()

