#follow from visualizing taxa distribution
##16s tau
View(datatau)
sum(datatau$n)#328 otus
nrow(datatau)#22 phyla
taupercent<-datatau%>%
  mutate(percent=(n/sum(n))*100)
View(taupercent)
print(taupercent)

##
##16s tme
sum(datatme$n)#279 otus
nrow(datatme)#20 phyla
tmepercent<-datatme%>%
  mutate(percent=(n/sum(n))*100)
View(tmepercent)
print(taupercent)

##
##16s tci
View(datatci)
sum(datatci$n)#384 otus
nrow(datatci)#23 phyla
tcipercent<-datatci%>%
  mutate(percent=(n/sum(n))*100)
View(tcipercent)
print(taupercent)

###without unclass
#follow from visualizing taxa distribution
##16s tau
xdatatau<-datatau%>%
  filter(Phylum!="Unclassified")
sum(xdatatau$n)#274 otus
nrow(xdatatau)#21 phyla
xtaupercent<-xdatatau%>%
  mutate(percent=(n/sum(n))*100)
print(xtaupercent)

##
##16s tme
xdatatme<-datatme%>%
  filter(Phylum!="Unclassified")
sum(xdatatme$n)#229 otus
nrow(xdatatme)#19 phyla
xtmepercent<-xdatatme%>%
  mutate(percent=(n/sum(n))*100)
print(xtmepercent)

##
##16s tci
xdatatci<-datatci%>%
  filter(Phylum!="Unclassified")
sum(xdatatci$n)#317 otus
nrow(xdatatci)#22 phyla
xtcipercent<-xdatatci%>%
  mutate(percent=(n/sum(n))*100)
print(xtcipercent)

head(tcix)
filt_tci<-tcix%>%
  filter(Phylum=="Proteobacteria")
View(filt_tci)
filt_tci%>%
  count(Class=="Gammaproteobacteria")
#64 alpha
#86 gamma
filt_tci<-tcix%>%
  filter(Phylum=="Proteobacteria")