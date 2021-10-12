rm(list=ls())
setwd("~/Documents")
library(fields); library(maps);library(ncdf4);library(abind)
library(mapdata);library(geosphere)
library(usmap);library(rgdal)
library(PBSmapping)
library(RColorBrewer)

source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

#=======================================
#========================================
ss0=load("/Users/lu/Documents/CH4_NorthAmerica/Figures/mask.Rdata")

ss=load("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors_uncertainty/All_posterior_bootstrap.Rdata")

all_regions=NULL
for(ifile in 1:dim(all_posteriors)[3]){
OG_posterior_final= all_posteriors[,,ifile]
temp=NULL	
emis.lon=GC.lon;emis.lat=GC.lat
dlon=diff(emis.lon)[1];dlat=diff(emis.lat)[1]

ind=(region_mask>=1 & region_mask<=2)
temp=abind(temp, 'US'=sum(OG_posterior_final[ind],na.rm=TRUE)/1000+0.131)

ind1=(GC.lon>=-105 & GC.lon <=-101)
ind2=(GC.lat>=30 & GC.lat <=34)
temp=abind(temp,'Permian'=sum(OG_posterior_final[ind1,ind2],na.rm=TRUE)/1000)

ind=(region_mask==3)
temp=abind(temp, 'Canada'=sum(OG_posterior_final[ind],na.rm=TRUE)/1000+0.032)

all_regions =cbind(all_regions , temp)
}

apply(all_regions,c(1),mean)
2*apply(all_regions,c(1),sd)
#US 12.7±2.1
#Permian 3.6±0.3
#Canada 2.1±0.6
#========

large_regions =  all_regions[1:3, ]

setwd("~/Documents/CH4_NorthAmerica/Figures/Figure 3")

pdf("Figure_3.pdf", width=8,height=2.4)
cols=c(1, 2, "#048a0a", "#3356FF", "#02A0B7", "#780FF8", 8,  6, "#B78202")

par(mar=c(3,2,1,3))
par(mfrow=c(1,3))
par(mai=c(0.3, 0.45, 0.20, 0.1), mgp=c(1.6, 0.4, 0), tcl=-0.2, ps=17)

bottom_up = c("EPA (2018)"=7.0, "EDGAR (2018)"=7.45, "Alvarez (2015)"=13.0)
top_down  = c("Maasakkers (2000-2015)" = 11.1+0.1, "Lu (2010-2017)"=14.5+0.1)
boxplot(large_regions[1,]+0.1,ylim=c(0,16),xlim=c(-0.2, 2),ylab=bquote("Emissions (Tg "*a^-1*")"))
points(rep(0.5, length(bottom_up)), bottom_up,pch=c(17,24,24),col=cols[c(1,2,3)],cex=c(1.8,1.5,1.5),lwd=2)

points(rep(1.05, length(top_down)), top_down,pch=25,col=cols[c(4,5)],cex=1.5,lwd=2)
mtext("Bottom-up",side=1,line=0.5,adj=0.12,cex=0.6)
mtext("Top-down",side=1,line=0.5,adj=0.7,cex=0.6)
title("(a) United States", font.main=1, cex.main=1)
# legend("bottomright",legend=bquote("Mean = 12.5±2.1 Tg "*a^-1*" ") ,bty="n",cex=0.9)
text(x=0.45,y=6.6,"EPA",pos=2,cex=0.8,col=1)
text(x=0.45,y=7.8,"EDGAR",pos=2,cex=0.8,col=cols[2])
text(x=0.45,y=13.0,"Alvarez",pos=2,cex=0.8,col=cols[3])
text(x=1.2,y=12.5+0.1,"This work",pos=4,cex=0.8,col=1,font=2)
text(x=1.2,y=10.7+0.1,"Maasakkers",pos=4,cex=0.8,col=cols[4])
text(x=1.2,y=14.2+0.1,"Lu",pos=4,cex=0.8,col=cols[5])

bottom_up = c("ECCC (2018)"=1.6, "EDGAR (2018)"=1.7)
top_down  = c("Maasakkers (2000-2015)"=2.3,"Lu (2010-2017)"=2.9, 
"Baray (2010-2015)"=3.6, "Chan (2010-2017)"=3.0)
boxplot(large_regions[3,],ylim=c(0,4),xlim=c(-0.2, 2),ylab=bquote("Emissions (Tg "*a^-1*")"))
points(rep(0.5, length(bottom_up)), bottom_up,pch=c(17,24,24),col=cols[c(6,2,7)],cex=c(1.8,1.5,1.5),lwd=2)
points(rep(1.05, length(top_down)), top_down,pch=25,col=cols[c(4,5,8,9)],cex=1.5,lwd=2)
text(x=0.45,y=1.5,"ECCC",pos=2,cex=0.8,col=cols[6])
text(x=0.45,y=1.8,"EDGAR",pos=2,cex=0.8,col=cols[2])
text(x=1.2,y=2.12,"This work",pos=4,cex=0.8,col=1,font=2)
text(x=1.2,y=2.37,"Maasakkers",pos=4,cex=0.8,col=cols[4])
text(x=1.2,y=2.7,"Lu",pos=4,cex=0.8,col=cols[5])
text(x=1.2,y=3.6,"Baray",pos=4,cex=0.8,col=cols[8])
text(x=1.2,y=3.0,"Chan",pos=4,cex=0.8,col=cols[9])


# legend("topleft",c("EDGAR", "This study", "ICF (2013)"),col=c(4,"#048a0a"),pch=16,bty="n",cex=0.8,text.col=c(4,"#048a0a"),y.intersp=1.5)
mtext("Bottom-up",side=1,line=0.5,adj=0.12,cex=0.6)
mtext("Top-down",side=1,line=0.5,adj=0.7,cex=0.6)
title("(b) Canada", font.main=1, cex.main=1)
# legend("bottomright",legend=bquote("Mean = 2.2±0.6 Tg "*a^-1*" ") ,bty="n",cex=0.9)

par(mai=c(0.2, 0.0, 0.20, 0.1), mgp=c(1.7, 0.4, 0), tcl=-0.2, ps=17)
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axe=F,xlab="",ylab="")
# legend("topleft",c("US EPA  (2018)", "Crippa (2018)", "Alvarez (2015)","Maasakkers (2010-2015)","Lu (2010-2017)","ECCC (2018)","ICF (2013)","Baray (2010-2015)", "Chan (2010-2017)"),col=cols,pch=c(17, 24, 24, 25, 25, 17, 24, 25,25),bty="n",cex=0.9,text.col=cols,y.intersp=1.2)
legend("topleft",c("EPA  (2018)", "EDGAR v6.0 (2018)", "Alvarez (2015)","Maasakkers (2010-2015)","Lu (2010-2017)","ECCC (2018)","Baray (2010-2015)", "Chan (2010-2017)"),col=cols[c(1:6,8:9)],pch=c(17, 24, 24, 25, 25, 17, 25,25),bty="n",cex=0.9,text.col=cols[c(1:6,8:9)],y.intersp=1.2)

legend("bottomleft",c("Prior in this study","Bottom-up estimates","Top-down estimates"),pch=c(17,24,25),bty="n",cex=0.9,text.col=1,y.intersp=1.2)

dev.off()