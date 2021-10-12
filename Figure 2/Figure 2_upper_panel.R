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

log_10=function(x){
	direction=sign(x)
	x2=abs(x)
	y=direction*log10(x2)
	y[x==0]=NA
	return(y)		
}

plot_canada=function(col=8, lwd=0.2){
	ss=load("/Users/lu/Documents/CH4_NorthAmerica/Figures/canada_states/output/Canada_states.Rdata")
	for(k in 1:length(canada_states)){
		lines(canada_states[[k]][,1], canada_states[[k]][,2],col=col,lwd=lwd)
	}
}

#=======================================
#========================================
setwd("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors/")
all_priors=NULL
all_scalings=NULL
all_posteriors=NULL
all_AK=NULL
files1=abind(Sys.glob("./Lcurve_noAdjBias2/*OG1.0*.Rdata"),Sys.glob("./Lcurve_noAdjBias3/*OG1.0*.Rdata"))
files3=abind(Sys.glob("./Lcurve_noAdjBias2/*OG1.5_v2.Rdata"),Sys.glob("./Lcurve_noAdjBias3/*OG1.5_v2.Rdata"))

files=abind(files1, files3)

for(ifile in 1:length(files)){
	met=load(files[ifile])
	all_priors=abind(all_priors, OG_prior_final, along=3)
	all_scalings =abind(all_scalings, OG_scaling_final, along=3)
	all_posteriors =abind(all_posteriors, OG_posterior_final, along=3)
	all_AK=abind(all_AK, OG_AK_final, along=3)
}
OG_prior_final=apply(all_priors,c(1,2),mean,na.rm=TRUE)
OG_scaling_final =apply(all_scalings,c(1,2),mean,na.rm=TRUE)
OG_posterior_final =apply(all_posteriors,c(1,2),mean,na.rm=TRUE)
OG_posterior_error =apply(all_posteriors,c(1,2),sd,na.rm=TRUE)
OG_AK = apply(all_AK,c(1,2),mean,na.rm=TRUE)
emis.lon=GC.lon;emis.lat=GC.lat
dlon=diff(emis.lon)[1];dlat=diff(emis.lat)[1]

ss=load(files[1])
prior_OG10= OG_prior_final


ss=load("/Users/lu/Documents/CH4_NorthAmerica/Data/Emissions_multiple_souces/mask/US_mask.Rdata")
US_mask=sp.dissolve(mask, mask.lon, mask.lat, GC.lon, GC.lat)
US_mask[is.na(US_mask)]=0
US_mask[,GC.lat>=40]=1
# plot.field(US_mask,GC.lon,GC.lat)

#figure_1a mask
ss=load("/Users/lu/Documents/CH4_NorthAmerica/Figures/Figure 1/fig1b.Rdata")
area=cal.area(fig1b.lon, fig1b.lat)
fig1b_emis= fig1b.data*1e9/1000/(area/1e6)/1000 #tons km-2

setwd("/Users/lu/Documents/CH4_NorthAmerica/Figures/Figure 2")
pdf("Figure_2_upper.pdf",width=6.5,height=2.2)

# dev.new(width=6.5,height=2.2)
par(mar=c(3,2,1,3))
par(mfrow=c(1,2))

xlim=c(-125,-65);ylim=c(25,58)
mai=c(0.05,0.1,0.2,0.2);ps=11
ind1=(GC.lon>=xlim[1] & GC.lon<=xlim[2])
ind2=(GC.lat>=ylim[1] & GC.lat<=ylim[2])

OG_scale= OG_posterior_final/prior_OG10
OG_scale[US_mask<=0.5]=NA
OG_scale[fig1b_emis<=0.1]=NA
plot.field(OG_scale[ind1,ind2],GC.lon[ind1],GC.lat[ind2],type="sign",zlim=c(-1,3),legend.mar=7, ps= ps, mai=mai,col=PN.cols(32))
map("state",add=TRUE,lwd=0.4,col=8)
plot_canada(col=8,lwd=0.4)
box()

setwd("/Users/lu/Documents/CH4_NorthAmerica/CPU/Basins/combined")
dir_names = list.dirs(path=".", recursive= FALSE, full.names=FALSE)
basin_names=NULL
categories=NULL
for(k in 1:length(dir_names)){
	ap=as.numeric(substr(dir_names[k],1,1))
	ap2=substr(dir_names[k],3,nchar(dir_names[k]))
	categories =abind(categories, ap)
	basin_names=abind(basin_names,ap2)
}
cols=c(1,2,4)
cols=c(4,4,4,4)
cols=c(2,2,2,2)
for(ifile in 1:length(dir_names)){
	files=Sys.glob(paste0(dir_names[ifile],"/*.csv"))
	NN=length(files)
	
	if(NN==1){		
		met=read.csv(files)
		lines(met[,1],met[,2],col= cols[categories[ifile]],lwd=0.5)
	}
	
	if(NN>1){	
		areas=array(NA,NN)	
		for(k in 1:length(files)){
			met=read.csv(files[k])
			areas[k]=areaPolygon(met)/1e6			
		}
		k=which.max(areas)
		met=read.csv(files[k])
		lines(met[,1],met[,2],col= cols[categories[ifile]],lwd=0.5)
	}	
	
}

title("(a) Correction factors", font.main=1, cex.main=1 )

# x=seq(-1,3, l=64)
# cols=PN.cols(length(x))
# ind=(x>=0)
# image.plot(legend.only=TRUE, col=cols[ind],zlim=c(0,3))

area=cal.area(GC.lon, GC.lat)
ap= OG_posterior_final*1e9/1000/(area/1e6)/1000 #kg km-2
# ap[ap<=0.1]=NA
ap[fig1b_emis<=0.1]=NA
plot.field(log10(ap[ind1,ind2]), GC.lon[ind1], GC.lat[ind2],type="def",zlim=c(-1,2),legend.mar=7, ps= ps, mai=mai)
map("state",add=TRUE,lwd=0.4,col=8)
plot_canada(col=8,lwd=0.4)
box()
title("(b) Posterior emissions", font.main=1, cex.main=1)
dev.off()
