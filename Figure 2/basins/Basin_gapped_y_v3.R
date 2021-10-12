rm(list=ls())
setwd("~/Documents")
library(fields); library(maps);library(ncdf4);library(abind)
library(mapdata);library(geosphere)
library(usmap);library(rgdal)
library(PBSmapping)
library(RColorBrewer)
library("plotrix")

source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

cal_grid_area=function(MAT, emis.lon, emis.lat){

aircraft <- cbind(rep(1,dim(MAT)[1]), POS=1:dim(MAT)[1],MAT)
colnames(aircraft)=c("PID","POS","X","Y")
aircraft =data.frame(aircraft)

loc_ij=NULL
for(k in 1:dim(aircraft)[1]){
  loc=find.lon.lat(aircraft[k,"X"], aircraft[k,"Y"], emis.lon, emis.lat)
  loc_ij=rbind(loc_ij,loc)
}
min_i=min(loc_ij[,1]);max_i=max(loc_ij[,1])
min_j=min(loc_ij[,2]);max_j=max(loc_ij[,2])

num=1
overlap_gris=list()
for(i in min_i:max_i){
	for(j in min_j:max_j){
    MAT2=cbind(c(emis.lon[i]-dlon/2,emis.lon[i]-dlon/2,emis.lon[i]+dlon/2,emis.lon[i]+dlon/2,emis.lon[i]-dlon/2),c(emis.lat[j]-dlat/2,emis.lat[j]+dlat/2,emis.lat[j]+dlat/2,emis.lat[j]-dlat/2,emis.lat[j]-dlat/2))
    p5 <- cbind(rep(1,dim(MAT2)[1]), POS=1:dim(MAT2)[1], MAT2)
    colnames(p5)=c("PID","POS","X","Y")
    p5=data.frame(p5)
    p6 <- joinPolys(aircraft,p5)
	if(!is.null(p6)){
overlap_gris[[num]]=list(i=i,j=j,center_lon=emis.lon[i],center_lat=emis.lat[j],coords=p6,area_overlap= sum(calcArea(p6)$area), area_grid=calcArea(p5)$area)
	num=num+1
	} 
}
}

return(overlap_gris)	
}

region_mean=function(metric){
fenzi=0;fenmu=0
for(k in 1:length(overlap_gris)){
   met=overlap_gris[[k]]
   if(is.na(metric[met$i,met$j]))next
   fenzi=fenzi+met$area_overlap* metric[met$i,met$j]
   fenmu=fenmu+met$area_overlap
 }
return(fenzi/fenmu)
}

region_sum=function(metric){
fenzi=0
for(k in 1:length(overlap_gris)){
   met=overlap_gris[[k]]
   if(is.na(metric[met$i,met$j]))next   
   fenzi=fenzi+met$area_overlap/met$area_grid*metric[met$i,met$j]
   #if (met$area_overlap/met$area_grid>=0.25)fenzi=fenzi+metric[met$i,met$j]
 }
return(fenzi)
}

region_area=function(){
fenzi=0
for(k in 1:length(overlap_gris)){
   met=overlap_gris[[k]]
   fenzi=fenzi+met$area_overlap
 }
return(fenzi)
}


#=======================================
#========================================
setwd("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors/")
all_priors=NULL;all_scalings=NULL;all_posteriors=NULL;all_AK=NULL
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
OG_prior_final= all_priors[,,1]
OG_scaling_final =apply(all_scalings,c(1,2),mean,na.rm=TRUE)
OG_posterior_final =apply(all_posteriors,c(1,2),mean,na.rm=TRUE)
OG_posterior_error =apply(all_posteriors,c(1,2),sd,na.rm=TRUE)
OG_AK = apply(all_AK,c(1,2),mean,na.rm=TRUE)
emis.lon=GC.lon;emis.lat=GC.lat
dlon=diff(emis.lon)[1];dlat=diff(emis.lat)[1]

all_emis=NULL;all_sites=NULL
setwd("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors/basin_tabulate/combined")
dir_names = list.dirs(path=".", recursive= FALSE, full.names=FALSE)

basin_names=NULL;categories=NULL
for(k in 1:length(dir_names)){
	ap=as.numeric(substr(dir_names[k],1,1))
	ap2=substr(dir_names[k],3,nchar(dir_names[k]))
	categories =abind(categories, ap)
	basin_names=abind(basin_names,ap2)
}

ss=load("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors_uncertainty/All_posterior_bootstrap.Rdata")

ibasin=5
for(ibasin in 1:length(dir_names)){

temp_emis=NULL
files=Sys.glob(paste0(dir_names[ibasin],"/*.csv"))

index=1
for(index in 1:length(files)){
	
filename=files[index]
site_name = basin_names[ibasin]
MAT=read.csv(filename)
if(max(MAT)>=1e3)next
overlap_gris=cal_grid_area(MAT, emis.lon, emis.lat)
ss=load("/Users/lu/Documents/CH4_NorthAmerica/CPU/altitudes_OG1.0/retrieval_metrics.Rdata")
# ss2=load("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors/uncertainty_contribution/uncertainty_ensemble.Rdata")
#===== calulcate emissions =======
# basin_area= region_area()*(100*100)
basin_area= areaPolygon(MAT)/1e6 #square kilometeres
basin_freq=region_sum(obs_freq)
basin_albedo=region_mean(albedo2)
basin_AOT2=region_mean(AOT2)
basin_altitude=region_mean(altitude)
basin_altitude_std=region_mean(altitude_std)
basin_blend_albedo =region_mean(blend_albedo)
basin_bias_std =region_mean(all_bias_std)
basin_prior = region_sum(OG_prior_final)#Gg
basin_posterior = region_sum(OG_posterior_final)#Gg
basin_AK = region_sum(OG_AK)#Gg

all_ensemble=array(NA, 2400)
for(iloop in 1:2400){
fenmu=0	
for(k in 1:length(overlap_gris)){
   met=overlap_gris[[k]]
   fenmu=fenmu+ met$area_overlap/met$area_grid*all_posteriors[met$i, met$j, iloop]
   #if(met$area_overlap/met$area_grid>0.25)fenmu=fenmu+ all_posteriors[met$i, met$j, iloop]
 }
all_ensemble[iloop] =fenmu
 }
basin_relative_uncertainty = sd(all_ensemble)/fenmu
basin_0025 = quantile(all_ensemble, probs=0.025)#Gg
basin_0975 = quantile(all_ensemble, probs=0.975)#Gg

temp_emis =rbind(temp_emis, c(basin_area, basin_freq, basin_albedo, basin_AOT2, basin_altitude, basin_altitude_std, basin_blend_albedo, basin_bias_std, basin_prior, basin_posterior, basin_AK, basin_relative_uncertainty, basin_0025, basin_0975 ))
}

all_emis=rbind(all_emis, apply(temp_emis,c(2),mean))
all_sites =rbind(all_sites, basin_names[ibasin])

}

rownames(all_emis)=all_sites
colnames(all_emis)=c("area","freq","albedo","AOT","altitude","altitude_std","blended_albedo",'bias_std',"prior","posterior","AK","relative_SD","basin_0025","basin_0975")


#========================================================
#========================================================
final = cbind(all_emis[,c("prior", "posterior", "relative_SD","AK","basin_0025","basin_0975")]) 
final[,c("relative_SD")] = 2*final[,c("relative_SD")]

rownames(final)[13]="Marcellus"
rownames(final)[19]="Delaware"

names=rownames(final)
names[names=="NE_PA"]="NE PA"
names[names=="San_Juan"]="San Juan"
names[names=="West_Arkoma"]="West Arkoma"
names[names=="Alberta_East"]="Alberta East"
names[names=="Alberta_West"]="Alberta West"
names[names=="SW_PA"]="SW PA"
rownames(final)=names

setwd("/Users/lu/Documents/CH4_NorthAmerica/Figures/Figure 2/basins")
save(final, file="basin.Rdata")

#=======================================================
pdf(file="./barchart_v2.pdf",width=6.5,height=2.5)

ss=load("./basin.Rdata")
ind=order(final[,"posterior"],decreasing=TRUE)
final=final[ind,]

# rownames(final)[8]="San Juan"
# rownames(final)[13]="Alberta West"
# rownames(final)[14]="West Arkoma"
# rownames(final)[15]="SW PA"
# rownames(final)[19]="NE PA"

figb=t(final[,1:2])/1000
standardErrors= final[,3]
error_0025=final[,"basin_0025"]/1e3
error_0975=final[,"basin_0975"]/1e3

from=0.98;delta=1.0
error_0025[error_0025> from]<-error_0025[error_0025> from]-delta
error_0975[error_0975> from]<-error_0975[error_0975> from]-delta

mean(final[1:9,3])
2*sd(final[1:9,3])
mean(final[10:19,3])
2*sd(final[10:19,3])

# dev.new(width=6.5,height=2.5)
par(mai=c(0.7, 0.7, 0.22, 0.2), mgp=c(1.8, 0.5, 0), tcl=-0.2, ps=10)
figb2= figb
figb2[figb2> from]<-figb2[figb2> from]-delta
barCenters=barplot(height= figb2,main="", ylab="Emissions (Tg)", beside=TRUE,col=rep(c("#E69F00", "#56B4E9"),7),bty="n",ylim=c(0,4.2-delta),border=0,axe=F,xaxt="n",yaxt="n",space=c(0,1.3))
width= barCenters[2,1]-barCenters[1,1]

rect(barCenters[1,1]-0.5*width,0,barCenters[1,1]+0.5*width,2.25-delta,col="#E69F00",border="#E69F00",density = 20)
rect(barCenters[2,1]-0.5*width,0,barCenters[2,1]+0.5*width,3.7-delta,col="#56B4E9",border="#56B4E9",density = 20)

box()
axis(side=2,las=2, at=c(0,0.5,0.93,1.07,1.5,2.0,2.5,3.0),labels=format(c(0,0.5,1.0,1.0+delta,1.5+delta, 2.0+delta,2.5+delta,3.0+delta), digits=1))
axis.break(2, from,style="gap",breakcol="grey80",brw=0.02)
axis.break(2, from*(1+0.02), breakcol="black", style="slash",brw=0.02)
axis.break(4, from*(1+0.02), breakcol="black", style="slash",brw=0.02)

text(x=(barCenters[1,]+barCenters[2,])/2, y=par("usr")[3] - 0.2,labels=colnames(figb),xpd=NA,srt=45,cex=0.9,adj=1.0)
means= figb2[2,]
ap=round(figb[2,],2)
ap= format(ap,digits=1,trim=TRUE, justify="centre")

arrows(barCenters[2,], means*(1+standardErrors), barCenters[2,], means*(1-standardErrors), lwd=1, angle=90, code=3,length=0.02)
arrows(barCenters[2,1]+0*width, 3.7+0.35-delta, barCenters[2,1]+0*width, 3.7-0.35-delta,lwd=1, angle=90, code=3,length=0.02,lty=2)

legend("topright", c("Gridded EPA/ECCC inventories","Posterior"), cex=0.9, bty="n",fill=c("#E69F00","#56B4E9"),border=1,density=c(NA,NA))
title("(c) Prior and posterior emissions in major oil/gas basins",font.main=1,cex.main=1)

dev.off()