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
   #if(met$area_overlap/met$area_grid>=0.25) fenzi=fenzi+metric[met$i,met$j]
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
OG_prior_final=apply(all_priors[,,1:12],c(1,2),mean,na.rm=TRUE)
OG_scaling_final =apply(all_scalings,c(1,2),mean,na.rm=TRUE)
OG_posterior_final =apply(all_posteriors,c(1,2),mean,na.rm=TRUE)
OG_posterior_error =apply(all_posteriors,c(1,2),sd,na.rm=TRUE)
OG_AK = apply(all_AK,c(1,2),mean,na.rm=TRUE)
emis.lon=GC.lon;emis.lat=GC.lat
dlon=diff(emis.lon)[1];dlat=diff(emis.lat)[1]

all_emis=NULL;all_sites=NULL
setwd("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors/basin_tabulate/combined")
dir_names = list.dirs(path=".", recursive= FALSE, full.names=FALSE)
dir_names = dir_names[c(1:9, 15:19)]

basin_names=NULL;categories=NULL
for(k in 1:length(dir_names)){
	ap=as.numeric(substr(dir_names[k],1,1))
	ap2=substr(dir_names[k],3,nchar(dir_names[k]))
	categories =abind(categories, ap)
	basin_names=abind(basin_names,ap2)
}

ss=load("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors_uncertainty/All_posterior_bootstrap.Rdata")

ibasin=2
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
   fenmu=fenmu+met$area_overlap/met$area_grid* all_posteriors[met$i, met$j, iloop]
 }
all_ensemble[iloop] =fenmu
 }
basin_relative_uncertainty = sd(all_ensemble)/fenmu
 
temp_emis =rbind(temp_emis, c(basin_area, basin_freq, basin_albedo, basin_AOT2, basin_altitude, basin_altitude_std, basin_blend_albedo, basin_bias_std, basin_prior, basin_posterior, basin_AK, basin_relative_uncertainty ))
}

all_emis=rbind(all_emis, apply(temp_emis,c(2),mean))
all_sites =rbind(all_sites, basin_names[ibasin])

}


#==========
rownames(all_emis)=all_sites
colnames(all_emis)=c("area","freq","albedo","AOT","altitude","altitude_std","blended_albedo",'bias_std',"prior","posterior","AK","relative_SD")

Alvarez_emis=c("Bakken"=27, "Barnett"=60,
"DJ"=19,"Fayetteville"=31,"Haynesville"=73,"NE PA"=18,
"San Juan"=57*1,"Uinta"=55*1,"West Arkoma"=26,
"Alberta_East"=6.7,"Alberta_West"=4.8, "Permian"=304, "SW_PA"=1.3, "Permian_sub"=176
)*24*365/1000

Alvarez_error=c("Bakken"=13,"Barnett"=11,
"DJ"=14,"Fayetteville"=30,"Haynesville"=54,"NE PA"=14,
"San Juan"=54,"Uinta"=31,"West Arkoma"=30,
"Alberta_East"=6.7,"Alberta_West"=6.4,"Permian"=109,  "SW_PA"=1.3, "Permian_sub" = 34
)*24*365/1000/2 #average relative 1xsd is 36%

final = cbind(all_emis[,c("prior", "posterior", "relative_SD")], Alvarez_emis, Alvarez_error) 
final[,c("relative_SD", "Alvarez_error")]=2*final[,c("relative_SD", "Alvarez_error")]#I am using 2*SD
final[12,"posterior"]= 2613 #whole Permian, part1
final[14,"posterior"]= 847.2#100x100 km Permian, part2

# final[12,"posterior"]= 2537 #whole Permian, part1
# final[14,"posterior"]= 718.7#100x100 km Permian, part2

final[12,"posterior"]= 2537 #whole Permian, part1
final[14,"posterior"]= 800#100x100 km Permian, part2
final[,c("prior","posterior","Alvarez_emis","Alvarez_error")]=final[,c("prior","posterior","Alvarez_emis","Alvarez_error")]/1000


cols=brewer.pal(9,"Paired")
cols[1]= "#5a5a5a"
cols[3]="#861FDE"
cols[7]="#80391e"
ap=cols[4];cols[4]=cols[5];cols[5]=ap
# cols=abind(rev(cols), cols[c(1,4,6,8)])
cols=abind(cols, cols[c(1,4,6,8,3)])
pchs=c(rep(16,9), rep(17,5))

legend_names= c("Bakken (2014)", "Barnett (2013)", "DJ (2012)", "Fayetteville (2015)", "Haynesville (2013)", "NE PA (2015)", "San Juan (2015)", "Uinta (2012)", "West Arkoma (2013)", "Alberta_East (2016)", "Alberta_West (2016)", "Permian (2018)", "SW PA (2015)", "Delaware (2020)")

# dev.new(width=8, height=3.5)
pdf("/Users/lu/Documents/CH4_NorthAmerica/Figures/Figure_4/Rplot_PermianX1.pdf", width=8, height=3.5)

# dev.new(width=8, height=3.5)
par(mar=c(3,2,1,3))
par(mfrow=c(1,2))

par(mai=c(0.5, 0.6, 0.25, 0.2), mgp=c(1.4, 0.4, 0), tcl=-0.2, ps=11)

plot(final[,"Alvarez_emis"],final[,"prior"],xlim=c(0,4),ylim=c(0,4),pch=pchs,xlab=bquote("Field campaign estimates (Tg "*a^-1*")"),ylab=bquote("Gridded national inventories (Tg "*a^-1*")"),col=cols,cex=1.0)
abline(a=0,b=1,lty=2)
arrows(final[,"Alvarez_emis"]-final[,"Alvarez_error"], final[,"prior"], final[,"Alvarez_emis"]+final[,"Alvarez_error"], final[,"prior"],length=0.03,angle=90,code=3,col=cols)
legend("topleft", legend_names,col=cols,text.col=cols,bty="n",lty=1,pch=pchs,cex=0.8,y.intersp=0.8)
title("Gridded EPA and ECCC vs. field estimates", font.main=1, cex.main=1)
legend("topright", legend=bquote(""*R^2*"=0.73"),text.col=1,bg="white")

plot(final[,"Alvarez_emis"],final[,"posterior"],xlim=c(0,4),ylim=c(0,4),pch=pchs,xlab=bquote("Field campaign estimates (Tg "*a^-1*")"),ylab=bquote("TROPOMI inversion estimates (Tg "*a^-1*")"),col=cols,cex=1.0)
abline(a=0,b=1,lty=2)
arrows(final[,"Alvarez_emis"]-final[,"Alvarez_error"], final[,"posterior"], final[,"Alvarez_emis"]+final[,"Alvarez_error"], final[,"posterior"],length=0.03,angle=90,code=3,col=cols)
arrows(final[,"Alvarez_emis"], final[,"posterior"]*(1-1*final[,"relative_SD"]), final[,"Alvarez_emis"], final[,"posterior"]*(1+1*final[,"relative_SD"]),length=0.03,angle=90,code=3,col=cols )

title("TROPOMI inversion vs. field estimates", font.main=1, cex.main=1)
legend("topright", legend=bquote(""*R^2*"=0.92"),text.col=1,bg="white")

dd=0.25
arrows(2.95-dd,0.9,3.25-dd,0.9,length=0.03,angle=90,code=3)
arrows(3.10-dd,0.15,3.10-dd,0.45,length=0.03,angle=90,code=3)

dev.off()

cor(final[,"Alvarez_emis"], final[,"posterior"])^2
cor(final[,"Alvarez_emis"], final[,"prior"])^2