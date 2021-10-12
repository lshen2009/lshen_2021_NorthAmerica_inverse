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
setwd("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors/")

all_regions=NULL

files=abind(Sys.glob("./Lcurve_noAdjBias/*OG1.0*.Rdata"),Sys.glob("./Lcurve_noAdjBias2/*OG1.0*.Rdata"),Sys.glob("./Lcurve_noAdjBias3/*OG1.0*.Rdata"))
allfiles=abind(Sys.glob("./Lcurve_noAdjBias/*.Rdata"),Sys.glob("./Lcurve_noAdjBias2/*.Rdata"),Sys.glob("./Lcurve_noAdjBias3/*.Rdata"))



for(ifile in 1:length(allfiles)){
setwd("/Users/lu/Documents/CH4_NorthAmerica/CPU/All_posteriors/")
temp=NULL	
print(ifile)	
met=load(allfiles[ifile])
OG_AK = OG_AK_final
emis.lon=GC.lon;emis.lat=GC.lat
dlon=diff(emis.lon)[1];dlat=diff(emis.lat)[1]

ind=(GC.lat<=49.1)
temp=abind(temp, 'US'=sum(OG_posterior_final[,ind],na.rm=TRUE)/1000)

ind1=(GC.lon>=-105 & GC.lon <=-101)
ind2=(GC.lat>=30 & GC.lat <=34)
temp=abind(temp,'Permian'=sum(OG_posterior_final[ind1,ind2],na.rm=TRUE)/1000)

ind=(GC.lat>49.1)
temp=abind(temp, 'Canada'=sum(OG_posterior_final[,ind],na.rm=TRUE)/1000)

#--- read each basin ----
all_emis=NULL
all_sites=NULL
setwd("/Users/lu/Documents/CH4_NorthAmerica/paper/Field_measurements/Alvarez_2018/combined")
files=Sys.glob("*.csv")

index=1
for(index in 1:length(files)){

filename=files[index]
site_name =strsplit(filename,"[_|.]")[[1]][1]
MAT=read.csv(filename)
if(max(MAT)>=1e3)next
#--- calculate overlap_grids -----
aircraft <- cbind(rep(1,dim(MAT)[1]), POS=1:dim(MAT)[1],MAT)
colnames(aircraft)=c("PID","POS","X","Y")
aircraft =data.frame(aircraft)
calcArea(aircraft)$area

loc_ij=NULL
for(k in 1:dim(aircraft)[1]){
  loc=find.lon.lat(aircraft[k,"X"], aircraft[k,"Y"],emis.lon,emis.lat)
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
overlap_gris[[num]]=list(i=i,j=j,center_lon=emis.lon[i],center_lat=emis.lat[j],coords=p6,area_overlap= calcArea(p6)$area, area_grid=calcArea(p5)$area)
	num=num+1
	} 
}
}

#===== calulcate emissions =======
 emis_orig=0;emis_inv1=0;emis_inv1_error=0
 fenmu=0
 for(k in 1:length(overlap_gris)){
   met=overlap_gris[[k]]
   emis=met$area_overlap/met$area_grid*OG_prior_final[met$i,met$j]
   emis_orig= emis_orig + emis
   emis=met$area_overlap/met$area_grid*OG_posterior_final[met$i,met$j]
   emis_inv1 = emis_inv1+ emis
 }

#================================
 all_emis=rbind(all_emis, c(emis_orig, emis_inv1))
 all_sites=abind(all_sites, site_name)
}

all_sites[1:7]="Barnett"
rownames(all_emis)=all_sites
colnames(all_emis)=c("prior","posterior")

basins=tapply(all_emis[,2],all_sites,mean)/1000

temp=abind(temp, basins)

all_regions =cbind(all_regions , temp)
}

apply(all_regions,c(1),sd)/apply(all_regions,c(1),mean)
#========

large_regions =  all_regions[1:3, ]
Others = all_regions[4:15, ]

setwd("/Users/lu/Documents/CH4_NorthAmerica/Figures/Figure 4")

pdf("Part_a.pdf", width=4,height=2.8)
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.4, 0.4, 0), tcl=-0.2, ps=11)
boxplot(t(large_regions),ylim=c(0,16),ylab="Emissions (Tg a-1)")
dev.off()

rownames(Others)=c("Alberta", "Bakken", "Barnett","DJ", "Faye" , "Hay","NE PA","Permian","sj" ,"SW PA","UI", "West Arkoma" )
pdf("Part_b.pdf", width=12,height=2.8)
par(mai=c(0.3, 0.6, 0.25, 0.2), mgp=c(1.4, 0.4, 0), tcl=-0.2, ps=11)
boxplot(t(Others),ylim=c(0,1),ylab="Emissions (Tg a-1)")
dev.off()