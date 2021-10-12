rm(list=ls())
library(fields); library(maps); library(ncdf4);library(abind)
library(geosphere);library(mapdata)
setwd("~/Documents")
source('Function/get_geo.R')
source('Function/get_met.R')
source('Function/read_met.R')
source('Function/read_method.R')

plot.box=function(loc,col=1,lwd=0.6,lty=1){
arrows(x0=loc[1],y0=loc[3],x1=loc[2],code=0,col=col,lwd= lwd,lty=lty)
arrows(x0=loc[1],y0=loc[4],x1=loc[2],code=0,col=col,lwd= lwd,lty=lty)
arrows(x0=loc[1],y0=loc[3],y1=loc[4],code=0,col=col,lwd= lwd,lty=lty)
arrows(x0=loc[2],y0=loc[3],y1=loc[4],code=0,col=col,lwd= lwd,lty=lty)
}


plot.field=function(spdata, lon.map, lat.map, type=NULL, same=FALSE, zlim=NULL, col=NULL, nlevel=32, mai=c(0.2, 0.2, 0.2, 0.2), mgp=c(1.4, 0.5, 0), tcl=-0.3, ps=12, legend.mar=3, legend.width=1.2, xaxt="n", yaxt="n", map.region='world', Pacific.centric=FALSE, custom.breaks=NULL, map.xlim=NULL, map.ylim=NULL) {
   
   # This function plots spatial field data using function "image.plot" from package "fields".
   # Packages "fields" and "maps" must be pre-loaded.
   # "spdata" is a matrix or an array of spatial data.
   # "type" is a vector of intended types of presentation, as explained below.
   # Five types of presentation are supported:
   # 1. "sign": variable that has both positive and negative values, good for comparing signs of correlations/effects/deviations.
   # 2. "frac": variable that is a fraction or percentage, good for comparing proportions.
   # 3. "abs": variable that has absolute values only, good for comparing magtitudes.
   # 4. "def": user-defined scale and z-limits. (If chosen, must define the same same scale/limits for all plots.) 
   # 5. Default: no specification for display.
   # If you want all plots to have the same type, simply enter a string scalar for "type".
   
	if (lat.map[length(lat.map)] < lat.map[1]) {
		# 'lat.map' is in decreasing order. Need to reverse it to proceed further.
		lat.map = rev(lat.map)
		spdata = spdata[,length(lat.map):1]
	}	
	par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
	if (is.null(custom.breaks)) {
		if (is.null(type)) {
			zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
		} else if (type == 'sign') {
			if (is.null(zlim)) zlim = c(-max(abs(na.omit(as.vector(spdata)))), max(abs(na.omit(as.vector(spdata))))) else zlim = zlim
			if (is.null(col)) col = rwb.colors(nlevel)
		} else if (type == 'frac') {
			zlim = c(0,1)
		} else if (type == 'abs') {
			zlim = c(0, max(na.omit(as.vector(spdata))))
		} else if (type == 'def') {
			zlim = zlim
		} else {
			zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
		}
		if (is.null(col)) col = tim.colors(nlevel)
		spdata[which(spdata > zlim[2])] = zlim[2]
		spdata[which(spdata < zlim[1])] = zlim[1]
		image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt)
	} else {
		nbreaks = length(custom.breaks)
		nlevel = nbreaks - 1
		if (!is.null(type)) {
			if (type == 'sign') col = rwb.colors(nlevel) else col = tim.colors(nlevel)
		} else col = tim.colors(nlevel)
		zval.breaks = 0:nlevel
		zval.center = 0:(nlevel - 1) + 0.5
		spdata.new = spdata
		for (n in 1:nlevel) spdata.new[which(spdata >= custom.breaks[n] & spdata <= custom.breaks[n + 1])] = zval.center[n]
		spdata.new[which(spdata < custom.breaks[1])] = zval.center[1]
		spdata.new[which(spdata > tail(custom.breaks, 1))] = tail(zval.center, 1)
		spdata = spdata.new
		zlim = c(0, nlevel)
		image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt, breaks=zval.breaks, lab.breaks=custom.breaks)
	}
	if (Pacific.centric) map('world2', add=TRUE, xlim=map.xlim, ylim=map.ylim,col=1,lwd=0.5) else map(map.region, add=TRUE, xlim=map.xlim, ylim=map.ylim,col=1,lwd=0.5)
}

plot_canada=function(col=8, lwd=0.2){
	ss=load("/Users/lu/Documents/CH4_NorthAmerica/Figures/canada_states/output/Canada_states.Rdata")
	for(k in 1:length(canada_states)){
		lines(canada_states[[k]][,1], canada_states[[k]][,2],col=col,lwd=lwd)
	}
}

#======================
setwd("/Users/lu/Documents/CH4_NorthAmerica/Figures/Figure 1")

pdf(file="Figure1.pdf",width=6.5,height=2.5)

# dev.new(width=6.5,height=2.5)
xlim=c(-130,-64);ylim=c(23,60)

par(mar=c(3,2,1,3))
par(mfrow=c(1,2))

ss=load("./fig1a.Rdata")
ind1=(fig1a.lon>=xlim[1] & fig1a.lon <=xlim[2])
ind2=(fig1a.lat>=ylim[1] & fig1a.lat <=ylim[2])

plot.field(fig1a.data[ind1,ind2], fig1a.lon[ind1], fig1a.lat[ind2],type="def",zlim=c(1830,1890),legend.mar=6,mai=c(0.3,0.3,0.5,0.2),ps=10,tcl=-0.2)
axis(1, at=c(-130,-110,-90,-70),labels=c("-130°", "-110°","-90°","-70°"), padj=-0.7,cex=0.8)
axis(2, at=c(25,35,45,55),labels=c("25°", "35°","45°","55°"), hadj=0.5, padj=+0.7,cex=0.8)
title("(a) TROPOMI XCH4,  elevation corrected",font.main=1,cex.main=0.9,line=+0.3)
text(x=-60,y=60,"ppb",xpd=NA,cex=0.9)
# text(x=-120,y=25,"May 2018 - Dec 2019",cex=0.8)
map("state",add=T,lwd=0.2,col=8)
plot_canada(col=8,lwd=0.1)
loc1=c(-124,-100,46,58);plot.box(loc1)
loc2=c(-113,-103,36,46);plot.box(loc2)
loc3=c(-105,-86,25,43);plot.box(loc3)
loc4=c(-86,-71,36,46);plot.box(loc4)
loc5=c(-124,-116,33,40);plot.box(loc5)

ss=load("./fig1b.Rdata")
ind1=(fig1b.lon>=xlim[1] & fig1b.lon <=xlim[2])
ind2=(fig1b.lat>=ylim[1] & fig1b.lat <=ylim[2])
area=cal.area(fig1b.lon, fig1b.lat)
ap= fig1b.data*1e9/1000/(area/1e6)/1000 #tons km-2
ap[ap<=0.1]=NA
plot.field(log10(ap[ind1,ind2]), fig1b.lon[ind1], fig1b.lat[ind2],legend.mar=5,type="def",zlim=c(-1,2),mai=c(0.3,0.3,0.5,0.2),ps=10,tcl=-0.2)
axis(1, at=c(-130,-110,-90,-70),labels=c("-130°", "-110°","-90°","-70°"), padj=-0.7,cex=0.8)
axis(2, at=c(25,35,45,55),labels=c("25°", "35°","45°","55°"), hadj=0.5, padj=+0.7,cex=0.8)
title("(b) Gridded EPA and ECCC oil/gas emissions in 2018",cex.main=0.9,font.main=1,line=+0.3)
map("state",add=T,lwd=0.2,col=8)
plot_canada(col=8, lwd=0.1)

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
cols=c(1,1,1,1)
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
# legend("topright",c("Alvarez 2018","US EIA ","Previous studies"),col=cols, text.col=cols,bty="n",lty=1,cex=0.7,lwd=1.0)
mtext("TROPOMI methane data and oil/gas emissions in the US and Canada",side=3,outer=T,line=-1.5)

dev.off()