go.pdf <- F
#copy default par()
par.def <- par()      

library("rgdal")
library("maptools")
library("mapplots")
library("raster")
library("plotrix")
# Generate a lancover map  overlap with TC tracks and effect size 
wrk.dir   <- c("/lfs/home/ychen/LAI_STUDY_EAsia/")
track.dir <- c("/TRACK_DATA/")
lc.dir    <- c("/LANDCOVER_DATA/")
table.dir <- c("/work/ychen/scripts/R/Rscripts/LAI_STUDY_EA/Rda_50/")

# load the land cover map 
load("/lfs/home/ychen/LAI_STUDY_EAsia/LANDCOVER_DATA/2015.esa.landcover.east.asia.rda")

#convert arrary to raster obj for ploting
xymn <- c(90, 150, 0, 60)

#set NA = 999 for the ocean 
esa.lc[is.na(esa.lc)] <- 250
#apply the tc occurence map as mask
#esa.lc[arr.for.plot < 0.05] <- 0 

lc.raster <- raster( x=t(esa.lc),
                         xmn=xymn[1], xmx=xymn[2], ymn=xymn[3], ymx=xymn[4],
                         crs=CRS("+proj=longlat +datum=WGS84"))

# variable name: esa.lc  
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30 
# set water(210) to NA
esa.lc[ (esa.lc == 210 ) ] <- NA 
esa.lc[ (esa.lc >= 10)  & (esa.lc <= 40 ) ] <- 1
esa.lc[ (esa.lc >= 50)  & (esa.lc <= 120) ] <- 2
#esa.lc[ (esa.lc == 160) | (esa.lc == 170) ] <- 2
esa.lc[ (esa.lc > 2) ]  <- 3 
# create the LC mask
lc.mask <- esa.lc
# forest mask
lc.for.mask <- lc.mask
lc.for.mask[lc.mask!=2] <- NA
lc.for.mask[lc.mask==2] <- 1


#only for tc accected area 
load("../TC_tack_mask/tc.occ.1999to2018.rda")
#only accounted the tc affected area 
lc.for.mask[arr.for.plot<=0.1] <- NA 

# variable name: esa.lc
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30
# set water(210) to NA
#generate the color palette 
lc.leg   <- read.csv("ESACCI-LC-Legend.csv",sep=";")
#lc.col    <- rgb(lc.leg$R/255, lc.leg$G/255, lc.leg$B/255)
#lc.breaks <- c(0, lc.leg$NB_LAB+0.1)

source("read_ncdf4.R")
nx=6722
ny=6722 

eff.acc.pos <- array( 0, dim=c( nx,ny) )  
eff.acc.neg <- array( 0, dim=c( nx,ny) )  


#nav.lat <- array(NA,dim=c(nx,ny)) 
#nav.lon <- array(NA,dim=c(nx,ny)) 
    
lai.1 <- fun_read_nc( arg1="/lfs/home/ychen/LAI_STUDY_EAsia/LAI_DATA/c_gls_LAI_202003200000_PROBAV_V2_EAST_ASIA.nc", var_st=2 )
#for (ix in 1:nx) {
#      for (iy in 1:ny) {
#         nav.lat[ix,iy] <- lai.1$lat[iy]
#         nav.lon[ix,iy] <- lai.1$lon[ix] 
#      }
#}


#load the coastlines data by readOGR function from sp package
coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_10m_coastline/ne_10m_coastline.shp")
#coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_110m_coastline/ne_110m_coastline.shp")
# Read this shape file with the rgdal library. 
#library(rgdal)
#land_mask <- readOGR( 
#  dsn= paste0("/lfs/home/ychen/GIS/world_border/") , 
#  layer="TM_WORLD_BORDERS_SIMPL-0.3",
#  verbose=FALSE
#)
#load river
#river <- readOGR("/lfs/home/ychen/GIS/ne_10m_rivers_lake/ne_110m_reviers_lake.shp")

#set work rda file
wrk.rda <- c("EA_19990220to20181231_table.combw_10_p_80_3D_pos_050_off_000_cut_0.rda")
# convert the land cover array to raster object 
# load the tc track data 

track.data <- data.frame()

for (iyr in 1999:2018) {
    
   tmp.data <- read.csv(file=paste(wrk.dir,track.dir,"/hourly_data/wpb_",iyr,"_hourly.txt",sep=""))
   track.data <- rbind(tmp.data,track.data)

}  

#ld.go <- F
#if (ld.go ) {
# initial arrays
#yr.id <- seq(1999, 2019,1) - 1945 + 1
#count.1999.to.2019.2d <- array( 0, dim=c(6722,6722))
#load TC occurence
load("../TC_tack_mask/tc.occ.1999to2018.rda")
#sum up all arraies in ann.land.frq
# sum up TC occurence
#for ( it in yr.id ) {
# print(paste("working on file:",it,"for 2D") )
#      count.1999.to.2019.2d <-  count.1999.to.2019.2d + ann.land.frq[,,it]
#    }
#count annual occurance of TC
#tc.ave.occ.2d <- count.1999.to.2019.2d / as.numeric(length(yr.id))
#tc.ave.occ.2d[tc.ave.occ.2d<=0.1] <- 0

#}else{
#load("tc.ave.occ.2d.rda")
#}#ld go 




#creat an array to account the forest area over each 10 by 10 degree 
xygd <- c(10)
nx10d <- (xymn[2]-xymn[1])/xygd
ny10d <- (xymn[4]-xymn[3])/xygd 
for.acc.10d <- array( 0, dim=c( nx10d,ny10d) )  
for.navlon.10d <- array( 0, dim=c( nx10d,ny10d) )  
for.navlat.10d <- array( 0, dim=c( nx10d,ny10d) )  

#array for accounting for positive or negtive effect size over 10 by 10 degree
eff.acc.10d.pos <- array( 0, dim=c( nx10d,ny10d) )  
eff.acc.10d.neg <- array( 0, dim=c( nx10d,ny10d) )  

# find forest pixel for each grid for 10 by 10 degree

for (ix in nx10d:1) {
    for (iy in ny10d:1) {
       # find pixel for the seach area 
        cn.lat <-  xymn[4]- (xygd/2) - ((ix-1)*xygd)
        cn.lon <-  xymn[2]- (xygd/2) - ((iy-1)*xygd)
        dl.lat <- xygd/2
        dl.lon <- xygd/2
        id.x.10d <- which ( (lai.1$lat > cn.lat-dl.lat) & (lai.1$lat <= cn.lat+dl.lat) )
        id.y.10d <- which ( (lai.1$lon > cn.lon-dl.lon) & (lai.1$lon <= cn.lon+dl.lon) )
        #count total pixels 
        for.navlon.10d[ix,iy] <- cn.lon 
        for.navlat.10d[ix,iy] <- cn.lat 
        for.acc.10d[ix,iy] <-  sum(lc.for.mask[id.x.10d,id.y.10d],na.rm=T)
    } # 
} #  
for.acc.10d <- t(for.acc.10d[nx10d:1,ny10d:1])


# load the effect size data
library("fields")
# load data table 
#load("./Rda//WP_20010110to20181231_table.combw_3_p_350.rda")

load(paste(table.dir,"/",wrk.rda, sep=""))

#QC1 mean difference  on LAI
#QC2 measurement uncertainty on LAI  
qc1.set <- 1.0
qc2.set <- 0.25
area.set <- 0

fun_table.qc <- function(wrk.table, qc1.set=0.5, qc2.set=0.2, area.set=1000) {
 
      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      #QC1 mean difference  on LAI
      #QC2 measurement uncertainty on LAI  
      #qc1.set <- 0.5
      #qc2.set <- 0.2
      wrk.table <- subset(wrk.table, (abs(wrk.table$qc1.score)<=qc1.set &  abs(wrk.table$qc2.score)>=qc2.set ))
      return(wrk.table)
      } 

#table quality check
table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set) 


# Generate a figure print plot with x-axis: months; y-axis: latitude
# define array 

nmon = 12
nlat = 20
fig.arr <- array(NA, dim=c(nmon,nlat))

mon.txt <- c("2.5","2.75","3.0","3.25","3.5","3.75","4.0","4.25","4.5","4.75","5.0","5.25")
lat.txt <- c(seq(10,50,2))

#create day column
table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)

ini.lat = 10.
del.lat = 2. 

for (ix in 1:12) {
    for (iy in 1:20) {

      lat1 <- ini.lat + del.lat * as.numeric(iy-1) 
      lat2 <- ini.lat + del.lat * as.numeric(iy)
      
      mon1 <- as.integer(ix)
      # subset the table to the condiction
      tmp.table <-  subset(table.comb, (table.comb$eve.lat.med.aff>lat1 ) & (table.comb$eve.lat.med.aff <=lat2) 
                                      & (table.comb$date.mon ==mon1 )   ) 
 
      #calculate the effect size in the table 
      tmp.eff.size <- mean(tmp.table$eff.size.for, na.rm=T)
 
      fig.arr[ix,iy] <- tmp.eff.size
                    
      print(paste("lat:",lat1,"to",lat2,"; ", "month:",mon1))
      print(tmp.table$eve.lat.med.aff)
      print(tmp.table$date.mon)
   }
}

#library(maptools)
library("raster")
library("rgdal")
# start to plot the figure print
#library("fields")

#my.col    <- colorRampPalette(c("darkblue","blue","cyan","lightgray","yellow","pink","red","brown","brown"))(25)
#eff.col    <- colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(16)
#eff.col    <- colorRampPalette(c("forestgreen","green","lightgreen","white","gray","orange","brown"))(16)
eff.col    <- rev(colorRampPalette(c("forestgreen","green","lightgreen","white","pink","red","brown"))(20) ) 
#eff.col    <- colorRampPalette(c("orange4","orange2","white","forestgreen","darkgreen"))(20)
eff.breaks <- c( seq(-1., 1., length.out = 19))
eff.breaks <- round(eff.breaks,2)

#assign the color base on effect size 
table.comb$col <-  eff.col[as.numeric(cut(table.comb$eff.size.for ,breaks = eff.breaks))]


#assign the maximum radius to the event
  for (ieve in 1:length(table.comb$date.time) ) {
      #subset track.data to event table 
      eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.comb$date.yr[ieve]) &
                                          (track.data$CY == table.comb$eve.id[ieve]) ))
      table.comb$RAD[ieve] <- max(eve.table$RAD)
  }  


#convert array  to raster obj for ploting 
r_p <-raster( t(fig.arr[,nlat:1]) ,
               xmn=2, xmx=13,
               ymn=10, ymx=50, 
               crs=CRS("+proj=longlat +datum=WGS84") )

#plot(r_p, xlab="month", ylab="Latitude",cex.lab=1.5,
#    xlim=c(1,24),ylim=c(8,38),axes=FALSE,
#    legend.args=list(text="Eff. Size", side=3, font=2, line=0.2, cex=1.0),legend.width=2,
#    col=my.col,breaks=my.breaks )

#x_pos <- as.numeric(table.comb$bef.for.lai.aff)
x_pos <- (table.comb$date.mon+ as.integer(table.comb$date.day)/30 )
y_pos <- as.numeric(table.comb$eve.lat.med.aff)
#eff_cex <- 2.0* sqrt(table.comb$eve.for.pix.aff/3.1416)/100 
#eff_cex <- 2.5 
eff_cex <-  table.comb$RAD

# calculate the accumulative effect size for the grid of 10 by 10 degree 
for (ieve in 1:length(table.comb$RAD) ) {
     #allocate array
     tmp.arr <- array(0,dim=c(nx,ny))
     #for positive event
     if (table.comb$eff.size.for[ieve] > 0.) {
        #find pixel for the events
        x.id <- which(  (lai.1$lon <= (table.comb$eve.lon.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lon >= (table.comb$eve.lon.med.aff[ieve]-eff_cex[ieve]))  ) 
 
        y.id <- which(  (lai.1$lat <= (table.comb$eve.lat.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lat >= (table.comb$eve.lat.med.aff[ieve]-eff_cex[ieve]))  ) 
        tmp.arr[c(x.id),c(y.id)] <- 1    
        eff.acc.pos <- tmp.arr + eff.acc.pos
     }#end if
     #for negtive event 
     if (table.comb$eff.size.for[ieve] < 0.) {
        x.id <- which(  (lai.1$lon <= (table.comb$eve.lon.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lon >= (table.comb$eve.lon.med.aff[ieve]-eff_cex[ieve]))  ) 
 
        y.id <- which(  (lai.1$lat <= (table.comb$eve.lat.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lat >= (table.comb$eve.lat.med.aff[ieve]-eff_cex[ieve]))  ) 
        tmp.arr[c(x.id),c(y.id)] <- 1    
        eff.acc.neg <- tmp.arr + eff.acc.neg
 
     }#end if

}#end for


# aggregate to 10 degree by 10 degree
for (ix in 1:nx10d) {
   for (iy in 1:ny10d) {
    #set spacinf for 10 degree
    x.del =xygd/2
    y.del =xygd/2
    #find xid
    x.id <- which (  (lai.1$lon >= for.navlon.10d[ix,iy]- x.del) &
                     (lai.1$lon <  for.navlon.10d[ix,iy]+ x.del) )  
    #find yid    
    y.id <- which (  (lai.1$lat >= for.navlat.10d[ix,iy]- y.del) &
                     (lai.1$lat <  for.navlat.10d[ix,iy]+ y.del) )  
    #aggregate by sum 
    eff.acc.10d.neg[ix,iy] <- sum(eff.acc.neg[c(x.id),c(y.id)],na.rm=T) 
    eff.acc.10d.pos[ix,iy] <- sum(eff.acc.pos[c(x.id),c(y.id)],na.rm=T) 
   }#end for
}#end for 

#eff.acc.10d.neg <- t(eff.acc.10d.neg[nx10d:1,ny10d:1])
#eff.acc.10d.pos <- t(eff.acc.10d.pos[nx10d:1,ny10d:1])




#
#points(x=x_pos , y=y_pos, col=table.comb$col, bg="NA",  cex=my_cex, pch=21)
#mtext(text=c(paste(mon.txt[1:12])), side=1, line=-1.5, at=seq(1,24,2)+0.5, las=2, cex=0.8)
#mtext(text=c(paste(lat.txt[1:15])), side=2, line=-2, at=seq(9,37,2), las=2, cex=0.8)
#box()

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    plot(x=c(2.,12.9), y=c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main='')
    axis(2, round(ticks,digits=1), las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

if (go.pdf) { pdf(file="acc.plot_b.pdf",width=8, height=8) }
#layoutsetting for combine plot 
par(oma = c(2, 0.5, 0.1, 0.1))
layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3, 1, 1), c(1, 3), TRUE)
op = par(no.readonly = TRUE)

out_col <- colorRampPalette(c("lightgray","gray","black"))(length(x_pos))
 
#plot 2 scater plot of the effect size
par(mar = c(4, 4, 0, 0))
   plot(x=x_pos , y=y_pos, xlab="Month", ylab="Latitude", cex.lab=1.5, cex.axis=1.5,
        col="black", bg=table.comb$col, xlim=c(2.0,12.9), ylim=c(5,50), cex=eff_cex, pch=21,lwd=1.5)
   grid(lwd=2)




#plot 4  averaged LAI pre-condition  vs  the effect size
par(mar = c(0, 4, 3, 0))
   #plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$date.mon)),
   #     col="lightgray", outline=FALSE,ylim=c(-0.8,0.8),main="",
   #     xlab="",ylab="Effect size",  xaxt = "n" )
    
   
   plot(y=table.comb$eff.size.for, x=(table.comb$date.mon+ as.integer(table.comb$date.day)/30 ),
         col="white",ylim=c(-1.,1.),xlim=c(2.0,12.9),
         xlab="", ylab="Effect size",main="", xaxt = "n" )
 

   allrun.table <- NA
   #get file-list 
   in.files <- list.files(path=table.dir, pattern="EA_*")
   in.files[1] <- in.files[7]
     
   for (id in 1:length(in.files)  ) { 
   #for (id in 1:1) { 
   load(paste(table.dir,in.files[id],sep=""))
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   # subset datataset (this is super critical!! we develop this QC flags)
   table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set,area.set=area.set)
   #filter outliner 
   table.comb <- subset(table.comb, (table.comb$eff.size.for) <= 1.5 )
   sm.ln <- smooth.spline(y=table.comb$eff.size.for, x=(table.comb$date.mon+ as.integer(table.comb$date.day)/30), df=10)
   lines(x=sm.ln$x, y=sm.ln$y, col = "darkgray",lwd=2,lty="dashed")
   #text(x=sm.ln$x[1], y=sm.ln$y[1],id)
   # combine all data to a table 
   allrun.table <- rbind(table.comb,allrun.table)
 }
   
   #add  all dataset avaerage information 
   allrun.table <- allrun.table[complete.cases(allrun.table), ]  
 
   sm.ln <- smooth.spline(y=allrun.table$eff.size.for,x=allrun.table$date.mon, df=10) 
   lines(x=sm.ln$x, y=sm.ln$y, col = "black",lwd=2)
   #grid(col="gray")
   abline(h=0,col="black",lty="dashed")
   #plot 1 zonal average of the effect size
   par(mar = c(4, 0, 0, 0))
   #plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$eve.lat.med.aff)),
   #     horizontal=TRUE, col="orange", ylim=c(-0.5,0.5), outline=FALSE,
   #     xlab="", ylab="Effect size",main="", yaxt = "n")
   
   plot(x=table.comb$eff.size.for, y=table.comb$eve.lat.med.aff,ylim=c(0,50), 
        col="white",xlim=c(-1.,1.) ,ylab="", xlab="Effect size",main="", yaxt = "n" )
  
   #scatter.smooth(x=table.comb$eff.size.for, y=table.comb$eve.lat.med.aff,
   #               col="gray", span=0.5, degree=1, xlab="", ylab="Effect size", 
   #               lpars=list(lwd=2,col="black"), xlim=c(-0.3,0.3))
   #grid(col="gray")
   allrun.table <- NA
   for (id in  1:length(in.files)  ) {
   load(paste(table.dir,"/",in.files[id],sep=""))
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   #call table quality check
   table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set) 
   #filter outliner 
   table.comb <- subset(table.comb, (table.comb$eff.size.for) <= 1.5 )
   sm.ln <- smooth.spline(y=table.comb$eff.size.for,x=table.comb$eve.lat.med.aff, df=10)
   lines(x=sm.ln$y, y=sm.ln$x, col = "darkgray",lwd=2, lty="dashed")
   #points(y=sm.ln$y, x=sm.ln$x, col="gray")
   #text(x=sm.ln$y[1], y=sm.ln$x[1],id)

   # combine all data to a table 
   allrun.table <- rbind(table.comb,allrun.table)
   }

   #add aveage of all dataset 
   allrun.table <- allrun.table[complete.cases(allrun.table), ]  
   sm.ln <- smooth.spline(y=allrun.table$eff.size.for,x=allrun.table$eve.lat.med.aff, df=15)
   #points(y=allrun.table$eve.lat.med.aff, x=allrun.table$eff.size.for,col="gray") 
   lines(x=sm.ln$y, y=sm.ln$x, col = "black",lwd=2)
  
   abline(v=0,col="black",lty="dashed")
   
 
# plot 3 legend on the top-right 
par(mar = c(0.5, 5, 3.5, 4), xaxt = "n")
    color.bar(lut=eff.col,min=-1.,title="Effect size scale")
#close devic()
if(go.pdf){dev.off()}



ld.go <- TRUE
if (ld.go) {
#open new device
go.pdf=T
if(go.pdf) {pdf(file="eff_plot_sub_a.pdf",width=8, height=8) } else{dev.new() }
#reset plot parameters
#par(mar = c(0, 0, 0, 0), mgp=c(3,1,0),xpd=FALSE)
par(mar = c(4, 4, 1, 2) , mgp=c(3,1,0))
#set graph parameter
#layout(matrix(c(1),1 , 1, byrow = TRUE), widths=c(1), heights=c(1), TRUE)
par(oma = c(1, 1, 1, 1))
#par(mar = c(4, 4.5, 0, 0), xaxs="i", yaxs="i")
# assign the table for plot 
load(paste(table.dir,"/",wrk.rda,sep=""))
table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
#call table quality check
table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set) 
#assign the color base on effect size 
table.comb$col <-  eff.col[as.numeric(cut(table.comb$eff.size.for ,breaks = eff.breaks))]
load("/lfs/home/ychen/LAI_STUDY_EAsia/LANDCOVER_DATA/2015.esa.landcover.east.asia.rda")
#
#asign lengend classes and texts
#x.lab   <- c(0,10,20,30,40,50,60,70,80,90,
#             100,110,120,130,140,150,160,170,180,190,200,210)
#txt.lab <- c("No-Obs","Agr-Rf","Agr-Ir","AF-M1","AF-M2","For-BE","For-BD","For-NE","For-ND","For-Mix",
#             "Mos-Tr","Mos-Hb","Shrub","Grass","Mosses","Spas-V","Fld-T1","Fld-T2","Fld-Hb","Urban","Bare","Water")
 x.lab   <- c(0, 40,90,120,180,190,200,210,230,250)
txt.lab  <- c("No-Obs","Agr-all","For-all","Grass","Bare","Urban","Water","Sea")
#lc.col  <- c("gray","yellow","forestgreen","darkolivegreen","orange","brown","orange","white")
lc.col  <- c("NA","NA","lightgray","NA","NA","NA","NA","NA")
lc.breaks <- c(0, x.lab+0.5)
  #replace x_pos as longitude
   x_pos <- as.numeric(table.comb$eve.lon.med.aff)   
  # set subplot margin
  par( usr=c(100,150,0,60)) #no buffer in x-y axis 
# par( xaxs="i", yaxs="i")
  plot(coastlines,col="white",lwd=0,ylim=c(0,60),xlim=c(100,150), ylab="", xlab="")
 #plot  raster land cover layer
 #plot legend
  #add coastaline
#  par( usr=c(100,150,0,60)) #no buffer in x-y axis 
  par( usr=c(100,150,0,60)) #no buffer in x-y axis 

  plot(lc.raster, col=lc.col, breaks=lc.breaks, legend=F, ylim=c(0.,60.),xlim=c(100,150),,add=T)
  par( usr=c(100,150,0,60), xpd=FALSE)  
  plot(coastlines, lwd=0.8,xlim=c(100,150), col="black", add=T) 
  grid(lwd=1,col="black")
  axis(side = 2, at = seq(0,60,length.out=7),
       labels = c("0","10","20","30","40","50","60"), tck = -0.02)
  mtext("Latitude" , side=2, line=2.5, cex=1.5)

  axis(side = 1, at = seq(90,150,length.out=7),
       labels = c("","100","110","120","130","140","150"), tck = -0.02)
  mtext("Longitude" , side=1, line=2.5, cex=1.5)
  box() 
   # add color bar
   colorbar.plot(x=135.,y=3,adj.y=.5,adj.x=0.0, strip=seq(-1,1,0.1), strip.width=.02, strip.length=0.02*10,
                horizontal=T,col=eff.col)
   text(x=140,y=1.0,srt=0,"Effect Size",cex=1.0)
   text(x=140,y=4.2,srt=0,"-1.0    -0.5    0     0.5     1.0 ",cex=.8) 
   #add land sea mask
   #add event-base TC tracks
   t.blue=rgb(red=0.1, green=.2, blue=.8, alpha=.5,  maxColorValue = 1)
   t.red=rgb(red=0.8, green=.2, blue=.1, alpha=.5,  maxColorValue = 1)
   t.green=rgb(red=0.1, green=.8, blue=.1, alpha=.5,  maxColorValue = 1)

   for (ieve in 1:length(table.comb$date.time) ) {
      #subset track.data to event table 
      eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.comb$date.yr[ieve]) &
                                          (track.data$CY == table.comb$eve.id[ieve]) ))
      eve.table  <- subset (eve.table, eve.table$VMAX >=35) 
      eve.stp <- length(eve.table$VMAX)
      eve.rad <- max(eve.table$RAD)
      print(paste("TC max RAD:",eve.rad,", ","TC_stps:",eve.stp,sep=""))
      # plot the tracks
      lines( x=eve.table$LonEW, y=eve.table$LatNS, lty="solid",
             col=table.comb$col[ieve],
             lwd=ifelse( abs(table.comb$eff.size.for[ieve])>0.4  ,0.25, 0.25))
  } 
#add coastline
  par( usr=c(100,150,0,60), xpd=FALSE)  
  plot(coastlines, lwd=0.8,xlim=c(100,150), col="black", add=T) 
 
  # add forest area that affected by TCs

    eff.acc.10d.neg[eff.acc.10d.neg==0]<- 1
    eff.acc.10d.pos[eff.acc.10d.pos==0]<- 1
    for (ix in 1:nx10d) {
      for (iy in 1:ny10d) {
          # add pie chart 
            if (for.acc.10d[ix,iy] >= 1000) {
                add.pie(z=c(eff.acc.10d.pos[ix,iy],eff.acc.10d.neg[ix,iy]),
                        x=for.navlon.10d[ix,iy], y=for.navlat.10d[ix,iy], 
                        col=c(t.green,t.red), radius=2.5, labels=NA)
                area.txt <- round(for.acc.10d[ix,iy]/10000, digits=1) 
                text(x=for.navlon.10d[ix,iy]+2.5, y=for.navlat.10d[ix,iy]-4.0,labels=paste(area.txt,"M",sep=""))
             }  
      } 
    }

    #add coastlines

#  plot(coastlines,lwd=1.2, add=T) 
 
 
   par(bty = "n")
   #add legend
   #plot(lc.raster, col=lc.col, breaks=lc.breaks, legend.only=T, horizontal=F,
   #     add=T, smallplot=c(0.93,0.96,0.2,0.9), legend.width=0.25, legend.shrick=0.5,
   #     axis.arg=list(at= x.lab,
   #                   labels= txt.lab ,
   #                    las=2, cex.axis=0.8, mgp=c(5,0.5,0)),  
   #     legend.arg= list(text="Land cover type",las=3, side=2, font=2, line=0.2, cex=0.8) )
#reset plot parameters
par(par.def)
if(go.pdf) {dev.off()}
}
#par(oma = c(0.5, 0.5, 0.5, 0.5))
