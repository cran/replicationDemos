# Rasmus Benestad
# Oppegård 07.09.2012
# ----------------------------------------
# README:
# This script runs on R which is freely available from
# http://cran.r-project.org
# R runs on many different platforms.
# This script requires Internet-access and data serves in service.
#
# To run this script start R and type the following lines:
# > download.file('http://www.realclimate.org/images//HSS2012.txt','hss2012.R')
# > source('hss2012.R')
#
# Purpose:
# Replication of: 'The phase relation between atmospheric carbon dioxide
#                  and global temperature'
# Ole Humlum, Kjell Stordahl, Jan-Erik Solheim
# http://dx.doi.org/10.1016/j.gloplacha.2012.08.008,
# http://www.sciencedirect.com/science/article/pii/S0921818112001658

# Cut-away parts of the data
# Time-scale issue: emphasises short-term variations - CO2 expected to
# matter on longer time scales.
# Focuses on older HadCRUT3, although other data are also cited. 
# log CO2 rather than CO2
# d.o.f.s? significance testing
# Expect a convoluted response, due to oceans. Naive simple.
# Natural/internal variations on inter-annual to decadal time scales
# Filtering/differenting - appropriate? Jan-Jan, Feb-Feb, etc.
# Smoothing with windows longer than 12 months - longer time scales.
# El Ninos: http://www.pmel.noaa.gov/tao/elnino/el-nino-story.html
# 

diff12 <- function(x,wfl=NULL) {
  yymm <- attr(x,'yymm')
  yy <- trunc(yymm)
  mm <- round(12*(yymm-yy) + 0.5)
  if (!is.null(wfl)) x <- filter(x,rep(1,wfl)/wfl)
  years <- as.numeric(row.names(table(yy)))
  fullyr <- as.numeric(table(yy))
  #print(years)
  years <- years[is.element(fullyr,12)]
  complete <- is.element(yy,years)
  mm <- mm[complete]
  yy <- yy[complete]
  x <- x[complete]
  #print(table(mm))
  #print(years)
  #print(fullyr)
  n <- length(years)-1
  diff12 <- matrix(rep(NA,n*12),12,n)
  for (m in 1:12) {
    ii <- is.element(mm, m) 
    #print(c(m,length(ii),sum(ii),length(diff12[m,]),length(diff(x[ii]))))
    diff12[m,] <- diff(x[ii])
  }
  diff12 <- c(diff12)
  attr(diff12,'yymm') <- yy[-(1:12)] + (mm[-(1:12)]-0.5)/12
  invisible(diff12)
}



Humlum.et.al.2012 <- function(wfl=12,forcing=FALSE,HadCRUT4=FALSE,HadSST3=FALSE) {
  
  ftp.co2="ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/co2_mm_mlo.txt"
  ftp.hc3="http://www.cru.uea.ac.uk/cru/data/temperature/hadcrut3gl.txt"
  ftp.hs2="http://www.cru.uea.ac.uk/cru/data/temperature/hadsst2gl.txt"
  ftp.CO2="ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/co2_mm_gl.txt"
  ftp.hs3="http://www.metoffice.gov.uk/hadobs/hadsst3/data/TS_all_realisations.zip"
  ftp.hc4="http://www.metoffice.gov.uk/hadobs/hadcrut4/data/time_series/hadcrut4_monthly_ns_avg.txt"

  print("This function requires Internet access to data sources.")
  print("Downloading can take some time.")
  if ((HadSST3) & (!file.exists("TS_all_realisations.txt"))) {
    print("Using HadSST# rather than HadSST2")
    download.file(ftp.hs3,"HadSST3.zip")
    unzip("HadSST3.zip")
    SST <- read.table("TS_all_realisations.txt",fill=TRUE)
    sst <- SST$V3
    attr(sst,'yymm') <- SST$V1 + (SST$V2 - 0.5)/12
  } else {
    SST <- read.table(ftp.hs2,fill=TRUE)
    sst <- c(t(as.matrix(SST)[is.finite(SST$V14),2:13]))
    attr(sst,'yymm') <- sort(rep(SST$V1[is.finite(SST$V14)],12)) +
                     (rep(1:12,sum(is.finite(SST$V14)))-0.5)/12

  }
  if (HadCRUT4) {
    HC <- read.table(ftp.hc4)
    t2m <- HC$V2
    attr(t2m,'yymm') <- as.numeric(substr(HC$V1,1,4)) +
                       (as.numeric(substr(HC$V1,6,7))-0.5)/12
  } else { 
    HC <- read.table(ftp.hc3,fill=TRUE)
  t2m <- c(t(as.matrix(HC)[is.finite(HC$V14),2:13]))
  attr(t2m,'yymm') <- sort(rep(HC$V1[is.finite(SST$V14)],12)) +
                            (rep(1:12,sum(is.finite(HC$V14)))-0.5)/12
  }
  
  CO2 <- read.table(ftp.co2)
  CO.2 <- read.table(ftp.CO2)

  co2 <- CO2$V5
  co.2 <- CO.2$V4
  if (forcing) {
    co2 <- log(co2)
    co.2 <- log(co.2)
  }
  attr(co2,'yymm') <- CO2$V1 + (CO2$V2 - 0.5)/12
  attr(co.2,'yymm') <- CO.2$V1 + (CO.2$V2 - 0.5)/12

# Remove entries not yet filled in - set to NA.
  remove.sst <- (sst == 0) & (attr(sst,'yymm') >2000)
  sst[remove.sst] <- NA
  remove.t2m <- (t2m == 0) & (attr(t2m,'yymm') >2000)
  t2m[remove.t2m] <- NA

  sco2 <- 70; bco2 <- -0.1
  yg <- (co2 - 300)/sco2 + bco2
  yG <- (co.2 - 300)/sco2 + bco2

  dev.new(width=8.75,height=6)
  par(mar=c(5.1, 4.1, 4.1, 4.1),xaxt="n",yaxt="n")
  plot(attr(t2m,'yymm'),t2m,type="l",lwd=1,col="red",xlab="",
       ylab="Original data",
       xlim=c(1980,2012),ylim=c(-0.4,1.3),main="HSS2012 - fig1",
       sub=paste("filtering:",wfl))
  par(las=2,cex.axis=0.8,xaxt="s",yaxt="s")
  rect(1982.5,-0.6,1983.5,1.6,density=30,col="grey80")
  rect(1986.5,-0.6,1987.5,1.6,density=30,col="grey80")
  rect(1994.5,-0.6,1995.5,1.6,density=30,col="grey80")
  rect(1997.5,-0.6,1998.5,1.6,density=30,col="grey80")
  axis(1,at=seq(1980,2012,by=2))
  axis(2,at=seq(-0.4,1.2,by=0.2))
  co2tcks <- seq(280,400,length=13)
  axis(4,at=(co2tcks - 300)/sco2 + bco2,labels=co2tcks)
  lines(attr(sst,'yymm'),sst,col="blue",lty=2)
  lines(attr(co.2,'yymm'),yG,col="darkgreen",lwd=2)
  lines(attr(co2,'yymm'),yg,col="green",lty=2)

  text(1983,1.3,'El Nino',col="grey",font=2,cex=0.8)
  text(1987,1.3,'El Nino',col="grey",font=2,cex=0.8)
  text(1998,1.3,'El Nino',col="grey",font=2,cex=0.8)

  legend(2002,-0.1,c(expression(paste("glob.",CO[2])),"Mauna Loa",
      "T(2m)","SST"),col=c("darkgreen","green","red","blue"),
      lty=c(1,2,1,2),lwd=c(2,1,1,1),bty="n",cex=0.7)   
                     
  
  dev2bitmap(file=paste("hss2012fig1-",wfl,"-",forcing,HadCRUT4,
               HadSST3,".png",sep=""),res=150)
  t2m.d <- diff12(t2m,wfl=wfl)
  co2.d <- diff12(co2,wfl=wfl)
  short <- attr(co2,'yymm')>1980
  co2s <- co2[short]
  attr(co2s,'yymm') <- attr(co2,'yymm')[short]
  co2s.d <- diff12(co2s,wfl=wfl)
  co.2.d <- diff12(co.2,wfl=wfl)

  sd12 <- 0.3; cd12 <- 0

  dev.new(width=8.75,height=6)
  par(mar=c(5.1, 4.1, 4.1, 4.1),xaxt="n",yaxt="n")
  plot(attr(t2m.d,'yymm'),t2m.d,type="l",lwd=2,col="red",
     xlim=c(1982,2012),ylim=c(-0.5,1),main="HSS2012 - fig2",
     sub=paste("filtering:",wfl),xlab="",
     ylab="Diff12- results")
  rect(1982.5+1,-0.6,1983.5+1,1.6,density=30,col="grey80")
  rect(1986.5+1,-0.6,1987.5+1,1.6,density=30,col="grey80")
  rect(1994.5+1,-0.6,1995.5+1,1.6,density=30,col="grey80")
  rect(1997.5+1,-0.6,1998.5+1,1.6,density=30,col="grey80")
  lines(attr(co.2.d,'yymm'),(co.2.d-cd12)*sd12,col="darkgreen",lwd=2)
  lines(attr(co2.d,'yymm'),(co2.d-cd12)*sd12,col="green",lty=2)
  par(las=1,cex.axis=0.8,xaxt="s",yaxt="s")
  axis(1,at=seq(1980,2012,by=2))
  axis(2,at=seq(-0.4,1.2,by=0.2))
  dco2tcks <- seq(-1,3,length=5)
  axis(4,at=(dco2tcks - cd12)*sd12 + bco2,labels=dco2tcks)

  text(1983+1,1.0,'El Nino',col="grey",font=2,cex=0.8)
  text(1987+1,1.3,'El Nino',col="grey",font=2,cex=0.8)
  text(1998+1,1.0,'El Nino',col="grey",font=2,cex=0.8)
    dev2bitmap(file=paste("hss2012fig2b-",wfl,"-",forcing,HadCRUT4,
               HadSST3,".png",sep=""),res=150)

  dev.new()
  i1 <- is.element(attr(t2m.d,'yymm'),attr(co.2.d,'yymm')) &
      attr(t2m.d,'yymm') >= 1982
  i2 <- is.element(attr(co.2.d,'yymm'),attr(t2m.d,'yymm')) &
      attr(co.2.d,'yymm') >= 1982
  lr1 <- ccf(co.2.d[i2],t2m.d[i1],na.action = na.pass,plot=FALSE) 
  i1 <- is.element(attr(t2m.d,'yymm'),attr(co2.d,'yymm'))
  i2 <- is.element(attr(co2.d,'yymm'),attr(t2m.d,'yymm'))
  lr2 <- ccf(co2.d[i2],t2m.d[i1],na.action = na.pass,plot=FALSE)
  i1 <- is.element(attr(t2m.d,'yymm'),attr(co2s.d,'yymm')) &
      attr(t2m.d,'yymm') >= 1982
  i2 <- is.element(attr(co2s.d,'yymm'),attr(t2m.d,'yymm')) &
      attr(co2s.d,'yymm') >= 1982
  lr3 <- ccf(co2s.d[i2],t2m.d[i1],na.action = na.pass,plot=FALSE)
  plot(lr1$lag,lr1$acf,type="l",lwd=2,
       xlab="lag",ylab="correlation",
     xlim=c(-12,24),ylim=c(-0.4,0.8),main="HSS2012 - fig4b")
  grid()
  lines(lr3$lag,lr3$acf,col="grey20",lwd=2,lty=2)
  lines(lr2$lag,lr2$acf,col="grey",lwd=2)

  lines(c(-12,24),rep(0.4,2),lty=2)
  lines(rep(10,2),c(-0.4,0.6),lty=2)

  legend(-12,0.8,c("HSS2012","Keeling 1980-","Keeling 1958-"),
         col=c("black","grey20","grey"),lty=c(1,2,1),lwd=2,bty="n")
  dev2bitmap(file=paste("hss2012fig4b-",wfl,"-",forcing,HadCRUT4,
               HadSST3,".png",sep=""),res=150)
}

diffdemo <- function(x=0.7*cos(seq(0,10*pi,length=1000))+0.4*rnorm(1000),
                     y=0.9*cos(seq(0,10*pi,length=1000))+0.3*rnorm(1000)) {

  dev.new(width=6,height=8)
  par(mfcol=c(2,1),xaxt="n",yaxt="n",bty="n")
  plot(x+0.5,type="l",main="diffdemo()",ylab="x & y",xlab="",
       sub="Highly correlated slow variations + noise",
       ylim=c(-4,3))
  lines(y+0.5,col="red")

  dx <- diff(x); dy <-diff(y)
  lines(dx - 2.5)
  lines(dy -2.5,col="red")
  text(1000,0.5,pos=4,srt=90,"Original",col="grey",cex=0.7)
  text(1000,-3,pos=4,srt=90,"Differentiated",col="grey",cex=0.7)

  par(xaxt="s",yaxt="s",bty="n")
  ccf(dx,dy,lwd=3,ylim=c(-1,1),
      sub="lag correlation of differentiated series")
  lines(rep(0,2),c(-1,1),lty=2,col="grey")
  grid()
  dev2bitmap(file="diffdemo.png",res=150)
}

diff12demo <- function(x=0.5*cos(seq(0,10*pi,length=1200))+rnorm(1200),
                      y=0.7*cos(seq(0,10*pi,length=1200))+rnorm(1200),
                       wfl=12) {
  attr(x,'yymm') <- sort(rep(1:100,12)) + (rep(1:12,100)-0.5)/12
  attr(y,'yymm') <- sort(rep(1:100,12)) + (rep(1:12,100)-0.5)/12

  dev.new(width=6,height=8)
  par(mfcol=c(2,1),xaxt="n",yaxt="n",bty="n")
  plot(x+0.5,type="l",main="diff12demo()",ylab="x & y",xlab="",
       sub="Highly correlated slow variations + noise",
       ylim=c(-4,3))
  lines(y+0.5,col="red")

  dx <- diff12(x,wfl=wfl)
  dy <- diff12(y,wfl=wfl)
  
  lines(dx - 2.5)
  lines(dy -2.5,col="red")
  text(1000,0.5,pos=4,srt=90,"Original",col="grey",cex=0.7)
  text(1000,-3,pos=4,srt=90,"Differentiated",col="grey",cex=0.7)

  par(xaxt="s",yaxt="s",bty="n")
  ccf(dx,dy,lwd=3,ylim=c(-1,1),na.action = na.pass,
      sub="lag correlation of differentiated series")
  lines(rep(0,2),c(-1,1),lty=2,col="grey")
  grid()
  dev2bitmap(file="diff12demo.png",res=150)
}

#diffdemo()
#diff12demo()
#hss2012()

