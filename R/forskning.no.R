forskning.no <- function(uah=TRUE,rss=TRUE,giss=TRUE,ncdc=TRUE,hadcrut3=TRUE,
                         lwd=1,ylim=c(-0.4,1),xlim=c(1995,2012),
                         base.period=1981:2010,type="s") {
 
# Replication of Humlum et al in forskning.no (Onsdag 14. mars 2012)
# RSS MSU(black), GISS(red), NCDC(violet) and HadCRUT3(green)
# http://www.forskning.no/artikler/2012/mars/316178

  anomaly <- function(X,yr,mo,base.period) {
    # Estimate anomalies with respect to given reference period.
    # each calendar month is done seperately
    anom <- X + NA
    iiref <- is.element(yr,base.period)
    for (im in 1:12) {
      iim <- is.element(mo,im)
      anom[iim] <- X[iim] - mean(X[iim & iiref],na.rm=TRUE)
    }
    anom
  }

  annual <- function(X,date) {
    # Estmimate annual mean from monthle data:
    yr <- trunc(date)
    yrs <- as.numeric(rownames(table(yr)))
    #print(yrs)
    anm <- yrs + NA
    for (i in 1:length(anm))
      anm[i] <- mean(X[is.element(yr,yrs[i])])
    anm              
  }

  L <- function(yrs) {
    # This function extract one single year from monthly data
    l=as.numeric(rownames(table(trunc(yrs))))
    l
  }

strleg <- c("UAH","RSS","GISS","NCDC","HadCRUT3")
colleg <- c("blue","black","red","violet","darkgreen")
ileg <- c(uah,rss,giss,ncdc,hadcrut3)

uahmsuurl <- "http://vortex.nsstc.uah.edu/data/msu/t2lt/uahncdc.lt"
rssmsuurl <- "http://www.remss.com/data/msu/monthly_time_series/RSS_Monthly_MSU_AMSU_Channel_TLT_Anomalies_Land_and_Ocean_v03_3.txt"
ncdcurl <- "ftp://ftp.ncdc.noaa.gov/pub/data/anomalies/monthly.land_ocean.90S.90N.df_1901-2000mean.dat"
hadcrut3url <- "http://www.metoffice.gov.uk/hadobs/hadcrut3/diagnostics/global/nh+sh/monthly"
gissurl <- "http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt"

ref <- paste(min(base.period),max(base.period),sep=" - ")

dev.new()  
par(bty="n",las=2)
plot(x=xlim,y=ylim,ylim=ylim,xlim=xlim,type="n",
     main="Replication of Humlum et al in forsking.no (March 2012)",
     sub="forskning.no (replicationDemos)",
     xlab="time",ylab=paste("anomaly wrt",ref,"(deg C)"))  
grid()
par(las=0)
mtext("http://www.forskning.no/artikler/2012/mars/316178",side=4,col="grey")

legend(xlim[1],ylim[2],strleg[ileg],col=colleg[ileg],
       bty="n",lty=1,lwd=lwd,cex=0.75)
       
if (uah) {
  # Reference period 1981-2010
  print("UAH")
  test <- readLines(uahmsuurl)
  lwy <- grep("Year",test)
  nrows <- lwy[2]-2
  uahmsu <- read.table(uahmsuurl,header=TRUE,nrows=nrows)
  uahtemp <- anomaly(uahmsu$Globe,uahmsu$Year,uahmsu$Mo,base.period)
  lines(uahmsu$Year + (uahmsu$Mo -0.5)/12, uahtemp,lwd=lwd,
        col="blue",type=type)
}

if (rss) {
  # Brightness temperature anomalies are the difference between the
  # monthly brightness temperatures and the average value for that
  # month (found by averaging that month from 1979 through 1998).
  print("RSS")
  rsscn <- c("year","mon","S70.0toN82.5","S20.0toN20.0","N20.0toN82.5",
             "S70.0toS20.0","N60.0toN82.5","S70.0toS60","Cont.USA",
             "EqtoN82.5","S70.0toEq")
  rssmsu <- read.table(rssmsuurl,skip=3,col.names=rsscn)
  rsstemp <- anomaly(rssmsu$S70.0toN82.5,rssmsu$year,rssmsu$mon,base.period)
  lines(rssmsu$year + (rssmsu$mon - 0.5)/12,rsstemp,lwd=lwd,type=type)
}

if (giss) {
  # GLOBAL Land-Ocean Temperature Index in 0.01 degrees Celsius
  # base period: 1951-1980
  print("GISS")
  test <- readLines(gissurl)
  cyrs <- substr(test,1,4)
  vdat <- is.element(cyrs,as.character(1880:2050))
  writeLines(con="giss.test.txt",text=test[vdat])
  widths=c(rep(5,13),7,rep(5,6))
  gisscn <- c("year","Jan","Feb","Mar","Apr","May","Jun","Jul",
              "Aug","Sep","Oct","Nov","Dec","Jan.Dec","Dec.Nov",
              "DJF","MAM","JJA","SON","Year")
  GISS <- read.fwf("giss.test.txt",widths=widths,col.names=gisscn,
                   na.strings = "*****")
  gisstemp <- c(t(GISS[,2:13])) * 0.01
  gissdate <- sort(rep(GISS[,1],12)) + (rep(1:12,length(GISS[,1]))-0.5)/12
  gistemp <- anomaly(gisstemp,sort(rep(GISS[,1],12)),
                     rep(1:12,length(GISS[,1])),base.period)
  lines(gissdate,gistemp,lwd=lwd,col="red",type=type)
}

if (ncdc) {
  # Base period: 1901-2000.
  print("NCDC")
  NCDC <- read.table(ncdcurl,col.names=c("year","mon","global"))
  NCDC$global[NCDC$global < -99] <- NA
  ncdctemp <- anomaly(NCDC$global,NCDC$year,NCDC$mon,base.period)
  lines(NCDC$year + (NCDC$mon - 0.5)/12, ncdctemp,col="violet",lwd=lwd,
        type=type)
}

if (hadcrut3) {
  # Base period: 1961-1990
  print("HadCRUT3")
  HadCRUT3 <- read.table(hadcrut3url)
  crudate <- as.numeric(substr(HadCRUT3$V1,1,4)) +
            (as.numeric(substr(HadCRUT3$V1,6,7))-0.5)/12
  crutemp <- anomaly(HadCRUT3$V2,as.numeric(substr(HadCRUT3$V1,1,4)),
                     as.numeric(substr(HadCRUT3$V1,6,7)),base.period)
  lines(crudate,crutemp,col="darkgreen",type=type)
}

dev.new()  
par(bty="n",las=2)
plot(x=xlim,y=ylim,ylim=ylim,xlim=xlim,type="n",
     main="Humlum et al in forsking.no: annual means",
     sub="forskning.no (replicationDemos)",
     xlab="time",ylab=paste("anomaly wrt",ref,"(deg C)"))  
grid()
par(las=0)
mtext("http://www.forskning.no/artikler/2012/mars/316178",side=4,col="grey")
trendleg <- rep(" ",length(strleg))
  
if (uah) {
  UAH <- data.frame(x=L(uahmsu$Year),y=annual(uahtemp,uahmsu$Year))
  lines(UAH,col="blue",type="b",pch=19,cex=1.2)
  uahtrend <- lm(y ~ x, data=UAH); uahstat <- summary(uahtrend)
  lines(UAH$x,predict(uahtrend),col="blue",lty=2)
  trendleg[1] <- paste(round(uahstat$coefficients[2]*10,2)," [",
                    round(2*uahstat$coefficients[4]*10,2),"]",sep="")
}
if (rss) {
  RSS <- data.frame(x=L(rssmsu$year),y=annual(rsstemp,rssmsu$year))
  lines(RSS,type="b",pch=19,cex=1.2)
  rsstrend <- lm(y ~ x, data=RSS); rssstat <- summary(rsstrend)
  lines(RSS$x,predict(rsstrend),col="black",lty=2)
  trendleg[2] <- paste(round(rssstat$coefficients[2]*10,2)," [",
                    round(2*rssstat$coefficients[4]*10,2),"]",sep="")  
}
if (giss) {
  x <- L(gissdate); y <- annual(gistemp,gissdate)
  GIS <- data.frame(x=x[(x >= xlim[1]) & (x <= xlim[2])],
                    y=y[(x >= xlim[1]) & (x <= xlim[2])])
  lines(GIS,col="red",type="b",pch=19,cex=1.2)
  gisstrend <- lm(y ~ x, data=GIS); gissstat <- summary(gisstrend)
  lines(GIS$x[is.finite(GIS$y)],predict(gisstrend),col="red",lty=2)
  trendleg[3] <- paste(round(gissstat$coefficients[2]*10,2)," [",
                    round(2*gissstat$coefficients[4]*10,2),"]",sep="")  
}
if (ncdc) {
  x <- L(NCDC$year); y <- annual(ncdctemp,NCDC$year)
  NCD <- data.frame(x=x[(x >= xlim[1]) & (x <= xlim[2])],
                    y=y[(x >= xlim[1]) & (x <= xlim[2])])
  lines(NCD,col="violet",type="b",pch=19,cex=1.2)
  ncdctrend <- lm(y ~ x, data=NCD); ncdcstat <- summary(ncdctrend)
  lines(NCD$x[is.finite(NCD$y)],predict(ncdctrend),col="violet",lty=2)
  trendleg[4] <- paste(round(ncdcstat$coefficients[2]*10,2)," [",
                    round(2*ncdcstat$coefficients[4]*10,2),"]",sep="")   
}
if (hadcrut3) {
  x <- L(crudate); y <- annual(crutemp,crudate)
  CRU <- data.frame(x=x[(x >= xlim[1]) & (x <= xlim[2])],
                    y=y[(x >= xlim[1]) & (x <= xlim[2])])  
  lines(CRU,col="darkgreen",type="b",pch=19,cex=1.2)
  crutrend <- lm(y ~ x, data=CRU); crustat <- summary(crutrend)
  lines(CRU$x[is.finite(CRU$y)],predict(crutrend),col="darkgreen",lty=2)
  trendleg[5] <- paste(round(crustat$coefficients[2]*10,2)," [",
                    round(2*crustat$coefficients[4]*10,2),"]",sep="")    
}


legend(xlim[1],ylim[2],paste(strleg[ileg],trendleg[ileg]),
       col=colleg[ileg],
       bty="n",lty=1,lwd=lwd,cex=0.75)
  
}
