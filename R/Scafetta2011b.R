hfit <- function(x,t,T1=60,T2=20) {
    cal <- data.frame(x=x,t=t,c1=cos(2*pi*t/T1),s1=sin(2*pi*t/T1),
                      c2=cos(2*pi*t/T2),s2=sin(2*pi*t/T2))
    hfit <- lm(x ~ t + c1 + s1 + c2 + s2,data=cal)
    hfit                   
}

Scafetta2011 <- function() {

  p <- function(t,P0=-0.30,P1=-0.0035, P2=0.000049) {
    pt <- P2*(t-1850)^2 + P1*(t - 1850) + P0
    pt
  }
  
  f <- function(t,C1=0.10, C2=0.04,T1=2000.8,T2=2000.8) {
    # Scafetta: eq. 3
    ft <- C1*cos(2*pi*(t-T1)/60) + C2*cos(2*pi*(t-T2)/20)
    ft
  }

  g <- function(t,C3=0.03,C4=0.05,T3=2002.7,T4=1997.7) {
    gt <- C3*cos(2*pi*(t-T3)/10.44) + C4*cos(2*pi*(t-T4)/9.07)
    gt
  }

  q <- function(t) {
    qt <- 0.009 * (t - 2000)
    qt
  }

  qerr <- function(t) {
    qerr <- rep(0,length(t))
    i <- (t > 2000)
    qerr[i] <- 0.004*(t[i] - 2000)
    qerr
  }
  
  h <- function(t,C3=0.03,T3=2002.7,C4=0.05,T4=1997.7) {
    ht <- rep(NA,length(t))
    i1 <- (t < 2000) & (t > 1850)
    i2 <- (t > 2000) & (t < 2100)
    ht[i1] <- f(t[i1]) + g(t[i1],C3=C3,C4=C4,T3=T3,T4=T4) + p(t[i1])
    ht[i2] <- f(t[i2]) + g(t[i2],C3=C3,C4=C4,T3=T3,T4=T4) + p(2000) + q(t[i2]) 
    ht
  }
    


  data(CRU,envir=environment())  
  yymm <- sort(rep(CRU[,1],12)) + (rep(1:12,length(CRU[,1]))-0.5)/12
  # Annual mean temperature:
  t2m <- c(t(CRU[,2:13])); N <- length(t2m)

  plot(yymm,t2m,type="l",col="red",xlim=c(1850,2050),ylim=c(-0.8,1.2),
       main="After Scafetta (2011): Fig. 5b",
       xlab="year",ylab="Temp. Anom (C)")
  lines(c(1850,2100),rep(0.6,2),lty=2,col="grey")
  grid()
  lines(yymm,ma.filt(t2m,4*12),lwd=3,col="grey40")

  t <- seq(1850,2100,by=1/12)
  lines(t,h(t,C3=0.03,T3=2003.0,C4=0.05,T4=1997.5),lwd=2,col="darkblue")
  lines(t,h(t,C3=0.04,T3=2002.1,C4=0.05,T4=1998.1),lwd=2)
  lines(t,h(t),lwd=2,lty=2,col="grey70")
  lines(t,h(t)+qerr(t),lty=2)
  lines(t,h(t)-qerr(t),lty=2)

  text(2000,1.20,"S2011: anomaly by",pos=4,cex=0.8)
  text(2000,1.09,"       2050= ~0.6C",pos=4,cex=0.8)
  text(1850,1.15,"f(t)= C1*cos(2*p*(t-T1)/60) + C2*cos(2*pi*(t-T2)/20)",
       pos=4,cex=0.8)
  text(1850,1.05,"g(t)= C3*cos(2*pi*(t-T3)/10.44) + C4*cos(2*pi*(t-T4)/9.07)",
       pos=4,cex=0.8)
  text(1850,0.90,"p(t) = P2*(t-1850)^2 + P1*(t - 1850) + P0",pos=4,cex=0.8)
  text(1850,0.80,"q(t) 0.009 * (t - 2000)",pos=4,cex=0.8)
  
  print("Test harmonic fit:")

  model <- hfit(t2m,yymm)
  print(summary(model))

  R2max <- 0; I <- NA; J <- NA
  for (i in seq(30,80,by=0.25)) {
    for (j in seq(10,30,by=0.25)) {
      model70.15 <- hfit(t2m,yymm,T1=i,T2=j)
      R2 <- summary(model70.15)$r.squared
      if (R2 > R2max) {
        R2max <- R2
        I <- i; J <- j
      }
    }
  }
  print(paste("Best harmonic fit: periodicities=",I,"&",J,"years"))
  text(1940,-0.7,paste("Best harmonic fit: periodicities=",I,"&",J,"years"),
       cex=0.8,pos=4,col="red")
  bestmodel <- hfit(t2m,yymm,T1=I,T2=J)
  print(summary(bestmodel))
  bestfit <- predict(bestmodel)
  
  
  lines(yymm,bestfit,col="green",lty=2,lwd=2)
}

Scafetta.tab1 <- function() {
  data(CRU,envir=environment())
  t0 <- sort(rep(CRU[,1],12)) + (rep(1:12,length(CRU[,1]))-0.5)/12

  # It seems they used annual data:
  t2m <- c(t(CRU[,2:13])); N <- length(t2m)  
  data(Scafetta2011.tab1,envir=environment())
  data(CMIP3.20c3m.sresa1b,envir=environment())
  CMIP3 <- CMIP3.20c3m.sresa1b
  d <- dim(CMIP3)
  A <- rep(NA,d[2]); B <- A; C <- A
  t <- attr(CMIP3,'year') + (attr(CMIP3,'month')-0.5)/12

  # Match the intervals and only use the interval when all have valid data:
  ensmean <- rowMeans(CMIP3)
  allok <- is.finite(ensmean)
  year <- attr(CMIP3,'year')[allok]
  month <- attr(CMIP3,'month')[allok]
  CMIP3 <- CMIP3[allok,]
  #print(table(year)); print(table(month)); stop("HERE")

  t <- t[allok]
  i0 <- is.element(t0,t)
  i1 <- is.element(t,t0)
  year <- year[i1];month=month[i1]
  t <- t[i1]; CMIP3 <- CMIP3[i1,]
  
  # Harmonic regression for CRU:
  hfit(t2m[i0],t0[i0],T1=60,T2=20) -> hmodel0
  fitstat0 <- summary(hmodel0)

  par(bty="n")
  plot(t0,t2m,pch=19,cex=0.75,col="grey",
       main="Global mean temperature: HadCRUT3v + CMIP3",xlab="time",
       ylab="Temperature anomaly from 1961-90 (deg C)",
       ylim=range(c(t2m*1.5),na.rm=TRUE))
  lines(t0[i0],predict(hmodel0),lwd=2)
  
  c0 <- fitstat0$coefficients[2]
  a0 <- sqrt(fitstat0$coefficients[3]^2 + fitstat0$coefficients[4]^2)
  b0 <- sqrt(fitstat0$coefficients[5]^2 + fitstat0$coefficients[6]^2)

  # Harmonic regression for GCMs:
  for (i in 1:d[2]) {
    y <- CMIP3[,i]
    for (im in 1:12) {
      ii <- is.element(year,1961:1990) & is.element(month,im)
      y.ref <- mean(y[ii],na.rm=TRUE)
      #print(c(sum(ii),y.ref,sum(is.element(month,im))))
      y[is.element(month,im)] <- y[is.element(month,im)] - y.ref
    }
    points(t,y,lty=2,col="pink",cex=0.5)
    hfit(y,t,T1=60,T2=20) -> hmodel
    lines(t,predict(hmodel),col="red")
    fitstat <- summary(hmodel)
    C[i] <- fitstat$coefficients[2]
    A[i] <- sqrt(fitstat$coefficients[3]^2 + fitstat$coefficients[4]^2)
    B[i] <- sqrt(fitstat$coefficients[5]^2 + fitstat$coefficients[6]^2)
 }
  points(t0,t2m,pch=19,cex=0.75,col="grey")
  lines(t0[i0],predict(hmodel0),lwd=2)

  
  results <- list(A=A,a=Scafetta2011.tab1$a[-1]*0.10,
                  B=B,b=Scafetta2011.tab1$b[-1]*0.040,
                  C=C*10,c=Scafetta2011.tab1$c[-1]*0.10)
  dev.new(); par(bty="n",col.axis="white")
  plot(c(0,7),c(-0.5,0.25),type="n",ylim=c(-0.05,0.25),
          main="Scafetta (2011): Table 1",xlab="coeffiecient",
          ylab="estimate (C1,C2)",col=c("wheat","grey"),lwd=2,
          sub=paste(min(trunc(t0[i0])),"-",max(trunc(t0[i0]))),
          cex=1.7)
  par(col.axis="black")
  grid()
  for (i in 1:6) rect(i-0.5,quantile(results[[i]],0.025),
                      i+0.5,quantile(results[[i]],0.975),
                      col="grey85",border="grey85")
  boxplot(results,add=TRUE,col=c("wheat","grey"))
  points(c(a0,Scafetta2011.tab1$a[1]*0.10,b0,Scafetta2011.tab1$b[1]*0.040,
           c0*10,Scafetta2011.tab1$c[1]*10/(2011-1850)),
         pch=c(19,21),col="red",cex=1.25)

  legend(0.4,0.25,c("Our replication","From Table 1"),
         lwd=10,col=c("wheat","grey"),bty="n")

  legend(3.4,0.25,c("Our estimate","From Table 1"),
         pch=c(19,21),col="red",bty="n")

  lines(c(0.5,2.5),rep(-0.03,2),lty=2,lwd=2,col="red")
  text(1.5,-0.04,"60 year",col="red")

  lines(c(2.5,4.5),rep(-0.03,2),lty=2,lwd=2)
  text(3.5,-0.04,"30 year")

  lines(c(4.5,6.5),rep(-0.03,2),lty=2,lwd=2,col="grey")
  text(5.5,-0.04,"trend",col="grey")

}
