Scafetta2010 <- function() {

y0 <- 0.000029*(seq(1850,2011,by=1/12) - 1850)^2 - 0.42

hadcrut <- read.table("http://climexp.knmi.nl/data/ihadcrut3_gl.dat",
                       col.names=c("year","month","value"))
y1 <- hadcrut$value; y1[y1 < -5] <- NA

gistemp <- read.table("http://climexp.knmi.nl/data/igiss_al_gl_m.dat")
y2 <- c(t(gistemp[,2:13])); y2[y2 < -5] <- NA

ncdc <- read.table("http://climexp.knmi.nl/data/incdc_gl.dat",
                    col.names=c("year","month","value"))
y3 <- ncdc$value; y3[y3 < -5] <- NA

s0 <- spectrum(y0)
s1 <- spectrum(y1[is.finite(y1)])
s2 <- spectrum(y2[is.finite(y2)])
s3 <- spectrum(y3[is.finite(y3)])


Rz <- sunspots()
y4 <- Rz$sunspotnumber
s4 <- spectrum(y4)

x11()
par(mfcol=c(2,1))
plot(hadcrut$year+(hadcrut$month-0.5)/12,hadcrut$value,type="l")
lines(seq(1850,2011,by=1/12),y0,lwd=2,col="grey40")
grid()

plot(1/(s1$freq*12),s1$spec,type="l",log="xy",lwd=2,
     xlim=c(1,100),ylim=c(1e-4,1e2))
grid()
lines(rep(9.1,2),c(1e-3,1e6),col="grey")
lines(rep(11,2),c(1e-3,1e6),col="grey")
lines(rep(22,2),c(1e-3,1e6),col="grey")
lines(rep(60,2),c(1e-3,1e6),col="grey")

lines(1/(s0$freq*12),s0$spec*1000,lwd=2,col="grey40")
lines(1/(s1$freq*12),s1$spec,lwd=2)
lines(1/(s2$freq*12),s2$spec,col="red",lwd=2)
lines(1/(s3$freq*12),s3$spec,col="blue",lwd=2)
lines(1/(s4$freq*12),s4$spec*0.0001,col="darkgreen",lwd=2)

Y1 <- hadcrut$value[1:1933] - y0; yr01 <- 1850
Y2 <- c(t(gistemp[,2:13]))[1:1573] - y0[1:1573]; yr03 <- 1880
Y3 <- ncdc$value[1:1573] - y0[1:1573]; yr02 <- 1880
x11()
par(mfcol=c(3,1))
nwin <- 8
plot(seq(yr01,2011,by=1/12),ma.filt(Y1,nwin*12),type="l",
     xlim=c(1850,2060),lwd=3,col="grey",main="HadCRUT3")
lines(seq(yr01,2011,by=1/12)+61.5,ma.filt(Y1,nwin*12))

plot(seq(yr02,2011,by=1/12),ma.filt(Y2,nwin*12),type="l",
     xlim=c(1850,2060),lwd=3,col="grey",main="GISTEMP")
lines(seq(yr02,2011,by=1/12)+61.5,ma.filt(Y2,nwin*12))

plot(seq(yr03,2011,by=1/12),ma.filt(Y3,nwin*12),type="l",
     xlim=c(1850,2060),lwd=3,col="grey",main="NCDC")
lines(seq(yr03,2011,by=1/12)+61.5,ma.filt(Y3,nwin*12))

x11()
nwin <- 8
calibr <- data.frame(y=hadcrut$value,
          t=hadcrut$year + (hadcrut$month-0.5)/12 -min(hadcrut$year)+1)
trend1 <- lm(y ~ t,data=calibr)
trend2 <- lm(y ~ t + I(t^2),data=calibr)
trend3 <- lm(y ~ t + I(t^2) + I(t^3),data=calibr)
trend4 <- lm(y ~ t + I(t^2) + I(t^3) + I(t^4),data=calibr)
nt <- min(length(y0),length(predict(trend1)))

a1 <- acf(trend1$residual)
a2 <- acf(trend2$residual)
a3 <- acf(trend3$residual)
a4 <- acf(trend4$residual)
wn <- rnorm(nt+1)
rn <- wn[2:length(wn)] + a1$acf[2]*wn[1:nt]
rn1 <- rn/sd(rn)*sd(y1,na.rm=TRUE) + mean(y1,na.rm=TRUE) +
        predict(trend1)[1:nt] - y0[1:nt];
rn1 <- ma.filt(rn1,nwin*12)

wn <- rnorm(nt+1)
rn <- wn[2:length(wn)] + a2$acf[2]*wn[1:nt]
rn2 <- rn/sd(rn)*sd(y1,na.rm=TRUE) + mean(y1,na.rm=TRUE) +
        predict(trend2)[1:nt] - y0[1:nt];
rn2 <- ma.filt(rn2,nwin*12)

wn <- rnorm(nt+1)
rn <- wn[2:length(wn)] + a3$acf[2]*wn[1:nt]
rn3 <- rn/sd(rn)*sd(y1,na.rm=TRUE) + mean(y1,na.rm=TRUE) +
        predict(trend3)[1:nt] - y0[1:nt];
rn3 <- ma.filt(rn3,nwin*12)

wn <- rnorm(nt+1)
rn <- wn[2:length(wn)] + a4$acf[2]*wn[1:nt]
rn4 <- rn/sd(rn)*sd(y1,na.rm=TRUE) + mean(y1,na.rm=TRUE) +
        predict(trend4)[1:nt] - y0[1:nt];
rn4 <- ma.filt(rn4,nwin*12)

plot(seq(yr01,2011,by=1/12),rn1,type="l",
     main="Different noise processes superimposed on different trends")
lines(seq(yr01,2011,by=1/12),rn2,col="red")
lines(seq(yr01,2011,by=1/12),rn3,col="blue")
lines(seq(yr01,2011,by=1/12),rn4,col="darkgreen")

x11()
wn <- rnorm(nt)
dw <- rep(0,nt)
for (i in 2:nt) dw[i] <- dw[i-1] + wn[i]
plot(dw,main="Random walk")
}
