ENSO.example <- function(interval=1980:1989) {

print("Demonstration why the argument in Scafetta (2011), p. 9 section 3")
print(" 'harmonic model calibrated during the period 1950-2010 is able to ")
print(" carefully reconstruct the decadal and ,ultidecadal oscillations of the")
print(" temperature record from 1850-1950.' does not distinguish from a mere ")
print(" curve fitting exercise: here shown on the ENSO Nino3.4 index:'")
data(enso,envir=environment())
#load("Debunking/data/enso.rda")

x=enso$year + (enso$month-0.5)/12
y=enso$nino3.4.3month

ii <- is.element(enso$year,interval)
s <- fft(y[ii])
n <- length(s)/2; f <- 0.5*seq(0,1,length=n)*12  
tau <- 1/f
#plot(tau,abs(s[1:n]),type="l")

srt <- order(abs(s[1:n]),decreasing = TRUE)
w <- abs(tau[srt])
print(w)

ii1 <- is.element(round(x),seq(min(round(x[ii])),median(round(x[ii])),by=1))
ii2 <- is.element(round(x),round(seq(median(round(x[ii])),max(round(x[ii])),by=1)))
cal1 <- data.frame(x=enso$year[ii1] + (enso$month[ii1]-0.5)/12, 
                  y=enso$nino3.4.3month[ii1],
                  c1=cos(2*pi/w[1]*x[ii1]), s1=sin(2*pi/w[1]*x[ii1]),
                  c2=cos(2*pi/w[2]*x[ii1]), s2=sin(2*pi/w[2]*x[ii1]),
                  c3=cos(2*pi/w[3]*x[ii1]), s3=sin(2*pi/w[3]*x[ii1]))
cal2 <- data.frame(x=enso$year[ii2] + (enso$month[ii2]-0.5)/12, 
                  y=enso$nino3.4.3month[ii2],
                  c1=cos(2*pi/w[1]*x[ii2]), s1=sin(2*pi/w[1]*x[ii2]),
                  c2=cos(2*pi/w[2]*x[ii2]), s2=sin(2*pi/w[2]*x[ii2]),
                  c3=cos(2*pi/w[3]*x[ii2]), s3=sin(2*pi/w[3]*x[ii2]))
cfit1 <- lm(y ~ c1 + s1 + c2 + s2+ c3 + s3,data=cal1)
cfit2 <- lm(y ~ c1 + s1 + c2 + s2+ c3 + s3,data=cal2)

#print(summary(cal1)); print(summary(cal2))



dev.new()
plot(x,y,type="l", lwd=2, col="grey20", main="NINO3.4 curve fit",xlim=range(x[ii]),
     sub="Scafetta (2011): ...distinguishes a mere curve fitting exercise.. p. 9")
lines(x[ii1],predict(cfit1),lty=1,col="red")
lines(x[ii2],predict(cfit1,newdata=cal2),lty=2,col="red")

lines(x[ii2],predict(cfit2),lty=1,col="blue")
lines(x[ii1],predict(cfit2,newdata=cal1),lty=2,col="blue")

legend(1980,-1.75,c("fitted","predicted"),lty=c(1,2),bg="grey95")

cal <- data.frame(x=enso$year[ii] + (enso$month[ii]-0.5)/12, 
                  y=enso$nino3.4.3month[ii],
                  c1=cos(2*pi/w[1]*x[ii]), s1=sin(2*pi/w[1]*x[ii]),
                  c2=cos(2*pi/w[2]*x[ii]), s2=sin(2*pi/w[2]*x[ii]),
                  c3=cos(2*pi/w[3]*x[ii]), s3=sin(2*pi/w[3]*x[ii]))
ext <- data.frame(x=enso$year + (enso$month-0.5)/12, 
                  y=enso$nino3.4.3month,
                  c1=cos(2*pi/w[1]*x), s1=sin(2*pi/w[1]*x),
                  c2=cos(2*pi/w[2]*x), s2=sin(2*pi/w[2]*x),
                  c3=cos(2*pi/w[3]*x), s3=sin(2*pi/w[3]*x))
cfit <- lm(y ~ c1 + s1 + c2 + s2+ c3 + s3,data=cal)

dev.new()
plot(x,y,type="l", lwd=2, col="grey40", main="NINO3.4: extending curve-fit (Scafetta, 2011, p. 8)",
     sub="Demonstrates how the curve-fit fails, even if it matches part of data.")
lines(x[ii],y[ii], lwd=2, col="grey20")
lines(x[ii],predict(cfit),col="red")
lines(x,predict(cfit,newdata=ext),lty=2,col="red")

}

