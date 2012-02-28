paleaoproxy <- function() {
# R.E. Benestad, Oslo 12.05.2005

stand <- function(x) {
 x <- (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
}
# Read the data from URL: 

url1 <- "ftp://cdiac.ornl.gov/pub/trends/co2/vostok.icecore.co2"
url2 <- "http://cdiac.esd.ornl.gov/ftp/trends/temp/vostok/vostok.1999.temp.dat"
url3 <- "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/greenland/summit/gisp2/cosmoiso/ber10.txt"
#a <- readLines(url1)

# Save the data in a local file: discard the header..
#a <- a[22:length(a)]
#writeLines(a,"vostoc_co2.dat")

#                Mean
#       Age of   age of    CO2
#Depth  the ice  the air concentration
# (m)   (yr BP)  (yr BP)  (ppmv)

#co2 <- read.table("vostoc.co2.dat",header=T,
#                  col.names=c("Depth","age.of.ice","age.of.air","co2.concentration"))
#attr(co2$Depth,"unit") <- "m"   
#attr(co2$age.of.ice,"unit") <- "yr BP"  
#attr(co2$age.of.air,"unit") <- "yr BP"  
#attr(co2$co2.concentration,"unit") <- "ppmv"
data(vostoc.co2,envir=environment())

# Save the data in a local file: discard the header..
#a <- readLines(url2)
#a <- a[59:length(a)]
#writeLines(a,"vostoc_tas.dat")

#                  Deuterium               
#         Age of     content   Temperature
# Depth   the ice   of the ice  Variation
#  (m)    (yr BP)   (delta D)    (deg C)

#tas <- read.table("vostoc_tas.dat",header=T,
#                  col.names=c("Depth","age.of.ice","deuterium","temperature"))
#attr(tas$Depth,"unit") <- "m"   
#attr(tas$age.of.ice,"unit") <- "yr BP"  
#attr(tas$deuterium,"unit") <- "delta D"  
#attr(tas$temperature,"unit") <- "deg C"
data(vostoc.temp,envir=environment())

#print("GISP2")
#be.10 <- read.table(url3,skip=41,header=FALSE,
# col.names=c("depth.top","depth.bottom","Be.10","error",
# "age.top","age.bottom"))
#attr(be.10$depth.top,"unit") <- "m"  
#attr(be.10$depth.bottom,"unit") <- "m"  
#attr(be.10$Be.10,"unit") <- "10^3 atom/g"  
#attr(be.10$error,"unit") <- "10^3 atom/g"  
#attr(be.10$age.top,"unit") <- "yr"  
#attr(be.10$age.bottom,"unit") <- "yr"  
#be.10$Be.10[be.10$Be.10>=999999] <- NA
#be.10$error[be.10$error>=999999] <- NA
data(Be.10,envir=environment())
be.age <- 0.5*(Be.10$age.top+Be.10$age.bottom)
#plot(co2$age.of.ice, co2$age.of.air)
#
#par(ask=TRUE)
#
#i1 <- is.element(co2$age.of.ice,tas$age.of.ice)
#i2 <- is.element(tas$age.of.ice,co2$age.of.ice)
#CO2 <- co2$co2.concentration[i1]
#T <- tas$temperature[i2]
#plot(1,1,type="n")
#acf(cbind(CO2,T))

plot(-vostoc.co2$age.of.ice, stand(vostoc.co2$co2.concentration), type="l",lwd=4,col="grey",
     xlab=attr(vostoc.co2$age.of.ice,"unit"),ylab="Standardised",ylim=c(-2,3),
     main="Historical CO2 Record from the Vostok Ice Core",sub=url1)
grid()
lines(-be.age,stand(Be.10$Be.10),lty=3,col="blue")
points(-be.age,stand(Be.10$Be.10),pch=2,col="blue",cex=0.5)
lines(-vostoc.co2$age.of.ice, stand(vostoc.co2$co2.concentration), lty=2,lwd=4,col="grey")
lines(-vostoc.temp$age.of.ice, stand(vostoc.temp$temperature),lty=2,lwd=2)
legend(-4e5,-1.5,c("CO2","TAS","Be-10  "),col=c("grey","black","blue"),
       lwd=c(4,2,1),lty=c(1,2,3),pch=c(NA,NA,2),cex=0.8)
dev2bitmap("paleaoproxy.png")


mills <- seq(1000*ceiling(min(-be.age)/1000),
             1000*ceiling(min(max(c(-be.age,-vostoc.co2$age.of.ice,-vostoc.temp$age.of.ice)))/1000),by=1000)
mills.be <-  1000*ceiling(-be.age/1000)
mills.co2 <- 1000*ceiling(-vostoc.co2$age.of.ice/1000)
mills.tas <- 1000*ceiling(-vostoc.temp$age.of.ice/1000)
y <- mills*NA; x1 <- y; x2 <- y
for (i in 1:length(mills)) {
  y[i]  <- mean(vostoc.temp$temperature[is.element(mills.tas,mills[i])],na.rm=TRUE)
  x1[i] <- mean(vostoc.co2$co2.concentration[is.element(mills.co2,mills[i])],na.rm=TRUE)
  x2[i] <- mean(Be.10$Be.10[is.element(mills.be,mills[i])],na.rm=TRUE)
}

x11()
plot(mills,stand(x1),type="l",lwd=4,col="grey",
     xlab=attr(vostoc.co2$age.of.ice,"unit"),ylab="Standardised millenium mean values",ylim=c(-2,3),
     main="Historical CO2 Record from the Vostok Ice Core",sub=url1)
lines(mills,stand(x2),lty=3,col="blue")
lines(mills,stand(y),lty=2,lwd=2)

good <- is.finite(y) & is.finite(x1) & is.finite(x2)
y <- y[good]; x1 <- x1[good]; x2 <- x2[good]

legend(-40000,3,c(paste("CO2",round(cor(y,x1),2)),"TAS    ",paste("Be-10",round(cor(y,x2),2))),
       col=c("grey","black","blue"),
       lwd=c(4,2,1),lty=c(1,2,3),pch=c(NA,NA,2),cex=0.8)
dev2bitmap("paleaoproxy_corr.png")

print("Correlation: TAS & CO2:")
print(cor.test(y,x1))
print("Correlation: TAS & Be-10:")
print(cor.test(y,x2))


#par(ask=FALSE)
}
