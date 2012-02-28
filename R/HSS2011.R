# This is an R-script reproduces the results from
# Humlum et al (2011),
# 'Identifying natural contributions to late Holocene climate change',
# Global and Planetary Change, vol 79, p. 145-156
# doi: 10.1016/j.gloplacha.2011.09.005
#
# http://cran.r-project.org
# Rasmus E. Benestad

Humlum.et.al.2011 <- function() {
  
#url <- "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/greenland/summit/gisp2/isotopes/gisp2_temp_accum_alley2000.txt"
#l <- readLines(url)
#Age <- grep("Age",l)
#last <- grep("49.981",l)

data(gisp2,envir=environment())
#gisp2 <- read.table(url,skip=Age[2],
#                    nrows=last - Age[2],
#                    col.names=c("age","temp"),as.is=TRUE)
#attr(gisp2,'url') <- url
#attr(gisp2,'date') <- "07.01.2012"
#save(file="Debunking/data/gisp2.rda",gisp2)
attach(gisp2)
# Age: thousand years before present
x <- age[age < 4]*1000; y <- temp[age < 4]

# Graphics: the GISP2 data
plot(-age*1000,temp,type="l",col="grey",lwd=3,
     xlab="Time (yrs before present)",ylab="Temperature (C)",
     xlim=c(-10000,2000),ylim=c(-34,-25),
     main="GISP2 temperature",
     sub="Replicating Humlum et al. (2011) and extending")
grid()
mtext("doi: 10.1016/j.gloplacha.2011.09.005", side=4)

# Only include the last 4000 yrs
phi <- c(2804, 1186, 556)
s1 <- sin(2*pi*x/phi[1]); c1 <- cos(2*pi*x/phi[1])
s2 <- sin(2*pi*x/phi[2]); c2 <- cos(2*pi*x/phi[2])
s3 <- sin(2*pi*x/phi[3]); c3 <- cos(2*pi*x/phi[3])
calibr <- data.frame(y=y,s1=s1,c1=c1,
                     s2=s2,c2=c2,s3=s3,
                     c3=c3,x=x)

# Include the lasty 10,000 yrs
X <- age[(age >= 4) & (age < 10)]*1000
S1 <- sin(2*pi*X/phi[1]); C1 <- cos(2*pi*X/phi[1])
S2 <- sin(2*pi*X/phi[2]); C2 <- cos(2*pi*X/phi[2])
S3 <- sin(2*pi*X/phi[3]); C3 <- cos(2*pi*X/phi[3])
extrap1 <- data.frame(s1=S1,c1=C1,
                      s2=S2,c2=C2,s3=S3,
                      c3=C3,x=X)
extrap2 <- data.frame(s1=S1,c1=C1,
                      s2=S2,c2=C2,s3=S3,
                      c3=C3,x=rep(min(X),length(X)))
# Calibrate the Humlum et al, 2011 model
fit <- lm(y ~ s1 + c1 + s2 + c2 + s3 + c3 + x,data=calibr)
lines(-x,y,lwd=3)

# Reproduce Humlum et al.'s model
lines(-x,predict(fit),col="red",lwd=2)

# Add the extended prediction 
lines(-X,predict(fit,newdata=extrap1),col="red",lty=2)
lines(-X,predict(fit,newdata=extrap2),col="red",lty=2)

#dev2bitmap("humlumetal2011-extended.png",res=150)

results <- list(t.hss=-x,fit.hss=predict(fit),
                t.hol=-X,fit.hol1=predict(fit,newdata=extrap1),
                fit.hol2=predict(fit,newdata=extrap1),
                model=fit)
invisible(results)
}
