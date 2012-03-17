# A reproduction of Lohle & Scafetta 2011
# Read and organise the data from cRU:

LoehleScafetta2011 <- function() {
  
#CRU <- read.table("CRU-T2m.mon.txt")
data(CRU,envir=environment())  
yymm <- sort(rep(CRU[,1],12)) + (rep(1:12,length(CRU[,1]))-0.5)/12

# It seems they used annual data:
t2m <- ma.filt(c(t(CRU[,2:13])),12); N <- length(t2m)

# Select the data before 1950
i1950 <- yymm < 1950
n <- sum(i1950)

# Construct similar harmonics
w20 <- pi/(10*12)*seq(1,sum(i1950),by=1)
w60 <- pi/(30*12)*seq(1,sum(i1950),by=1)
c20 <- cos(w20); s20 <- sin(w20)
c60 <- cos(w60); s60 <- sin(w60)
W20 <- pi/(10*12)*seq(1,N,by=1)
W60 <- pi/(30*12)*seq(1,N,by=1)
C20 <- cos(W20); S20 <- sin(W20)
C60 <- cos(W60); S60 <- sin(W60)
cal60 <- data.frame(y=t2m[i1950],x1=c60,x2=s60,x3=c20,x4=s20)
pre60 <- data.frame(y=t2m,x1=C60,x2=S60,x3=C20,x4=S20)

# Regression analysis & prediction
fit60m <- lm(y ~ x1 + x2 + x3 + x4,data=cal60)
fit60 <- predict(fit60m,newdata=pre60)

# Plot the data
plot(yymm,t2m,type="l",col="grey")
lines(yymm[i1950],t2m[i1950])
lines(yymm,fit60,col="red",lwd=2)
#dev2bitmap("loehlesscafetta2011redux.png")

# Explore other frequencies too, keeping the best-fit amplitudes:
dev.new()
plot(yymm,t2m,type="l",col="grey")
lines(yymm[i1950],t2m[i1950])
A <- rep(NA,50); A1 <- A; A2 <- A
for (i in 1:50) {
  wX <- pi/(i*12)*seq(1,sum(i1950),by=1)
  cX <- cos(wX); sX <- sin(wX)
  calX <- data.frame(y=t2m[i1950],x1=cX,x2=sX)
  fitXm <- lm(y ~ x1 + x2,data=calX)
  fitXs <- summary(fitXm)
  fitX <- predict(fitXm,newdata=calX)
  lines(yymm[i1950],fitX,col="blue",lwd=2)
  A[i] <- sqrt(fitXs$coefficients[2]^2 + fitXs$coefficients[3]^2)
  A1[i] <- sqrt((fitXs$coefficients[2]-fitXs$coefficients[5])^2 +
                (fitXs$coefficients[3]-fitXs$coefficients[6])^2)
  A2[i] <- sqrt((fitXs$coefficients[2]+fitXs$coefficients[5])^2 +
                (fitXs$coefficients[3]+fitXs$coefficients[6])^2)
}
lines(yymm,fit60,col="red",lwd=2)
#dev2bitmap("loehlesscafetta2011-manycycles.png")

# Plot the amplitudes:
dev.new()
plot(c(2,100),range(c(A1,A2)),type="n",lwd=2,
     main="Amplitude of best-fit cycle",
     xlab="Periodicity (Year)",ylab="Amplitude (C)")
polygon(c(seq(2,100,by=2),reverse(seq(2,100,by=2))),
        c(A1,reverse(A2)),col="grey",border="grey")
lines(seq(2,100,by=2),A,lwd=3)
grid()
#dev2bitmap("loehlesscafetta2011amplitudes.png")

#-----------------------
# Construct a synthetic time series: 10.000 year long
# Set up data matrix representing random amplitudes, and split between
# cosine and sine:
A <- rnorm(1000); dim(A) <- c(500,2)

# Set up matrix and vector holding the results; x holds the random curve and
# y the different cosine+sinus components.
x <- rep(0,10000*12); # y <- matrix(rep(x,5000),5000,10000*12)

# Define the basis frequency of the cosine and sine curves: the largest
# frequency corresponds to one cycle for the whole curve:
w <- pi*seq(0,2,length=10000*12)

# Carry out the construction of the random curves:
for (i in 1:500) {
#  y[i,] <- ( A[i,1]*cos(i*w) + A[i,2]*sin(i*w))/sqrt(i)
#  x <- x + y[i,]
  y <- ( A[i,1]*cos(i*w) + A[i,2]*sin(i*w))/sqrt(i)
  x <- x + y
}

x <- ma.filt(x,12)

dev.new()
plot(x,type="l",lwd=2,
     main="10.000yr synthetic series")
#dev2bitmap("loehlesscafetta2011syntheticseries.png",res=150)

dev.new()
plot(c(1,n),range(x,na.rm=TRUE),type="n",
     main=paste("Looking at",n/12,"year seqments"),
     xlab="month",ylab="")

sequence <- seq(1,10000*12,b=n)
ns <- length(sequence)
a1 <- rep(NA,ns); b1 <- a1
ii <- 1
for (i in sequence) {
  CAL60 <- data.frame(y=x[i:(i+n-1)],x1=c60,x2=s60,x3=c20,x4=s20)
# Regression analysis & prediction
  FIT60m <- lm(y ~ x1 + x2 + x3 + x4,data=CAL60)
  FIT60 <- predict(FIT60m,newdata=CAL60)
  lines(x[i:(i+n-1)],lwd=2,col="grey")
  lines(FIT60,col="red",lwd=2)
  a1[ii] <- summary(FIT60m)$coefficients[2]^2 +
            summary(FIT60m)$coefficients[3]^2
  b1[ii] <- summary(FIT60m)$coefficients[4]^2 +
            summary(FIT60m)$coefficients[5]^2
  ii <- ii + 1
}
#dev2bitmap("loehlescafetta2011sequences.png",res=150)

dev.new()
plot(sqrt(a1),type="l",lwd=3,ylim=c(0,3),
     main="Fitted 20yr & 60yr amplitudes",
     xlab="sequence",ylab="amplitude",
     sub="test: how the amplitudes vary between different sequences")
mtext("Same underlying harmonics for all sequences",side=4)
lines(sqrt(b1),lwd=3,col="blue")
legend(0,3,c("20year","60year"),col=c("blue","black"),lty=1,cex=0.8)
}
