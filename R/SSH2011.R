# For the Walker test: see 
Walker.test <- function(N,alpha=0.05) 1 - (1- alpha)^(1/N)

DJF <- function(obs,ssh2011.tab1=NULL) {
# Manually obtained from eKlima (14.03.2012):
# Svalbard Longyearbyen  
#  St.no   Month     TAM    TAMA     TAN    TANM     TAX    TAXM
#     99840 12.2011    -6.6     6.8   -19.8    -9.4     3.1    -3.8 
#     99840 01.2012    -3.4    11.9   -15.0    -6.2     4.8    -0.9          
#     99840 02.2012    -5.6    10.6   -19.1    -8.3     7.0    -2.5
  eklima <- mean(c(-6.6,-3.4,-5.6))
  X <- as.matrix(obs)[,3:14]
  I <- 1:length(obs$year)
  ii <- I[is.element(obs$year,obs$year[-1])]
  iii <- I[is.element(obs$year,obs$year[-1]-1)]
  #print(ii); print(iii)
  djf <- c(rowMeans(cbind(X[iii,12],X[ii,c(1,2)])),eklima)
  attr(djf,'year') <- c(obs$year[iii],2011)
  plot(attr(djf,'year'),djf,type="b",pch=19)
  polygon(c(2008,2020,2020,2008,2008),
          c(-20.5,-20.5,-14,-14,-20.5),
          col="grey",border="grey")
  lines(rep(2008,2),c(-50,50),col="red",lty=2)
  lines(c(2008,2020),rep(-17.2,2),lwd=3,col="grey30")
  if (!is.null(ssh2011.tab1)) {
    if (ssh2011.tab1$yr1[4]==1924) ssh2011.tab1$yr1[4] <- 1934
                                        # The year in Table 1 must be wrong
    n <- length(ssh2011.tab1$djf)
    for (i in 1:n) {
      lines(c(ssh2011.tab1$yr1[i],ssh2011.tab1$yr2[i]),
            rep(ssh2011.tab1$djf[i],2),lty=3,col="red")
    }
    lines(attr(djf,'year'),filter(djf,rep(1,7)/7),col="pink",lwd=3)
    # SCL=11 years: -17.2 = ( M(2008,2009,2010,2011)*4 + M(2012:2018)*7 ) /11
    M <- round(mean(djf[is.element(attr(djf,'year'),2008:2011)]),2)
    x <- (-17.2*11/7 - M*4/7 )
    print(paste("mean over 2008-2011=",M,"expect over 2008-2018 is -17.2",
                "hence the mean of 2012-2018 must be",x))
    lines(c(1900,2018),rep(x,2),lty=3,col="red")
  }
  points(attr(djf,'year'),djf,type="b",pch=19)
  invisible(djf)
}

obs2tab1 <- function(obs,tab,winter="djf") {
# Turn observations into table similar to SSH 2011 Tab1:
  print(summary(tab))
  attach(tab)
  X <- as.matrix(obs)[,3:14]
  #print(dim(X))
  maat[] <- NA; djf[] <- NA; mam[] <- NA; jja[] <- NA; son[] <- NA
  s.maat <- maat; s.djf <- maat; s.mam <- maat; s.jja <- maat
  s.son <- maat; n <- maat
  NN <- length(yr.min)
  I <- 1:length(obs$year)
   for (i in 1:NN) {
    #print(c(i,yr1[i],yr2[i]))
    if (is.finite(yr1[i]) & is.finite(yr2[i])) {
       ii <- I[is.element(obs$year,yr1[i]:yr2[i])]
      iii <- I[is.element(obs$year,(yr1[i]:yr2[i])-1)]
      maat[i] <- mean(X[ii,],na.rm=TRUE)
      mam[i] <- mean(X[ii,3:5],na.rm=TRUE)
      jja[i] <- mean(X[ii,6:8],na.rm=TRUE)
      son[i] <- mean(X[ii,9:11],na.rm=TRUE)
      n[i] <- sum(is.finite(rowMeans(X[ii,])))
      s.maat[i] <- sd(rowMeans(X[ii,]),na.rm=TRUE)
      s.mam[i] <- sd(rowMeans(X[ii,3:5]),na.rm=TRUE)
      s.jja[i] <- sd(rowMeans(X[ii,6:8]),na.rm=TRUE)
      s.son[i] <- sd(rowMeans(X[ii,9:11]),na.rm=TRUE)
      if (winter=="jfd") {
        djf[i] <- mean(X[ii,c(1,2,12)],na.rm=TRUE)
        s.djf[i] <- sd(rowMeans(X[ii,c(1,2,12)]),na.rm=TRUE)
      } else {
         djf[i] <- mean(cbind(X[iii,12],X[ii,c(1,2)]),na.rm=TRUE)
         s.djf[i] <- sd(rowMeans(cbind(X[iii,12],X[ii,c(1,2)])),na.rm=TRUE)
      }
    }
  }
  
  Tab <- tab
  #print(summary(Tab))
  Tab$maat <- round(maat,2); Tab$djf <- round(djf,2)
  Tab$mam <- round(mam,2); Tab$jja <- round(jja,2); Tab$son <- round(son,2)
  attr(Tab,"sigma.maat") <- s.maat
  attr(Tab,"sigma.djf") <- s.djf
  attr(Tab,"sigma.mam") <- s.mam
  attr(Tab,"sigma.jja") <- s.jja
  attr(Tab,"sigma.son") <- s.son
  attr(Tab,"n") <- n
  detach(tab)
  invisible(Tab)

  attr(Tab,"description") <- paste("obs2tab1:",attr(obs,"description"))
  invisible(Tab)
}

check.table1 <- function(ssh2011.tab1=NULL) {
  # Svalbard temperature.
  if (is.null(ssh2011.tab1)) {
    print("use the Svalbard temperature:")
    data(ssh2011.tab1,envir=environment())
    print(summary(ssh2011.tab1))
    #load("Debunking/data/ssh2011.tab1.rda")
    ssh2011.tab1$yr1[4] <- 1934 # The year in Table 1 must be wrong
  }
  #load("Debunking/data/svalbard.rda")
  data(svalbard,envir=environment())
  obs2tab1(svalbard,ssh2011.tab1) -> Ssh2011.Tab1
  n <- sqrt(attr(Ssh2011.Tab1,'n'))
  attach(ssh2011.tab1)
  plot(range(yr1,yr2,na.rm=TRUE),range(djf,jja,na.rm=TRUE),type="n",
       main="Check Table 1",xlab="Time",ylab="Temperature")
  NN <- length(yr1)
  for (i in 1:NN) {
   if (is.finite(yr1[i]) & is.finite(yr2[i])) {
     rect(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$maat[i]-attr(Ssh2011.Tab1,'sigma.maat')/n,
          Ssh2011.Tab1$yr2[i],Ssh2011.Tab1$maat[i]+attr(Ssh2011.Tab1,'sigma.maat')/n,
          col="grey80",border="grey70",lty=1)
     rect(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$djf[i]-attr(Ssh2011.Tab1,'sigma.djf')/n,
          Ssh2011.Tab1$yr2[i],Ssh2011.Tab1$djf[i]+attr(Ssh2011.Tab1,'sigma.djf')/n,
          col="grey80",border="grey70",lty=1)
     rect(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$jja[i]-attr(Ssh2011.Tab1,'sigma.jja')/n,
          Ssh2011.Tab1$yr2[i],Ssh2011.Tab1$jja[i]+attr(Ssh2011.Tab1,'sigma.jja')/n,
          col="grey80",border="grey70",lty=1)
     rect(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$mam[i]-attr(Ssh2011.Tab1,'sigma.mam')/n,
          Ssh2011.Tab1$yr2[i],Ssh2011.Tab1$mam[i]+attr(Ssh2011.Tab1,'sigma.mam')/n,
          col="grey80",border="grey70",lty=1)
     rect(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$son[i]-attr(Ssh2011.Tab1,'sigma.son')/n,
          Ssh2011.Tab1$yr2[i],Ssh2011.Tab1$son[i]+attr(Ssh2011.Tab1,'sigma.son')/n,
          col="grey80",border="grey70",lty=1)

     lines(c(yr1[i],yr2[i]),rep(maat[i],2),col="black",lwd=3)
     lines(c(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$yr2[i]),rep(Ssh2011.Tab1$maat[i],2),col="grey40",
           lwd=3,lty=2)
     lines(c(yr1[i],yr2[i]),rep(djf[i],2),col="blue",lwd=3)
     lines(c(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$yr2[i]),rep(Ssh2011.Tab1$djf[i],2),col="steelblue",
           lwd=3,lty=2)
     lines(c(yr1[i],yr2[i]),rep(mam[i],2),col="darkgreen",lwd=3)
     lines(c(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$yr2[i]),rep(Ssh2011.Tab1$mam[i],2),col="green",
           lwd=3,lty=2)
     lines(c(yr1[i],yr2[i]),rep(jja[i],2),col="wheat",lwd=3)
     lines(c(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$yr2[i]),rep(Ssh2011.Tab1$jja[i],2),col="yellow",
           lwd=3,lty=2)
     lines(c(yr1[i],yr2[i]),rep(son[i],2),col="darkred",lwd=3)
     lines(c(Ssh2011.Tab1$yr1[i],Ssh2011.Tab1$yr2[i]),rep(Ssh2011.Tab1$son[i],2),col="red",
           lwd=3,lty=2)
   }
  }
  
}


Solheim.et.al.2011 <- function(obs=NULL,ssh2011.tab1=NULL,N.tests=30000) {
#ssh2011.tab1 <- NULL; N.tests=1000  # Testing the script.
# R-package which provides the SCL results from Benestad (2005), GRL:  
require(cyclones)

# R-package providing the Durbin-Watson test:
require(lmtest)

# Retrieve the result from Table 1 in Solheim et al (2011): 
#ssh2011.tab1 <- read.table("SSH201b-table1.txt",header=TRUE)
if (is.null(ssh2011.tab1))  {
  #load("Debunking/data/ssh2011.tab1.rda")
  data(ssh2011.tab1,envir=environment())
  ssh2011.tab1$yr1[4] <- 1934 # The year in Table 1 must be wrong
}

# If other station data is provided, replace the Svalbard temperatures
# with new temperatures.
if (!is.null(obs)) {
  print("Use other station data")
  obs2tab1(obs,ssh2011.tab1) -> Ssh2011.Tab1
  rm(ssh2011.tab1)
  Ssh2011.Tab1 -> ssh2011.tab1
  print(ssh2011.tab1)
} else {
# Make a copy of the table and estimate sigmas simultaneously
  #load("Debunking/data/svalbard.rda")
  data(svalbard,envir=environment())
  obs2tab1(svalbard,ssh2011.tab1) -> Ssh2011.Tab1
}

attach(ssh2011.tab1)
NN <- length(scl)


#SCL from Benestad (2005) based on fitting harmonics and hence derived
# using all the sunspot data rather than just the sunhspot numbers around
# the time os solar min.
SCL <- Benestad2005()
scl.grl <- SCL$scl.max[SCL$yymm.max[,1]>1900,1]
yr.grl <- SCL$yymm.max[SCL$yymm.max[,1]>1900,1]
scl.grl <- scl.grl[is.finite(scl.grl)]
yr.grl <- yr.grl[is.finite(yr.grl)]

R <- sunspots(plot=FALSE)

# Plot the sunspot numbers, the SCL, and the start/end of solar cycles:
dev.new()
plot(R$year,R$sunspotnumber,type="n",ylab="Sunspots",
     xlab="Time",
     main="SCL, GCR, 10.7cm flux & sunspot number",
     xlim=c(1900,max(R$year)))
grid()
polygon(c(R$year,max(R$year),min(R$year)),c(R$sunspotnumber,0,0),
        col="grey70",border="grey70",density=30)

for (i in 1:NN) {
  if (is.finite(yr1[i]) & is.finite(yr2[i])) {
   lines(c(yr1[i],yr2[i]),10*rep(scl[i],2)+50,col="red",lwd=2)
   lines(rep(yr.min[i],2),c(0,300),col="red",lty=2)
 }
}

M <- length(scl.grl)-1
for (i in 1:M) {
   lines(c(yr.grl[i],yr.grl[i+1]),
         10*rep(scl.grl[i],2)+50,col="blue",lwd=2)
   lines(rep(yr.grl[i],2),c(0,300),col="blue",lty=2)
   lines(rep(yr.grl[i+1],2),c(0,300),col="blue",lty=2)
}

# Compare SCL 
dev.new()
plot(yr.min,scl,pch=19,ylim=c(18,6),type="b",
     maion="SCL from SSH2011, B2005, F-C&L1991")
points(yr.grl,scl.grl,pch=19,col="red",type="b")
grid()
mtext(attr(ssh2011.tab1,"description"),side=4,cex=0.8,col="grey")

data(fl1991,envir=environment())
i1 <- (fl1991$Cntr.Year > 1900)
y <- fl1991$L[i1] # plot with same year as SSH
points(yr.min,y,type="b",pch=21,lty=2,col="blue")
legend(1980,16,c("SSH2011","B2005","F-C&L1991"),
       col=c("black","red","blue"),pch=c(19,19,21),
       lty=c(1,1,2),bg="grey95")

# Estimate errors based on differences in SCL estimates:
scl.err <- scl[is.finite(scl)] - scl.grl
yr.err <- yr.min[is.finite(scl)] - yr.grl
scl.sd <- sd(scl.err)
yr.sd <- sd(yr.err)
# Plot 2 x st.dev. for ~95% conf. int.
for (i in 1:length(scl.err)) {
  rect(yr.min-2*yr.sd,scl-2*scl.sd,yr.min+2*yr.sd,scl+2*scl.sd,
       lty=2,border="grey")
}

# SSH lag the SCL record with respect to the solar cycle number in Fig. 1
x <- yr.min[2:NN]; y <- scl[1:(NN-1)]

# Reproduce Solheim et al. (2011) Fig. 1:
dev.new()
plot(x,y,xlim=c(1900,2020),ylim=c(15,8),
     pch=19,type="b",main="Solar Cycle Length",
     ylab="Length (yr)", xlab="Year")
grid()
mtext(attr(ssh2011.tab1,"description"),side=4,cex=0.8,col="grey")
t <- seq(-1,1,length=length(y))
scl.trend <- lm(y ~ t)
abline(scl.trend,lty=3,col="grey")
N <- length(t)-1
y <- scl.trend$residual[1:N]

# De-trend to get a more objective correlation estimate:
# The coincidence of long-term trends which are unrelated
# is non-negligible and may affect the correlation estimate:
maat <- maat[1:(NN-1)]; djf <- djf[1:(NN-1)]
maat.trend <- lm(maat ~ t)
djf.trend <- lm(djf ~ t)
maat.dt <- maat.trend$residual
djf.dt <- djf.trend$residual

# Correlations:
z <- maat.dt
print("Correlation: de-trended MAAT & SCL")
print(cor(y,z))

z <- djf.dt
print("Correlation: de-trended DJF & SCL")
print(cor(y,z))

dev.new()
acf(y,main="detrended SCL autocorrelation") -> scl.ar
mtext(attr(ssh2011.tab1,"description"),side=4,cex=0.8,col="grey")

dev.new()
acf(z,main="detrended DJF autocorrelation") -> djf.ar
mtext(attr(ssh2011.tab1,"description"),side=4,cex=0.8,col="grey")

# Monte-Carlo simulations to estimate the null-distribution for correlation:
cor.null <- rep(NA,N.tests)
for (i in 1:N.tests) {
  # generate red noise with the same autocorrelation as SCL
  wn <- rnorm(length(z))
  rn <- wn
  for (ii in 2:length(rn)) rn[ii] <- scl.ar$acf[2]*rn[ii-1] +
       (1-abs(scl.ar$acf[2]))*rnorm(1)
  # Correlaiton random red noise
  cor.null[i] <- cor(z,rn)
}
p.cor <- round(100*sum( abs(cor.null) > abs(cor(y,z)) )/N.tests,2)
h.null <- hist(cor.null,breaks=seq(-1,1,length=30),plot=FALSE)

# Estimate the error bars for the correlation estimates based on errors
# estimated from differences with Benestad (2005)
# Also repeat the Durbin-Watson test through Monte-Carlo simulations
cor.err <- rep(NA,N.tests)
cor.ssh <- cor.err
dw.null <- rep(NA,N.tests)

print("Monte-Carlo simulations:")

# Compute a null-distributions for correlation and DW-test

for (i in 1:N.tests) {
  ## generate red noise with the same autocorrelation as SCL
  wn.s <- rnorm(length(y),mean=0,sd=scl.sd)
  wn.t <- rnorm(length(z),mean=rep(0,length(z),
                            sd=attr(Ssh2011.Tab1,'sigma.djf')/sqrt(attr(Ssh2011.Tab1,'n'))))
    
## perform two correlation tests
# based on error estimates of temp & SCL:  
  cor.err[i] <- cor(z+wn.t,y+wn.s)

# Repeat SSH strategy  
  scramble <- trunc(9*runif(9))+1
  cor.ssh[i] <- cor(z[scramble],y[scramble])

## generate red noise with the same autocorrelation as SCL for DW-test:
  wn <- rnorm(length(z))
  rn <- wn
  for (ii in 2:length(rn)) rn[ii] <- scl.ar$acf[2]*rn[ii-1] +
       (1-abs(scl.ar$acf[2]))*rnorm(1)
    
## perform Durbin-Watson test
  dwtest(z ~ rn) -> test.results
  dw.null[i] <- test.results$statistic
}

P.cor <- round(100*sum( abs(cor.err) > abs(cor(y,z)) )/N.tests,2)
p.ssh <- round(100*sum( abs(cor.ssh) > abs(cor(y,z)) )/N.tests,2)
h.err <- hist(cor.err,breaks=seq(-1,1,length=30),plot=FALSE)
h.ssh <- hist(cor.ssh,breaks=seq(-1,1,length=30),plot=FALSE)

print(paste("95%CI=",quantile(cor.err,0.025),"-",quantile(cor.err,0.975)))
print("Reported 95% CI in Solheim et al. (2011): 0.54--0.96")

dev.new(ces.sub=0.75)
plot(h.null$mids,h.null$density/max(h.null$density),type="n",
     ylim=c(0,max(h.null$density,h.err$density)),ylab="probability density",
     xlab="Correlation",main="Correlation analysis: de-trended SCL & DJF",
     sub=paste("Pr(R>r)=",p.cor,"r(SCL|SSH2011) > ",
       P.cor,"% of M-C errors &",p.ssh," % SSH2011 errors"))
polygon(c(h.null$mids,max(h.null$mids),min(h.null$mids)),
        c(h.null$density/max(h.null$density),0,0),
        col="grey70",border="grey70",lwd=3)
lines(h.err$mids,h.err$density/max(h.err$density),lwd=3)
lines(h.ssh$mids,h.ssh$density/max(h.ssh$density),lwd=1,lty=3,col="grey45")
lines(rep(cor(y,z),2),c(0,10),lwd=2,col="red",lty=2)
lines(rep(cor(scl.grl[2:(N+1)],z),2),c(0,10),lwd=2,col="blue",lty=2)
lines(rep(mean(cor.err),2),c(0,10),lwd=1,col="grey30",lty=2)

# Table 2 in SSH2011 djf r 95% conf. int.:
lines(rep(mean(cor.ssh),2),c(0,10),lwd=1,col="grey45",lty=3)
lines(rep(-0.97,2),c(0,10),lwd=1,lty=3)
lines(rep(-0.52,2),c(0,10),lwd=1,lty=3)
grid()
mtext(attr(ssh2011.tab1,"description"),side=4,cex=0.8,col="grey")
legend(0.3,1.7,c("Monte-Carlo","bootstrap","r(SSH2011)","r(B2005)",
                 "95%(SSH2011)"),
        lwd=c(3,1,2,2,1),col=c("black","grey45","red","blue","black"),
        lty=c(1,3,2,2,2),bg="grey95")

# Finished with correlation...

dwtest(z ~ y) -> test.results
p.dw <- sum( dw.null > test.results$statistic )/N.tests
h.dw <- hist(dw.null,plot=FALSE,breaks=seq(0,4,length=40))

dev.new()
plot(h.dw$mids,h.dw$density,lwd=3,type="l",
     main="Durbin-Watson",col="grey")
lines(rep(test.results$statistic,2),c(0,10),lwd=2,col="red")
mtext(attr(ssh2011.tab1,"description"),side=4,cex=0.8,col="grey")

# The Walker test to assess the field significance: 10 different tests
# (4 seasons + 1 annual mean) x (0 & 1 lag)

p.W <- round(Walker.test(10),4)
print(paste( "p(r)=",p.cor,"p(DW)=",p.dw,"  p(Walker,N=10)=",p.W) )
print(paste( "Field signf. correlation at 0.05% level:",(p.cor<p.W) ) )
print(paste( "Field signf. Durbin-Watson test at 0.05% level:",(p.dw<p.W) ) )
detach(ssh2011.tab1)
}

do.vardo <- function() {
  # Vardø temperature:
  #load("Debunking/data/vardo.rda")
  data(vardo,envir=environment())
  Solheim.et.al.2011(obs=vardo)
}
