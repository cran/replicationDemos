Douglas2007 <- function() {

  df2m <- function(X) {
# Convert the data.frame into a matrix:
    #print("df2m:")
    v <- names(X)[-(1:2)]
    d <- dim(X)
    #print(d)
    d[2] <- length(v)
    M <- matrix(rep(NA,d[1]*d[2]),d[1],d[2])
    for (i in 1:d[2])
      eval(parse(text=paste("M[,i]<-X$",v[i],sep="")))
    colnames(M) <- substr(v,2,nchar(v))
    rownames(M) <- X$runs
    #print("M:"); print(M)
    invisible(M)
  }
  

  p.hydrostatic <- function (h, p0 = 1000, Temp = 288, g = 9.81,
                             k = 1.38e-23, M = 0.027/6.022e+23) 
{
    p <- p0 * exp(-(M * g * h)/(k * Temp))
    p
}

  cat("Reproduction of results in Fig 1. of Douglas et al. (2007)")
  cat("'A comparison of tropical temperature trends with model predictions'")
  cat("INTERNATIONAL JOURNAL OF CLIMATOLOGY")
  cat("Int. J. Climatol. (2007)")
  cat("Published online in Wiley InterScience")
  cat("(www.interscience.wiley.com) DOI: 10.1002/joc.1651")
  cat("")
  cat("Based on Tables I & II in the paper. The values have been") 
  cat("copied from the on-line PDF through acroreader.")
  cat("(the negative sign of the values had to be set to '-')")

  data(Douglasetal.tab1,envir=environment())
  #load("Debunking/data/Douglasetal.tab1.rda")
  data(Douglasetal.tab2,envir=environment())
  #load("Debunking/data/Douglasetal.tab2.rda")
  X1 <- df2m(Douglasetal.tab1)/1000
  lev1 <- attr(Douglasetal.tab1,'levels')
  X2 <- df2m(Douglasetal.tab2)/1000
  dim(X2) <- c(22,13)
  #print(class(X2))
  lev2 <- attr(Douglasetal.tab2,'levels')
  plot(range(100,1000),c(-0.5,1.5),type="n",
       xlab="level",ylab="Temperature trend",
       xlim=c(950,100),ylim=c(-0.25,1.0),
       main="Trends from Table II [Douglas et al. 2007]")
  rect(1000,-0.2,100,0.4,col="grey95",border="grey95")
  grid()
  nobs <- dim(X1)[1]
  ngcm <- dim(X2)[1]
  nz <- dim(X2)[2]
  lightred <- rgb(1,0.5,0.5)

  for (i in 1:nz) lines(rep(lev2[i],2),range(X2[,i]),lwd=7,col="lightpink")
  for (i in 1:ngcm) {
    points(lev2,X2[i,],pch=19,type="b",lty=3,col=lightred)
  }

  Dm <- attr(Douglasetal.tab2,'Average')/1000
  Ds <- 2*attr(Douglasetal.tab2,'Std. Dev.')/(sqrt(21)*1000)
  lines(lev2,Dm,type="b",cex=1.2,lwd=2)
  lines(lev2,Dm+Ds,type="b",cex=0.7)
  lines(lev2,Dm-Ds,type="b",cex=0.7,pch=4)
  
  outside <- lev2*0
  for (i in 1:nz) outside[i] <- sum( (X2[i,]>Dm[i]+Ds[i]) |
                                     (X2[i,]<Dm[i]-Ds[i]) )
  text(1000,1,"GCMs outside mean +- 'uncertainty':",pos=4)
  for (i in 1:nz) text(lev2[i],0.9,outside[i],cex=1.2,font=2)
  for (i in 1:nz) text(lev2[i],0.8,paste(round(100*outside[i]/22),"%"),
                       cex=0.8,srt=90)

  for (i in 1:nobs) {
      points(lev1,X1[i,],pch=19,col="blue",type="b",lty=3)
  }

#print(Douglasetal.tab1)
#print(X1)
}
