FL1991 <- function(dmi=FALSE,url="http://web.dmi.dk/fsweb/solarterrestrial/sunclimate/SCL.txt") {
  if (dmi) {
    L <- readLines(url)
    s1 <- grep("Solar cycle data for minima",L)+2
    meta <- L[1:s1]
    s2 <- grep("-----------",L)
    col.names <- c("Epoch","Notes","L","L121","L12221","Cntr.Year")
    widths <- c(8,4,6,6,7,10)
    scl1=read.fwf(url,skip=s1,nrows=s2[1]-s1-1,
      col.names=col.names,widths=widths,as.is=TRUE)
    s3 <- grep("Solar cycle data for maxima",L)
    scl2=read.fwf(url,skip=s2[1],nrows=s3-s2[1]-1,
      col.names=col.names,widths=widths,as.is=TRUE)
    scl3 <- read.fwf(url,skip=s3+1,nrows=s2[2]-s3-7,
      col.names=col.names,widths=widths,as.is=TRUE)
    scla <- cbind(scl1,rep("tab1",19))
    sclb <- cbind(scl2,rep("tab2",25))
    colnames(scla) <- c(col.names,'table')
    colnames(sclb) <- c(col.names,'table')
    fl1991 <- rbind(scla,sclb)
    attr(fl1991,'meta') <- meta
    attr(fl1991,'url') <- url
    attr(fl1991,'solmax') <- scl3
  } else {
    data(fl1991,envir=environment())
  }

  plot(fl1991$Cntr.Year,fl1991$L,type="b",pch=19)
  lines(fl1991$Cntr.Year,fl1991$L121,lwd=2,col="grey")
  lines(fl1991$Cntr.Year,fl1991$L12211,lwd=2,col="pink")
  
  invisible(fl1991)
}
