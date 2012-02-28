MM2004 <- function() {

# R.E. Benestad, Oslo, Norway, July 8 2004
# The Norwegian Meteorological Institute
# rasmus.benestad@met.no
#
# This is R-script 'mcitrickCR.R'
# R is an analysis environment/statistical package (a "GNU version of Splus") freely available from 
# http://cran.r-project.org. (R runs on Linux, Windows and Mac.)
# 
# downloaded data from http://www.uoguelph.ca/~rmckitri/research/gdptemp.html
# -> gdptemp03.dif
# Loaded data using Excel andsaving the data as text:
# -> gdptemp03.txt
# Data file edited using emacs: "D-" replaced by "E-" to avoid 'NA's.
stand <- function(x) {  # Standardize the series
  x <- (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  invisible(x)
}

print("R-script for the analysis in Benestad, R.E. (2004) ")
print("'Are temperature trends affected by economic activity? ")
print("Comment on McKitrick & Michaels', Climate Research CR 27:171-173")

data(gdptemp03,envir=environment()); mm <- gdptemp03
##mm <- read.table("~/data/gdptemp03.txt",header=TRUE,as.is=TRUE)
lat <- mm$LAT
ysrt <- order(lat)
mm$STREND <- as.numeric(mm$STREND[ysrt])
mm$MTREND <- as.numeric(mm$MTREND[ysrt])
mm$GDPGROW <- as.numeric(mm$GDPGROW[ysrt])
mm$COAL80 <- as.numeric(mm$COAL80[ysrt])/1000
mm$COALGROW <- as.numeric(mm$COALGROW[ysrt])
mm$ABSCOSLAT <- abs(cos(pi*mm$LAT/180)[ysrt])
mm$PRESS <- mm$PRESS[ysrt]
mm$DPRESS <- mm$DPRESS[ysrt]
mm$WATER <- mm$WATER[ysrt]
mm$POP <- mm$POP[ysrt]
mm$INC79 <- mm$INC79[ysrt]
mm$SURFMISS <- mm$SURFMISS[ysrt]
mm$SOVIET <- mm$SOVIET[ysrt]
mm$LIT79 <- mm$LIT79[ysrt]
mm$SCALE79 <- mm$INC79[ysrt] * mm$POP[ysrt]

N <- length(mm$STREND)
mod1 <- lm(STREND ~ PRESS + WATER + ABSCOSLAT + POP + SCALE79 + COAL80 + 
           COALGROW + INC79 + GDPGROW + SOVIET + SURFMISS + LIT79,data=mm)
print(summary(mod1))

ii <- is.element(mm$SOVIET,1)
mod2 <- step(mod1,trace=0)
print(summary(mod2))

print("Correlation between STREND & SOVIET:")
print(cor.test(mm$STREND,mm$SOVIET))
print("Correlation between STREND & GDPGROW:")
print(cor.test(mm$STREND,mm$GDPGROW))

#mod1 <- lm(STREND ~ PRESS + DPRESS + WATER + ABSCOSLAT + POP + INC79 + GDPGROW + 
#           COAL80 + COALGROW + SURFMISS + SOVIET + LIT79,data=mm)
#
#print("Correlation between STREND & COAL80:")
#print(cor.test(mm$STREND,mm$COAL80))

# Random pick - not very good if the values are dependent...
#i <- rep(TRUE,N)
#rnd <- rnorm(N); isrt <- order(rnd)
#i[isrt[1:floor(N/2)]] <- FALSE
i <- c(rep(TRUE,floor(N/2)),rep(FALSE,ceiling(N/2)))

cal <- data.frame(STREND=mm$STREND[i], PRESS=mm$PRESS[i], DPRESS=mm$DPRESS[i], WATER=mm$WATER[i],
                  ABSCOSLAT=mm$ABSCOSLAT[i], POP=mm$POP[i], INC79=mm$INC79[i], GDPGROW=mm$GDPGROW[i], 
                  COAL80=mm$COAL80[i], COALGROW=mm$COALGROW[i], SURFMISS=mm$SURFMISS[i], 
                  SOVIET=mm$SOVIET[i], LIT79=mm$LIT79[i],SCALE79=mm$SCALE79[i],ind=(1:N)[i])
pre <- data.frame(STREND=mm$STREND[!i], PRESS=mm$PRESS[!i],DPRESS=mm$DPRESS[!i], WATER=mm$WATER[!i], 
                  ABSCOSLAT=mm$ABSCOSLAT[!i], POP=mm$POP[!i], INC79=mm$INC79[!i], GDPGROW=mm$GDPGROW[!i], 
                  COAL80=mm$COAL80[!i], COALGROW=mm$COALGROW[!i], SURFMISS=mm$SURFMISS[!i],
                  SOVIET=mm$SOVIET[!i], LIT79=mm$LIT79[!i],SCALE79=mm$SCALE79[!i],ind=(1:N)[!i])
mod1 <- lm(STREND ~ PRESS + WATER + ABSCOSLAT + POP + SCALE79 + COAL80 + COALGROW + INC79 + GDPGROW + SOVIET + SURFMISS + LIT79,data=cal)
mod2 <- lm(STREND ~ PRESS + WATER + ABSCOSLAT,data=cal)
mod3 <- lm(STREND ~ PRESS + DPRESS + WATER + ABSCOSLAT + POP,data=cal)
mod4 <- lm(STREND ~ POP + SCALE79 + INC79 + GDPGROW + COAL80 + COALGROW + SURFMISS + SOVIET + LIT79,data=cal)
mod5 <- lm(STREND ~ COAL80 + INC79 + GDPGROW + SOVIET + LIT79,data=cal)

#pre1 <- predict(mod1,newdata=pre)
#pre2 <- predict(mod2,newdata=pre)
#pre3 <- predict(mod3,newdata=pre)
#pre4 <- predict(mod4,newdata=pre)
#pre5 <- predict(mod5,newdata=pre)
pre1 <- predict(step(mod1,trace=0),newdata=pre)
pre2 <- predict(step(mod2,trace=0),newdata=pre)
pre3 <- predict(step(mod3,trace=0),newdata=pre)
pre4 <- predict(step(mod4,trace=0),newdata=pre)
pre5 <- predict(step(mod5,trace=0),newdata=pre)

x11()   # use 'windows()' on window systems
par(cex.main=0.8)
plot(mm$STREND,ylab="Trend estimates (degC/decade)",xlab="Station index",col="grey70",pch=20,
     main="Prediction using a subsample for calibration and remainding data for validation",
     sub="Sorted according to latitude & using an AIC-based step-wise model",ylim=c(-3,2),xlim=c(0,230))
grid()
points(cal$ind,predict(mod1),pch=21,col="blue")
points(pre$ind,pre1,pch=21,col="red",cex=0.8)
points(pre$ind,pre2,pch=2,col="grey40",cex=0.8)
points(pre$ind,pre3,pch=3,col="grey50",cex=0.8)
points(pre$ind,pre4,pch=4,col="darkgreen",cex=0.8)
points(pre$ind,pre5,pch=5,col="black",cex=0.8)
legend(10,-1.8,c("STREND","Calibration","Indep. all","PRESS+WATER+ABSCOSLAT","PRESS+DPRESS+WATER+ABSCOSLAT+POP",
                  "POP+SCALE79+INC79+GDPGROW+COAL80+COALGROW+SURFMISS+SOVIET+LIT79","COAL80+INC79+GDPGROW+SOVIET+LIT79"),
       pch=c(20,21,21,2,3,4,5),col=c("grey70","blue","red","grey40","grey50","darkgreen","black"),cex=0.75)
dev2bitmap(file="mckitric3.jpg",res=200,height=6,width=6,type="jpeg")

x11()   # use 'windows()' on window systems
par(cex.main=0.8)
plot(mm$STREND,ylab="Trend estimates (degC/decade)",xlab="Station index",col="grey70",pch=20,
     main="Prediction using a subsample for calibration and remainding data for validation",
     sub="Sorted according to latitude & using an AIC-based step-wise model",ylim=c(-3,2),xlim=c(0,230))
grid()
points(cal$ind,predict(mod1),pch=21,col="black")
points(pre$ind,pre1,pch=21,col="grey60",cex=0.8)
points(pre$ind,pre2,pch=2,col="grey40",cex=0.8)
points(pre$ind,pre3,pch=3,col="grey50",cex=0.8)
points(pre$ind,pre4,pch=4,col="grey10",cex=0.8)
points(pre$ind,pre5,pch=5,col="black",cex=0.8)
legend(10,-1.8,c("STREND","Calibration","Indep. all","PRESS+WATER+ABSCOSLAT","PRESS+DPRESS+WATER+ABSCOSLAT+POP",
                  "POP+SCALE79+INC79+GDPGROW+COAL80+COALGROW+SURFMISS+SOVIET+LIT79","COAL80+INC79+GDPGROW+SOVIET+LIT79"),
       pch=c(20,21,21,2,3,4,5),col=c("grey70","black","grey60","grey40","grey50","grey10","black"),cex=0.75)
dev.copy2eps(file="mckitric-retrial.eps")

}
