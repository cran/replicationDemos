# This is an R-script that demonstrates how a random curve consists
# of a number of cosine and sine curves with random amplitudes.
# The script can be run in the R-environment, freely available from
# http://cran.r-project.org (Mac/Linux/Windows).
# To run, save the file as 'FT.R' in a path that R reads from
# and type 'source("FT:R")'
# Rasmus E. Benestad, Boulder, Aug. 14, 2001.

decomposeFT <- function(N=1000) {

# Set up data matrix representing random amplitudes, and split between
# cosine and sine:
A <- rnorm(100); dim(A) <- c(50,2)

# Set up matrix and vector holding the results; x holds the random curve and
# y the different cosine+sinus components.
x <- rep(0,N); y <- matrix(rep(x,50),50,N)

# Define the basis frequency of the cosine and sine curves: the largest
# frequency corresponds to one cycle for the whole curve:
w <- pi*seq(0,2,length=N)

# Colour scheme for the graphics:
redblue <- rgb( c(seq(1,0.1,length=30),rep(0,20)),
               rep(0,50),
               c(rep(0,20),seq(0.1,1,length=30)) )

# Carry out the construction of the random curves:
for (i in 1:50) {
  y[i,] <- ( A[i,1]*cos(i*w) + A[i,2]*sin(i*w))/sqrt(i)
  x <- x + y[i,]
}

# Produce the graphics to show the results
plot(x,type="l",lwd=3,ylim=c(-30,10))
text(500,-7.5,"= sum of ...",cex=1.75,font=2)
for (i in 1:50) lines(y[i,] -10 -i/2.5,col=redblue[i])

attr(x,'history') <- 'decomposeFT()'
attr(x,'desription') <- 'random series'
attr(x,'Fourier amplitudes') <- A
invisible(x)
}
