resonance <- function(x0=0,v0=0,h=0.1,main=NULL,
                      F.ext=rnorm(5000) + 0.1*cos(2*pi*(1:5000)/500),
                      f=0.1,m=0.1,w0=1) {
# It's a long time since I last used 4th order RungeKutta to simulate an ODE.  
# But I think I still remember the essence, and I used the information about
# forced damped damped oscillator from URL
# http://physics.tamuk.edu/~suson/html/4390/DiffEq.html

  require(deSolve)

  dvdt <- function(t,x,parms) {
    with(as.list(parms), {
      dv <- -f/m*v - w0^2*x + F.ext/m
      list(dv)
    } )

  }

  dxdt <- function(t,v,parms) {
    with(as.list(parms), {
      dxdt=v
      list(dxdt)
    } )
  }

# rk4(y, times, func, parms)

  n <- length(F.ext); x <- rep(NA,n); v <- x
  x[1] <- x0; v[1] <- v0; time <- c(0,h)

  for (i in 2:n) { 
    print(paste("i=",i,"x=",round(x[i-1],3),"v=",round(v[i-1],3)))
    X <- c(N = x[i-1]); V <- c(N= v[i-1])
    parmsv <- c(v=v[i-1],x=x[i-1],f=f,w0=w0,F.ext=F.ext[i],m=m)
    parmsx <- c(v=v[i-1])
    outv <- as.data.frame(rk4(V,time,dvdt,parmsv))
    outx <- as.data.frame(rk4(X,time,dxdt,parmsx))
    v[i] <- outv$N[2]; x[i] <- outx$N[2]
  }
  if(is.null(main)) main <- "Forced damped oscillator"
  plot(1:n,x,type="l",lwd=2,main=main,
       sub=paste("m=",m,"w0=",w0,"f=",f,"h=",h))
  lines(1:n,F.ext,col="red",lty=2)

  results <- list(x=x,v=v,F.ext=F.ext)
  invisible(results)
}

resonanceTest <- function() {
  dev.new()
  resonance(x0=1,v0=0,h=0.01,main="Unforced oscillator",
          F.ext=rep(0.00,1000),f=0.00,m=1,w0=1)
  
  dev.new()
  resonance(x0=1,v0=0,h=0.1,main="Unforced damped oscillator",
          F.ext=rep(0,1000),f=0.03,m=0.1,w0=1)

  dev.new()
  resonance(x0=1,v0=0,h=0.1,main="Damped oscillator forced w. sinusoid",
          F.ext=sin(2*pi*seq(0,1000,by=1)/500),f=0.3,m=0.1,w0=1)

  dev.new()          
  F.ext <- decomposeFT(N=10000)
  
  dev.new()
  resonance(x0=1,v0=0,h=0.1,main="Damped oscillator forced w. random signal",
          F.ext=F.ext,f=0.3,m=0.1,w0=1)

  dev.new()
  y <- resonance()

  dev.new()
  spectrum(y)

  data(co2,envir=environment())
  F.ext <- log(co2)
  spectrum(F.ext,main="log|CO2|")
  dev.new()
  
  resonance(x0=0,v0=0,h=0.1,main="Damped oscillator forced w. ln|CO2|",
          F.ext=F.ext,f=3,m=10,w0=1)
}
