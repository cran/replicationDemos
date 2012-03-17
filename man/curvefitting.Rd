\name{Replication demo scripts.}
\alias{forskning.no}
\alias{Humlum.et.al.2011}
\alias{decomposeFT}
\alias{LoehleScafetta2011}
\alias{Scafetta2010}
\alias{Scafetta2011}
\alias{hfit}
\alias{ENSO.example}
\alias{resonance}
\alias{Scafetta.tab1}
\alias{resonanceTest}
\alias{Solheim.et.al.2011}
\alias{DJF}
\alias{do.vardo}
\alias{Scafetta2006}
\alias{check.table1}
\alias{obs2tab1}
\alias{Walker.test}
\alias{MM2004}
\alias{Douglas2007}
\alias{paleaoproxy}
\alias{FL1991}
\title{Demonstrations which debunk some methods and analytical set-ups.}
\description{
This set of function provide demonstrations showing why the choice of
methods and analytrical set-ups in a number of papers are
inappropriate. This computer code is also used to do the analysis on
which the paper 'Bad science papers' by Benestad et al. is based.

Some of these analyses are indeed very trivial, but carried out
nevertheless. The objective with this R-package is partly to show that
code sharing and open-source can be an effective means of resolving
differences. Methods and analytical set-up should be tested with
surrogate data for which the anwers are known a priori (a kind of method
calibration and evaluation). The spirit of this is very much like the
replication carried out by Benestad and Schmidt (2009) \url{http://pubs.giss.nasa.gov/abs/be02100q.html}

\code{Humlum.et.al.2011()} exmines the results of Humlum et al. (2011)
\url{http://www.sciencedirect.com/science/article/pii/S0921818111001457}
and extends their analysis to the data that they cut off.

\code{LoehleScafetta2011()}. Reference: Loehle and Scafetta (2011) \url{http://benthamscience.com/open/toascj/articles/V005/74TOASCJ.htm}

\code{Scafetta2010()} Some replication and some demonstrations showing
how different tempertures change results, and how different noise
processes may seem to contain long-term cycles. 

\code{ENSO.example()} shows how the strategy adopbed in Scafetta (2011)
fails when applied to ENSO (the NINO3.4 index). \code{resonance()} shows
that a system with resonance will pick up the resonant frequency from
any noisy forcing - this is analogous to the whine from the wind blowing
aroundcorners, and how musical pipes/trumpets work (the function
\code{resonanceTest} provides some
demonstrations). \code{Scafetta.tab1()} replicates Table 1 and the
harmonic regressions using \code{hfit} fits harmonics according to
eq. (3)in Scafetta (2011). \code{decomposeFT()} 
shows how any curve can be represented as a sum of harmonics - Fourier
series. REference: Scafetta,N.,Testing an astronomically based
decadal-scale empirical harmonic climate model versus the
IPCC (2007) general...., Journal of Atmospheric and Solar-Terrestrial
Physics (2011) doi:10.1016/j.jastp.2011.12.005.
Also on \url{http://arxiv.org/abs/1201.1301} 

\code{Solheim.et.al.2011()} replicates the analysis carried out by
solheim et al. (2011), and shows how their bootstrapping is biased, and
that their results are not really significant after
all. \code{Walker.test()} is the function for the Walker test to test
the significance when many tests are made (problem of multiplicity/field
significance). \code{check.table1()} and \code{obs2tab1()} are
supporting functions for this replication. \code{do.vardo} repeats the
analysis for Vardo. Reference: Solheim et al (2011)
\url{http://arxiv.org/pdf/1112.3256}. +code{DJF} extracts the
December-February mean.

\code{Scafetta2006()} is a modification of the R-scripts used in
Benestad & Schmidt (2008) replicating and testing the analysis in the
papers Scafetta (2006).

\code{MM2004()} The R-script used to carry out the analysis of Benestad
(2004) \url{http://www.int-res.com/articles/cr2004/27/c027p171.pdf}

\code{Douglas2007()} shows how the conficende intervals used in Douglas
et al. (2007) excludes many of the GCM results, on which it was
based. Reference: 
Douglas et al. (2007) \url{http://www.pas.rochester.edu/~douglass/papers/Published\%20JOC1651.pdf}

\code{paleaoproxy()} examines some assertions made by Veizer (2005):
Comparison between cosmogenic Be-10 isotope proxies, CO2, and
temperature from the Vostoc ice core. \url{http://www.realclimate.org/index.php/archives/2005/05/on-veizers-celestial-climate-driver/}

\code{forskning.no} replicates a figure (again, by Solheim et al.) shown in
\url{http://www.forskning.no/artikler/2012/mars/316178} that claims that
the warming has stopped. The figure gives a wrong impression, as it's
shown on monthly scales where the monthly anomalies swamp the signal. By
plottingthe monthly mean values instead, the warming trends becomes
visible. 
 
These papers are also discussed on \url{http://www.RealClimate.org}: \url{http://www.realclimate.org/index.php/archives/2011/12/curve-fitting-and-natural-cycles-the-best-part/} and \url{http://www.realclimate.org/index.php/archives/2011/08/an-exercise-about-meaningful-numbers-examples-from-celestial-attribution-studies/}
}

\usage{
Humlum.et.al.2011()
forskning.no(uah=TRUE,rss=TRUE,giss=TRUE,ncdc=TRUE,hadcrut3=TRUE,
             lwd=1,ylim=c(-0.4,1),xlim=c(1995,2012),
             base.period=1981:2010,type="s")
decomposeFT(N=1000)
LoehleScafetta2011()
Scafetta2010()
Scafetta2011()
ENSO.example(interval=1980:1989)
resonance(x0=0,v0=0,h=0.1,main=NULL,
          F.ext=rnorm(5000) + 0.1*cos(2*pi*(1:5000)/500),f=0.1,m=0.1,w0=1)
Scafetta.tab1()
resonanceTest()
Walker.test(N,alpha=0.05)
check.table1(ssh2011.tab1=NULL)
Solheim.et.al.2011(obs=NULL,ssh2011.tab1=NULL,N.tests=30000)
DJF(obs,ssh2011.tab1=NULL)
do.vardo()
Scafetta2006(GISS.temp = TRUE, do.MonteCarlo = TRUE, test.bp = TRUE, 
    lag = 0, stepwise = TRUE, interval = 1958:2000, same.interval = TRUE, 
    all.data = FALSE, SW06.coefs.only = FALSE, wavelets.only = FALSE, 
    bivariate = TRUE, figures = TRUE, tables = TRUE, wavelet = TRUE, 
    boundary = "reflection")
MM2004()
Douglas2007()
paleaoproxy()
hfit(x,t,T1=60,T2=20)
FL1991(dmi=FALSE,url="http://web.dmi.dk/fsweb/solarterrestrial/sunclimate/SCL.txt")
}
\arguments{
    \item{obs}{observations: station data.}
    \item{GISS.temp}{TRUE: use GISTEM, otherwise hadCRUT3}
    \item{test.bp}{pand-pas test}
    \item{lag}{lag between forcing and temperature response (years)}
    \item{stepwise}{TRUE: step-wise regression}
    \item{interval}{Time interval to analyse}
    \item{same.interval}{TRUE: same.interval as SW2006a}
    \item{all.data}{TRUE: use all the data data}
    \item{SW06.coefs.only}{TRUE: only estimate the coefficients in SW2006}
    \item{wavelets.only}{TRUE: only do the wavelet}
    \item{bivariate}{TRUE: only include soalr forcing and GHG in regression}
    \item{figures}{TRUE: plot figures}
    \item{tables}{TRUE: print tables}
    \item{ssh2011.tab1}{Table as in \code{ssh2011.tab1} for defining
      solar cycle epochs, correponding to Table 1 in Solheim et al. (2011).}
    \item{wavelet}{TRUE: do wav}
    \item{boundary}{Boundary for the wavelet: mra}
    \item{x0}{Initial condition for x of oscillator}
    \item{v0}{Initial condition for v = dx/td}
    \item{h}{time step size}
    \item{main}{Main title in plot}
    \item{F.ext}{External forcing for oscillator}
    \item{f}{friction term for damped oscillator}
    \item{m}{inertia term for oscillator}
    \item{w0}{frequency term}
    \item{alpha}{level of statistical significance}
    \item{N.tests}{Number of Monte-Carlo simulations}
    \item{x}{curve to which harmonics are fitted}
    \item{t}{times in year}
    \item{T1}{periodicity for first harmonic}
    \item{T2}{periodicity for second harmonic}
    \item{do.MonteCarlo}{TRUE: Carry out Monte-Carlo simulations}
    \item{N}{length of time series}
    \item{url}{URL address of FL(1991) solar cycle length data.}
    \item{dmi}{Flag: if TRUE read from DMI web site}
    \item{uah}{Univ. Alabama Huntsville}
    \item{rss}{REmote Sensing Systems}
    \item{giss}{NASA/GISS}
    \item{ncdc}{NOAA NCDC}
    \item{hadcrut3}{Hadley Centre and Climate Research Unit, U.K.}
    \item{lwd}{see \code{\link{plot}}}
    \item{ylim}{see \code{\link{plot}}}
    \item{xlim}{see \code{\link{plot}}}
    \item{base.period}{Reference for anomalies.}
    \item{type}{see \code{\link{plot}}}

}
\value{A table or lists containing the relevant data.}
\author{R.E. Benestad}
\examples{
forskning.no()
Humlum.et.al.2011()
decomposeFT()
LoehleScafetta2011()
Walker.test(N=10)
check.table1()
Solheim.et.al.2011(N.tests=3000)
}
\keyword{manip}


