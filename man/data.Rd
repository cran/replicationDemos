\name{svalbard}
\alias{AKRIM}
\alias{CRU}
\alias{crutemp}
\alias{F10.7cm}
\alias{forcings}
\alias{GISS.GCMs}
\alias{gistemp}
\alias{Lean1995}
\alias{Lean2004}
\alias{svalbard}
\alias{ssh2011.tab1}
\alias{vardo}
\alias{Douglasetal.tab1}
\alias{Douglasetal.tab2}
\alias{Mauna.Loa}
\alias{gdptemp03}
\alias{gisp2}
\alias{enso}
\alias{vostoc.co2}
\alias{vostoc.temp}
\alias{Be.10}
\alias{CMIP3.20c3m.sresa1b}
\alias{Scafetta2011.tab1}
\alias{fl1991}
\title{Data for demonstrations of replication and testing.}
\description{
Various data sets used in the demonstrations. Several of these are
'standard' data sets (CRU, Lean2004, AKRIM, crutemp, F10.7cm, forcings,
gistemp, Lean1995, GISP2, Mauna.Loa). Some are from tables in papers
(tab1, Douglasetal.tab1, Douglasetal.tab2,Scafetta2011.tab1).

The tables were copied digitally from the PDF-version in acroreader
(copy text) and then saves as ASCII-files, read in R, and then re-saved as
rda-files. The negative signs ('-') had to be set to '-' since the ASCII
code for the signs in the tables did not correspond to the ASCII code
used by R. Once these minor issues were fixed, these should be exact
reproductions of the tables in the papers.  

\code{ssh2011.tab1} is the data from Table 1 in Solheim et al. (2011)
\code{Douglasetal.tab1} and \code{Douglasetal.tab1} are from Douglas et
al.

The other data sets have been taken from the same sources as stated in
the papers. The URL from where these were obtained are given in the data
attributes (e.g. type \code{names(attributes(gisp2))}).

By copying the numbers in published tables, and providing these together
with the source code, we hopefully should be able to prove the inappropriateness
in a way that even deniers find it hard to deny. This is how science
should work - solid piece of work will mean that the results should be
reproduced over and over again...

All the data have the attribute 'url' which provides information about
the source from which the data were obtained.

\code{FL1991} reads the solar cycle lengths from the DMI web site:
Friis-Christensen, E., and K. Lassen, Length of the Solar Cycle: An
indicator of Solar Activity Closely Associated with Climate, Science
254, 698-700, (1991) 
}
\usage{
data(gisp2)
data(ssh2011.tab1)
data(svalbard)
data(vardo)
data(CRU)
data(Lean2004)
data(AKRIM)
data(crutemp)
data(F10.7cm)
data(forcings)
data(GISS.GCMs)
data(gistemp)
data(Lean1995)
data(Douglasetal.tab1)
data(Douglasetal.tab2)
data(Mauna.Loa)
data(enso)
data(gdptemp03)
data(vostoc.co2)
data(vostoc.temp)
data(Be.10)
data(CMIP3.20c3m.sresa1b)
data(Scafetta2011.tab1)
data(fl1991)
}

\value{Lists, data.frames, or vecors containing the relevant data.}
\author{R.E. Benestad}
\keyword{data}


