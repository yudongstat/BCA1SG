\name{adapt_duser}
\alias{adapt_duser}
\docType{data}
\title{
A data set adapted from the data set "duser" in the package "FHtest"
}
\description{
Data set of 763 drug users in Badalona (Spain). The data come from the detoxification unit of Hospital Universitari Germans Trias i Pujol in Badalona, Spain. See Gómez et al. (2000) for details.
}
\usage{data("adapt_duser")}
\format{
  A data frame with 763 observations on the following 4 variables.
  \describe{
    \item{\code{left}}{left endpoint of time to HIV-infection.}
    \item{\code{right}}{right endpoint of time to HIV-infectio.n}
    \item{\code{zgen}}{gender (0: male; 1: female).}
    \item{\code{age}}{patient's age.}
  }
}
\source{
The "FHtest" package.
}
\references{
Gómez, G., Calle, M. L., Egea, J. M., and Muga, R.(2000). Risk of HIV infection as a function of the duration of intravenous drug use: A non-parametric Bayesian approach. Statistics in Medicine 19(19), 2641–2656.
}
\keyword{datasets}
