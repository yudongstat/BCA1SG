\name{adapt_skiTum}
\alias{adapt_skiTum}
\docType{data}
\title{
A data set adapted from the data set "skiTum" in the package "spef"
}
\description{
A data frame containing the recurrence of skin tumors. See Chiou et al. (2017) for details.
}
\usage{data("adapt_skiTum")}
\format{
  A data frame with 894 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{patient id (repeated for each recurrence).}
    \item{\code{time}}{observation time.}
    \item{\code{count}}{cumulative number of tumors till the current observation time.}
    \item{\code{age}}{patient's age at enrollment; age = 1 if greater than 65, age = 0 otherwise.}
    \item{\code{male}}{gender; male = 1, female = 0.}
    \item{\code{dfmo}}{treatment (DFMO) group = 1; placebo = 0.}
    \item{\code{priorTumor}}{number of prior tumor from diagnosis to randomization.}
  }
}

\source{
The "spef" package.
}
\references{
Chiou S. H., Xu G., Yan J., and Huang C.(2018). Semiparametric Estimation of the Accelerated Mean Model with Panel Count Data under Informative Examination Times. Biometrics, 74(3), 944-953.
}
\keyword{datasets}
