\name{BCA1SG_degradation}
\alias{BCA1SG_degradation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BCA1SG algorithm for degradation data
}
\description{
This function implements the BCA1SG algorithm on the semiparametric random-effects inverse Gaussian process model for degradation data to solve the ML estimates of the model parameters.
}
\usage{
BCA1SG_degradation(input_data, initial_delta, initial_r, initial_Lambda = function(x){x},
threshold = 1e-05, max_iter = 5000, max_stepsize = 1e+05, xi = 0.3, contraction = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_data}{
An object of class data.frame.
The structure of the data frame must be
\{subject ID, time of measurement, measurement\}..
This data frame cannot contain missing values. See the dataset "liner" for an example.
}
  \item{initial_delta}{
The initial value of the shape parameter of the gamma distributed scale parameter in the random-effects inverse Gaussian process. See Wang and Xu (2010) for details.
}
  \item{initial_r}{
The initial value of the rate parameter of the gamma distributed scale parameter in the random-effects inverse Gaussian process. See Wang and Xu (2010) for details.
}
  \item{initial_Lambda}{
An R function which serves as the initial value of the baseline mean function.
}
 \item{threshold}{
Convergence threshold. The algorithm is terminated when the infinity norm of the difference between successive iterates is less than the convergence threshold.
}
  \item{max_iter}{
Maximum number of iterations allowed.
}
  \item{max_stepsize}{
Maximum stepsize allowed.
}
  \item{xi}{
The xi parameter in the inexact backtracking line search algorithm. See Wang et al. (2020) for details.
}
  \item{contraction}{
The contraction parameter in the inexact backtracking line search algorithm. See Wang et al. (2020) for details.
}
}
\details{
Details about the BCA1SG algorithm can be found in Wang et al. (2020), and the details concerning the semiparametric random-effects inverse Gaussian process model can be found in Section 3 of Wang and Xu (2010).
}
\value{
  \item{distinct_time}{The set of distinct observation time points.}
  \item{est_Lambda}{The estimated baseline mean function at the set of distinct observation time points.}
  \item{est_delta}{The estimated shape parameter of the gamma distributed scale parameter in the random-effects inverse Gaussian process.}
  \item{est_r}{The estimated rate parameter of the gamma distributed scale parameter in the random-effects inverse Gaussian process.}
  \item{iteration}{The number of iterations.}
  \item{timecost}{The computational time in seconds.}
}
\references{
Wang Y., Ye, Z.-S., and Cao, H.(2020). On Computation of Semi-Parametric Maximum Likelihood Estimators with Shape Constraints. Submitted.

Wang X. and Xu, D.(2010). An Inverse Gaussian Process Model for Degradation Data. Technometrics, 52(2), 188-197.
}
\author{
Wang Y., Ye Z., and Cao H.
}

\examples{
data(liner)
res<-BCA1SG_degradation(liner,initial_delta = 1,initial_r = 1)
res$est_delta
res$est_r
res$iteration
res$timecost
plot(res$distinct_time,res$est_Lambda,type="s",lwd=3,xlab="time",ylab="Baseline mean function")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Degradation data}% use one of  RShowDoc("KEYWORDS")
