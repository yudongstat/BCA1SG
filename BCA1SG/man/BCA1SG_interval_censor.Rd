\name{BCA1SG_interval_censor}
\alias{BCA1SG_interval_censor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BCA1SG algorithm for interval-censored survival data
}
\description{
This function implements the BCA1SG algorithm on the semiparametric proportional hazard model for interval-censored data to solve the ML estimates of the model parameters.
}
\usage{
BCA1SG_interval_censor(input_data, initial_beta, initial_Lambda = function(x){x},
threshold = 1e-05, max_iter = 5000, max_stepsize = 10000, xi = 0.3, contraction = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_data}{
An object of class data.frame.
The structure of the data frame must be
\{lower bound of the survival time,upper bound of the survival time,covariate_1,...,covariate_p\}.
This data frame cannot contain missing values. See the dataset "duser" for an example.
}
  \item{initial_beta}{
The initial value of the regression coefficients.
The dimension of this input should comply with the dimension of the covariates.
}
  \item{initial_Lambda}{
An R function which serves as the initial value of the baseline cumulative hazard function.
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
Details about the BCA1SG algorithm can be found in Wang et al. (2020), and the details concerning the semiparametric proportional hazard model can be found in Section 4 of Huang and Wellner (1997).
}
\value{
  \item{distinct_time}{The set of distinct observation time points.}
  \item{est_Lambda}{The estimated baseline cumulative hazard function at the set of distinct observation time points.}
  \item{est_beta}{The estimated regression coefficients.}
  \item{iteration}{The number of iterations.}
  \item{timecost}{The computational time in seconds.}
}
\references{
Wang Y., Ye, Z.-S., and Cao, H.(2020). On Computation of Semi-Parametric Maximum Likelihood Estimators with Shape Constraints. Submitted.

Huang J. and Wellner, J.A.(1997). Interval-Censored Survival Data: A Review of Recent Progress. Proceedings of the Fifth Seattle Symposium in Biostatistics, 123-169.
}
\author{
Wang Y., Ye Z., and Cao H.
}
\note{
If we directly run this function on the data set "duser", we may get a different result from that presented in Section 6.2 of Wang et al. (2020). This is because the settings about the initial values of the nonparametric baseline cumulative hazard function are different.
}


\examples{
data(adapt_duser)
res<-BCA1SG_interval_censor(adapt_duser,initial_beta = rep(0,2))
res$est_beta
res$iteration
res$timecost
plot(res$distinct_time,res$est_Lambda,type="s",lwd=3,
xlab="time",ylab="Baseline cumulative hazard function")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Interval-censored data}% use one of  RShowDoc("KEYWORDS")
