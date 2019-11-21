#' A data set adapted from the data set "skiTum" in the package "spef"
#'
#' A data frame containing the recurrence of skin tumors. See Chiou et al.(2017) for details.
#'
#' @format A data frame with columns:
#' \describe{
#'  \item{id}{patient id (repeated for each recurrence).}
#'  \item{time}{observation time.}
#'  \item{count}{cumulative number of tumors till the current observation time.}
#'  \item{age}{patient's age at enrollment; age = 1 if greater than 65, age = 0 otherwise.}
#'  \item{male}{gender; male = 1, female = 0.}
#'  \item{dfmo}{treatment (DFMO) group = 1; placebo = 0.}
#'  \item{priorTumor}{number of prior tumor from diagnosis to randomization.}
#' }
"adapt_skiTum"



#' A data set adapted from the data set "duser" in the package "FHtest"
#'
#' Data set of 763 drug users in Badalona (Spain). The data come from the detoxification unit of Hospital Universitari Germans Trias i Pujol in Badalona, Spain. See Gomez et al.(2000) for details.
#'
#' @format A data frame with columns:
#' \describe{
#'  \item{left}{left endpoint of time to HIV-infection.}
#'  \item{right}{right endpoint of time to HIV-infectio.n}
#'  \item{zgen}{gender (0: male; 1: female).}
#'  \item{age}{patient's age.}
#' }
"adapt_duser"
