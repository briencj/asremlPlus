\name{addto.test.summary}
\alias{addto.test.summary}
\title{Adds a row to a \code{test.summary} \code{data.frame}.}
\description{A row that summarizes the result of a proposed change to a model is added to a 
             \code{test.summary} \code{data.frame}. Only the values of those arguments for whch there 
             are columns in \code{test.summary} will be included in the row.}

\usage{addto.test.summary(test.summary, terms, DF = 1, denDF = NA, 
                          p = NA, AIC = NA, BIC = NA,
                          action = "Boundary")}
\arguments{
 \item{test.summary}{A \code{\link{data.frame}} whose columns are a subset of \code{terms}, 
                     \code{DF}, \code{denDF}, \code{p}, \code{AIC}, \code{BIC} and \code{action}. 
                     Each row summarizes the results of proposed changes to the fitted model. 
                     See \code{\link{asrtests.object}} for more information.}
 \item{terms}{A \code{\link{character}} giving the name of a term that might be added to or removed 
                     from the model or a label indicating a change that might be made to the model.}
 \item{DF}{A \code{\link{numeric}} giving the numerator degrees of freedom for a Wald F-statistic 
                     or the number of variance parameters in the current model minus the number 
                     in the proposed model.}
 \item{denDF}{A \code{\link{numeric}} giving the denominator degrees of freedom for a 
                     Wald F-statistic.}
 \item{p}{A \code{\link{numeric}} giving the p-value for a Wald F-statistic or REML ratio test.}
 \item{AIC}{A \code{\link{numeric}} giving Akiake Information Criterion (AIC) for a model or the 
                     difference between the AIC values for the current and proposed models.}
 \item{BIC}{A \code{\link{numeric}} giving Bayesian (Schwarz) Information Criterion for a model or the 
                     difference between the AIC values for the current and proposed models.}
 \item{action}{A \code{\link{character}} giving what action was taken with respect to the 
                     proposed change. See \code{\link{asrtests.object}} for more 
                     information.}
}
\value{A \code{\link{data.frame}}.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{asrtests.object}},  \code{\link{print.test.summary}}}
\examples{
  \dontrun{
data(Wheat.dat)

## Fit an autocorrelation model
ar1.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                  random = ~ Row + Column + units,
                  residual = ~ ar1(Row):ar1(Column), 
                  data=Wheat.dat)
ar1.asrt <- as.asrtests(ar1.asr, NULL, NULL, 
                        label = "Autocorrelation model")
ar1.asrt <- rmboundary.asrtests(ar1.asrt)

## Fit a tensor spline
Wheat.dat <- within(Wheat.dat, 
                    {
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn - mean(unique(cColumn))
                    })
ts.asr <- asreml(yield ~ Rep + cRow + cColumn + WithinColPairs + 
                         Variety, 
                  random = ~ spl(cRow) + spl(cColumn) + 
                             dev(cRow) + dev(cColumn) + 
                             spl(cRow):cColumn + cRow:spl(cColumn) + 
                             spl(cRow):spl(cColumn),
                  residual = ~ Row:Column, 
                  data=Wheat.dat)
ts.asrt <- as.asrtests(ts.asr, NULL, NULL, 
                       label = "Tensor spline model")
ts.asrt <- rmboundary.asrtests(ts.asrt)
ar1.ic <- infoCriteria(ar1.asrt$asreml.obj)
ts.ic <- infoCriteria(ts.asrt$asreml.obj)
if (ar1.ic$AIC < ts.ic$AIC)
{
  ic.diff <- ar1.ic - ts.ic
  new.asrt <- ar1.asrt 
  new.asrt$test.summary <- addto.test.summary(ar1.asrt$test.summary, 
                                              terms = "Compare ar1 to ts", 
                                              DF = ic.diff$varDF, 
                                              AIC = ic.diff$AIC, BIC = ic.diff$BIC, 
                                              action = "Chose ar1")
} else
{
  ic.diff <- ts.ic - ar1.ic
  new.asrt <- ts.asrt
  new.asrt$test.summary <- addto.test.summary(ts.asrt$test.summary, 
                                              terms = "Compare ar1 to ts", 
                                              DF = ic.diff$varDF, 
                                              AIC = ic.diff$AIC, BIC = ic.diff$BIC, 
                                              action = "Chose ts")
}
}}

\keyword{asreml}
\keyword{htest}