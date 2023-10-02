\name{simulate.asreml}
\alias{simulate.asreml}
\title{Produce sets of simulated data from a multivariate normal distribtion and save 
quantites related to the simulated data}
\description{Produce in parallel sets of simulated data corresponding to an \code{asreml} model,
             along with its fitted values and residuals. A variance matrix \code{V},
             corresponding to the \code{random} and \code{residual} models must be 
             supplied. What to save is specified by the \code{which} argument.}
\usage{
\method{simulate}{asreml}(object, nsim=100, seed = NULL, means=NULL, V, tolerance = 1E-10, 
         update = TRUE, trace = FALSE, which="data", units = "ignore", 
         ncores = 2, ...)}
\arguments{
 \item{object}{An \code{asreml} object from a call to \code{asreml} in which the 
                \code{data} argument has been set.}
 \item{means}{The \code{vector} of means to be used in generating simulated data sets. If
             it is \code{NULL}, the fitted values based on \code{object} are used. 
             It must be the same length as the response variable for \code{object}.}
 \item{V}{The fitted variance \code{matrix}, i.e. having the pattern and values that 
          conform to the model fit stored in the supplied \code{object}.}
 \item{nsim}{The number of data sets to be simulated.}
 \item{seed}{A single value, interpreted as an integer, that specifies the 
              starting value of the random number generator. The "L'Ecuyer-CMRG" random 
              generator is used and \code{nextRNGStream} is used to seed each core from the
              original \code{seed}.}
 \item{tolerance}{The value such that eigenvalues less than it are considered to be zero.}
 \item{update}{If \code{TRUE }then the arguments \code{R.param} 
              and \code{G.param} are set to those in the \code{asreml} object supplied 
              in \code{object} so that the values from the original model are used as 
              starting values. If \code{FALSE} then \code{asreml} \code{calls} are 
              evaluated, the only changes from the previous call being that 
              (i) the model is fitted to simulated data and (ii) modifications 
              specified via \code{\dots} are mode, except that changes cannot be made to 
              any of the models.}
 \item{trace}{If \code{TRUE} then partial iteration details are displayed when ASReml-R 
              functions are invoked; if \code{FALSE} then no output is displayed.}
 \item{which}{The quantites from the simulated data set to be stored. Any combination of
              \code{"response"}, \code{"residuals"} and \code{"fitted"}, or \code{"all"}. 
              If \code{residuals} and/or \code{fitted} is specified, those for the 
              analysis stored in \code{object} will be added to the \code{data.frame} 
              nominated in the \code{data} argument of \code{object} and the modified 
              \code{data.frame} added as a component named \code{data} in the \code{list} 
              that is the \code{value} returned by the \code{function}.}
 \item{units}{A \code{character} indicating whether the BLUPs for \code{units} are 
              added to the residuals when this reserved factor is included in the 
              \code{random} model. Possible values are \code{addtoresiduals} and 
              \code{ignore}.}
 \item{ncores}{A \code{numeric} specifying the number of cores to use in doing the 
              simulations. In choosing a value for \code{ncores}, it is necessary to 
              take into account other processes that are using parallel processing at 
              the same time.}
 \item{\dots}{Other arguments that are passed down to the function asreml. Changes to the 
            models are not allowed. Other changes are dangerous and generally should be 
            avoided.}
}
 \value{A \code{list} with the following components whose presence depends on the setting of 
       \code{which}:
\enumerate{
\item{\bold{observed:} present if \code{which} includes \code{residuals} or \code{fitted}, 
                       in which case it will be the \code{data.frame} on which the fit in 
                       \code{object} is based, with \code{residuals} and/or \code{fitted}.} 
\item{\bold{data:} present if \code{which} includes \code{data}, a \code{data.frame} 
                  containing the simulated data sets.}
\item{\bold{fitted:} present if \code{which} includes \code{fitted}, a \code{data.frame} 
                  containing the fitted values from the analyses of the simulated data sets.}
\item{\bold{residuals:} present if \code{which} includes \code{residuals}, a \code{data.frame} 
                  containing the residuals from the analyses of the simulated data sets.}}
}
\details{Generate \code{nsim} sets of data and analyse them using \code{asreml} 
         using the model in \code{object}, performing the generation and analysis of several 
         sets in parallel. Note, if the analysis for a data set does 
         not converge in \code{maxiter} iterations, it is discarded and 
         a replacement data set generated. The value of \code{maxiter} can be specified 
         in the call to \code{simulate.asreml}. The fitted values and residuals are extracted 
         as required. If \code{aom = TRUE} when the simulated data are analysed, standardised 
         conditional residuals are stored. If \code{which} includes \code{residuals} or
         \code{fitted}, the specified quantities for the observed data are added to the
         \code{data.frame} on which the fit in \code{object} is based.}
\author{Chris Brien}
\seealso{\code{asreml}, \code{\link{newfit.asreml}}, \code{\link{variofaces.asreml}}, \code{\link{plotVariofaces.data.frame}}, 
         \code{\link{set.seed}}.}
\examples{\dontrun{
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary.asrtests(current.asrt)
# Form variance matrix based on estimated variance parameters
s2 <- current.asr$sigma2
gamma.Row <- current.asr$gammas[1]
gamma.unit <- current.asr$gammas[2]
rho.r <- current.asr$gammas[4]
rho.c <- current.asr$gammas[5]
row.ar1 <- mat.ar1(order=10, rho=rho.r)
col.ar1 <- mat.ar1(order=15, rho=rho.c)
V <- gamma.Row * fac.sumop(Wheat.dat$Row) + 
     gamma.unit * diag(1, nrow=150, ncol=150) + 
     mat.dirprod(col.ar1, row.ar1)
V <- s2*V

#Produce residuals from 100 simulated data sets
resid <- simulate(current.asr, V=V, which="residuals", ncores = parallel::detectCores())
}}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
