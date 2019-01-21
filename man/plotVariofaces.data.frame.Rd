\name{plotVariofaces.data.frame}
\alias{plotVariofaces.data.frame}
\alias{plotVariofaces}
\title{Plots empirical variogram faces, including envelopes, from supplied residuals as 
described by Stefanova, Smith & Cullis (2009).}
\description{Produces a plot for each face of an empirical 2D 
             \code{variogram} based on supplied residuals from both an observed data set 
             and simulated data sets. Those from simulated data sets are used to
             produce confidence envelopes If the data consists of sections, such as separate
             experiments, the two variogram faces are produced for each section. This 
             function is less efficient in storage terms than \code{\link{variofaces.asreml}},
             because here the residuals from all simulated data sets must be saved, in 
             addition to the values for the variogram faces; in 
             \code{\link{variofaces.asreml}}, the residuals for each simulated data set are 
             discarded after the variogram has been calculated. On the other hand, the 
             present function is more flexible, because there is no restriction on how the 
             residuals are obtained.}
\usage{\method{plotVariofaces}{data.frame}(data, residuals, restype="Residuals", ...)}
\arguments{
 \item{data}{A \code{data.frame} with either 3 or 4 columns. Only if there are 4 columns, the
             first should be a factor indexing sections for which separate variogram plots are 
             to be produced. In either case, the other 3 columns should be, in order, 
             (i) a factor indexing the x-direction, (ii) a factor indexing the y-direction, 
             and (iii) the residuals for the observed response.}
 \item{residuals}{A \code{data.frame}, with either 2 or 3 initial columns followed by 
                  columns, each of which are the residuals from a simulated 
                  data set.}
 \item{restype}{A \code{character} describing the type of residuals that have been supplied. 
                It will be used in the plot titles.}
 \item{\dots}{Other arguments that are passed down to the function \code{asreml.variogram}.}
 }
\value{A \code{list} with the following components:
\enumerate{
\item{\bold{face1:} a \code{data.frame} containing the variogram values on 
                    which the plot for the first dimension is based.} 
\item{\bold{face2:} a \code{data.frame} containing the variogram values on 
                    which the plot for the second dimension is based.}}}
\details{For each set of residuals, \code{asreml.variogram} is used to obtain the empirical 
         variogram, from which the values for its faces are obtained.  Plots are produced for
         each face and include the observed residuals and the 2.5\%, 50\% & 97.5\% quantiles.}
\references{Stefanova, K. T., Smith, A. B. & Cullis, B. R. (2009) Enhanced diagnostics for the 
            spatial analysis of field trials. \emph{Journal of Agricultural, Biological, 
            and Environmental Statistics}, \bold{14}, 392--410.}
\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{asreml}, \code{asreml.variogram},
        \code{\link{variofaces.asreml}}, \code{\link{simulate.asreml}}.}
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

#Produce variogram faces plot (Stefanaova et al, 2009)
resid <- simulate(current.asr, V=V, which="residuals")
resid$residuals <- cbind(resid$observed[c("Row","Column")],
                         resid$residuals)
plotVariofaces(data=resid$observed[c("Row","Column","residuals")],
               residuals=resid$residuals, 
               restype="Standardized conditional residuals")
}}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
