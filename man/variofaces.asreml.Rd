\name{variofaces.asreml}
\alias{variofaces.asreml}
\alias{variofaces}
\title{Plots empirical variogram faces, including envelopes, as described by 
Stefanova, Smith & Cullis (2009).}
\description{A function that produces a plot for each face of an empirical 2D 
             \code{variogram} based on residuals produced after the fitting of a model 
             using the function \code{asreml}. 
             It also adds envelopes to the plot by simulating data sets in parallel 
             from a multivariate normal distribution with expectation equal to the 
             fitted values obtained from the fixed and spline terms and variance 
             matrix equal to the fitted variance matrix 
             (Stefanova, Smith & Cullis, 2009). The plot is constrolled by the 
             \code{residual} model, which must consist of two factors corresponding to 
             the two physical dimensions underlying the data. It can also have a third 
             term involving the \code{at} or \code{dsum} function that defines sections 
             of the data, such as experiments in different environments. 
             In this case, the two variogram faces are produced for each section.}
\usage{
\method{variofaces}{asreml}(asreml.obj, means=NULL, V=NULL, nsim=100, seed = NULL, 
           extra.matrix = NULL, ignore.terms = NULL, fixed.spline.terms = NULL, 
           bound.exclusions = c("F","B","S","C"), tolerance=1E-10, 
           units = "ignore", update = TRUE, trace = FALSE, 
           graphics.device=NULL, ncores = parallel::detectCores(), ...)}
\arguments{
 \item{asreml.obj}{An \code{asreml} object from a call to \code{asreml} in which the 
                \code{data} argument has been set.}
 \item{means}{The \code{vector} of means to be used in generating simulated data sets. If
             it is \code{NULL}, the fitted values based on \code{object} are used. 
             It must be the same length as the response variable for \code{object}.}
 \item{V}{The fitted variance \code{matrix}, i.e. having the appropriate pattern and values 
          given the model fitted to the observed data and the estimates of the 
          parameters obtained. If \code{V} is \code{NULL} then \code{\link{estimateV.asreml}} 
          is called to obtain it from \code{asreml.obj}}
 \item{nsim}{The number of data sets to be simulated in obtaining the envelopes.}
 \item{seed}{A single value, interpreted as an integer, that specifies the 
              starting value of the random number generator. The "L'Ecuyer-CMRG" random 
              generator is used and \code{nextRNGStream} is used to seed each core from the
              original \code{seed}.}
 \item{extra.matrix}{A \code{matrix} of order equal to the number of observations that is to 
               be added to the varaince \code{matrix}, the latter based 
               on the information in \code{asreml.obj}. It is assumed that the sigma-parameterized 
               values of the variance parameter estimates, such as is given in the \code{varcomp} 
               component of \code{summary.asreml}, have been used in calculating 
               \code{extra.matrix}; the values in the \code{vparameters} component of 
               \code{G.param} and \code{R.param} may be either gamma- or sigma-parameterized. 
               The argument \code{extra.matrix} can be used in conjunction with 
               \code{ignore.terms} as a workaround to include components of the variance matrix 
               for variance functions that have not been implemented in \code{estimateV}.}
 \item{ignore.terms}{A \code{character} giving terms from either the \code{random} or 
               \code{residual} models that are to be ignored in that their contributions to 
               the variance is not to be included in the estimated matrix. The term names are those 
               given in the \code{vparameters} component of the \code{asreml} object or the 
               \code{varcomp} component produced by \code{summary.asreml}, but only up to the 
               first exclamation mark (\code{!}). This can be used 
               in conjunction with \code{estimateV.asreml} as a workaround to include components 
               of the variance matrix for variance functions that have not been implemented                
               in \code{estimateV}.}
 \item{fixed.spline.terms}{A \code{character} vector giving one or mor spline terms in the 
               \code{random} model that are regarded as fixed  and so are to be ignored because 
               they are not regarded as contributing to the variance. The term names are those
               given in the \code{vparameters} component of the \code{asreml} object or the 
               \code{varcomp} component produced by \code{summary.asreml}, but only up to the
               first exclamation mark (\code{!}).}
 \item{bound.exclusions}{A \code{character} specifying one or more bound codes that 
                    will result in a variance parameter in the \code{random} model being excluded 
                    from contributing to the variance. If set to \code{NULL} then none will
                    be excluded.}
 \item{tolerance}{The value such that eigenvalues less than it are considered to be zero.}
 \item{units}{A \code{character} indicating whether the BLUPs for \code{units} are 
              added to the residuals when this reserved factor is included in the 
              \code{random} model. Possible values are \code{addtoresiduals} and 
              \code{ignore}. If standardized conditional residuals are plotted and 
              the BLUPs for \code{units} are to be added then it is the standardized 
              BLUPs that are added.}
 \item{update}{If \code{TRUE }then the arguments \code{R.param} 
              and \code{G.param} are set to those in the \code{asreml} object 
              supplied in \code{object} so that the values from the original 
              model are used as starting values. If \code{FALSE} then calls are made to 
              \code{asreml} in which the only changes from the previous call are 
              (i) the model is fitted to simulated data and (ii) modifications 
              specified via \code{\dots} are made, except that changes cannot be made to 
              any of the models.}
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{graphics.device}{A \code{character} specifying a graphics device for plotting. 
              The default is \cr\code{graphics.device = NULL}, which will result 
              in plots being produced on the current graphics device. Setting it to 
              \code{"windows"}, for example,  will result in a windows graphics 
              device being  opened.}
 \item{ncores}{A \code{numeric} specifying the number of cores to use in doing the 
             simulations.}
 \item{\dots}{Other arguments that are passed down to the function \code{asreml}. Changes 
              to the models are not allowed. Other changes are dangerous and generally 
              should be avoided.}
 }
\value{A \code{list} with the following components:
\enumerate{
\item{\bold{face1:} a \code{data.frame} containing the variogram values on 
                    which the plot for the first dimension is based.} 
\item{\bold{face2:} a \code{data.frame} containing the variogram values on 
                    which the plot for the second dimension is based.}}}
\details{The \code{residual} model is scanned to ensure that it involves only two factors 
         not included in the \code{at} function, and to see if it has a third factor in
         an \code{at} function. If so, the faces of the 2D variogram, each based on one 
         of the two non-\code{at} factors, are derived from the residuals in the 
         supplied \code{asreml} object using \code{asreml.variogram}, this yielding the observed 
         \code{variogram} faces. If \code{aom} was set to \code{TRUE} for the \code{asreml} 
         object, the standardized consitional residuals are used. 
         Then \code{nsim} data sets are generated by 
         adding the \code{fitted.values}, extracted from the \code{asreml} object,
         to a vector of values randomly generated from a normal distribution with 
         expectation zero and variance matrix \code{V}. Each data set is analyzed 
         using the model in \code{object} and several sets are generated and analyzed 
         in parallel. The variogram values for the faces are 
         obtained using \code{asreml.variogram} stored. Note, if the analysis for a 
         data set does not converge in \code{maxiter} iterations, it is discarded and 
         a replacement data set generated. The value of \code{maxiter} can be specified 
         in the call to \code{variofaces.asreml}. Plots are produced for each face and 
        include the observed values and the 2.5\%, 50\% & 97.5\% quantiles.}
\references{Stefanova, K. T., Smith, A. B. & Cullis, B. R. (2009) Enhanced diagnostics for the 
            spatial analysis of field trials. \emph{Journal of Agricultural, Biological, 
            and Environmental Statistics}, \bold{14}, 392--410.}
\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{asreml}, 
         \code{\link{plotVariofaces.data.frame}}, \code{\link{simulate.asreml}}, \code{\link{set.seed}}.}
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
variofaces(current.asr, V=V, ncores = 2)
}}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
