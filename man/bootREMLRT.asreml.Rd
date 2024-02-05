\name{bootREMLRT.asreml}
\alias{bootREMLRT.asreml}
\alias{bootREMLRT}
\title{Uses the parametric bootstrap to calculate the p-value for a REML ratio test to compare two models.}
\description{Extracts the REML log likelihood for two \code{asreml} objects 
  and forms the observed REML ratio statistic. It assumes that the second 
  \code{asreml} object is the result of fitting a model that is a reduced 
  version of the model for the first object and is considered to the null model. 
  Using the \code{mean} and \code{V}, \code{nboot} bootstrap samples of simulated 
  response values are generated in parallel; that is, \code{ncores} cores are used 
  and each is used to generate and analyse a sample.  The full and reduced models 
  are fitted to the data and if either analysis fails to converge another sample 
  is generated and analysed using the current core, with a maximum of 
  \code{max.retries} attempts to obtain a sample that converges for both analysis. 
  Thus the maximum number of data sets that will be generated is 
  \code{nboot} * \code{max.retries}. If a bootstrap sample converges for both analyses, 
  the REML ratio test statistic is formed for it. The p-value is then calculated 
  as \eqn{(k + 1) / (b + 1)}{(k + 1) / (b + 1)} where \eqn{k}{k} is the number 
  of simulated ratio test statistics greater than the observed test statistic and 
  \eqn{s}{b} is the number of bootstrap samples that were returned.
  
  The function checks that the models do not differ in either their fixed 
  or sparse models. It also check the difference in the number of variance 
  parameters between the two fits to the models, taking into account the 
  \code{bound.exclusions}.}
  
\usage{
\method{bootREMLRT}{asreml}(h0.asreml.obj, h1.asreml.obj, 
           nboot = 100, max.retries = 5, seed = NULL, 
           means=NULL, V = NULL, extra.matrix = NULL, ignore.terms = NULL, 
           fixed.spline.terms = NULL, 
           bound.exclusions = c("F","B","S","C"), 
           tolerance = 1E-10, update = TRUE, trace = FALSE, 
           ncores = 2, ...)
}
\arguments{
 \item{h0.asreml.obj}{\code{asreml} object containing the fit under the model for the 
             null hypothesis.}
 \item{h1.asreml.obj}{\code{asreml} object containing the fit under the model for the 
             alternative hypothesis.}
 \item{nboot}{The number of bootstrap samples to be generated.}
 \item{max.retries}{The maximum number of attempts to generate a sample whose analyses 
              converge for both models.}
 \item{seed}{A single value, interpreted as an integer, that specifies the 
              starting value of the random number generator. The "L'Ecuyer-CMRG" random 
              generator is used and \code{nextRNGStream} is used to seed each core from the
              original \code{seed}.}
 \item{means}{The \code{vector} of means to be used in generating simulated bootstrap 
             samples. If it is \code{NULL}, the fitted values based on \code{object} 
             are used. It must be the same length as the response variable for 
             \code{object}.}
 \item{V}{The fitted variance \code{matrix}, i.e. having the pattern and values that 
          conform to the model fit stored in the supplied \code{object}. If it is 
          \code{NULL}, \code{\link{estimateV.asreml}} is used to estimate the 
          variance matrix for the observations from the variance parameter estimates 
          from the \code{reduced.asreml.obj}.}
 \item{extra.matrix}{A \code{matrix} of order equal to the number of observations that is to 
               be added to the variance \code{matrix}, the latter based 
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
 \item{fixed.spline.terms}{A \code{character} vector giving one or more spline terms in the 
               \code{random} model that are regarded as fixed  and so are to be ignored 
               because they are not regarded as contributing to the variance. The term 
               names are those given in the \code{vparameters} component of the 
               \code{asreml} object or the \code{varcomp} component produced by 
               \code{summary.asreml}, but only up to the first exclamation mark 
               (\code{!}).}
 \item{bound.exclusions}{A \code{character} specifying one or more bound codes that 
                    will result in a variance parameter in the \code{random} model being 
                    excluded from contributing to the variance. If set to \code{NULL} 
                    then none will be excluded.}
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
 \item{ncores}{A \code{numeric} specifying the number of cores to use in doing the 
              simulations. In choosing a value for \code{ncores}, it is necessary to 
              take into account other processes that are using parallel processing at 
              the same time.}
 \item{\dots}{Other arguments that are passed down to the function asreml. Changes to the 
            models are not allowed. Other changes are dangerous and generally should be 
            avoided.}
}
\value{A \code{list} with the following components:
\enumerate{
\item{\bold{REMLRT:} the observed REML ratio statistic.} 
\item{\bold{p:} the bootstrap p-value for the observed test statistic.}
\item{\bold{DF:} the calculated difference in DF for the variance parameters in the two models.}
\item{\bold{totalunconverged:} the total number of unconverged analyses over the simulations.}
\item{\bold{REMLRT.sim:} a \code{numeric} containing the values of the ratio statistics 
                         for the simulated data. It has an attribute called \code{na.action} that 
                         can be retrieved using \code{attr(REMLRT.sim, which = "na.action")}; it 
                         contains a list of the simulation numbers that were abandoned because 
                         \code{max.retries} failed to converge for both models.}
\item{\bold{nunconverged:} the number of unconverged analyses for each bootstrap sample, the 
                  maximum being \code{max.retries}.}}
}

\note{A bootstrap sample is generated using a multivariate normal distribution with expected value as 
      specified by \code{means} and variance matrix given by \code{V}. Each simulated sample 
      is analysed according to the reduced model and, provided this analysis converges, 
      according to the full.model. If one of these analyses fails to converge, it is 
      abandoned and another sample is generated for this simulation. As many as 
      \code{max.retries} attempts are made to generate a data set for which both analyses 
      converge. If data set that converges for both analyses is not generated for a 
      simulation, \code{NA} is returned for that bootstrap sample. Hence, the maximum number 
      of data sets that will be generated is \code{nboot} * \code{max.retries} and less than 
      \code{nboot} samples will be generated if a data set that converges for both analyses 
      is not obtained within \code{max.retries} attempts. 

      If a bootstrap sample converges for both analyses, the REML ratio test statistic is calculated 
      as \eqn{2(log(REML)_F - log(REML)_R)}{2(log(REML)_F - log(REML)_R)}.

      The \code{DF} is calculated from the information in \code{full.asreml.obj} and 
      \code{reduced.asreml.obj}. The degrees of freedom are computed as the difference 
      between the two models in the number of variance parameters whose estimates do 
      not have a code for \code{bound} specified in \code{bound.exclusions}.

      If ASReml-R version 4 is being used then the codes specified in 
      \code{bound.exclusions} are not restricted to a subset of the default codes, but 
      a warning is issued if a code other than these is specified. 
      For ASReml-R version 3, only a subset of the default codes are allowed: 
      \code{F} (\code{Fixed}), \code{B} (\code{Boundary}), \code{C} (\code{Constrained}) 
      and \code{S} (\code{Singular}). 
}

\author{Chris Brien}
\seealso{\code{\link{REMLRT.asreml}}, \code{\link{infoCriteria.asreml}}, \code{\link{newfit.asreml}}, 
         \code{\link{testranfix.asrtests}}}
\examples{\dontrun{
    bootREMLRT(ICV.max, ICV.red, ncores = parallel::detectCores())
}}
\keyword{htest}
\keyword{asreml}