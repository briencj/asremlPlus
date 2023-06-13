\name{addSpatialModel.asrtests}
\alias{addSpatialModel.asrtests}
\alias{addSpatialModel}
\title{Adds, to a supplied model, a spatial model that accounts for local spatial variation.}
\description{Adds either a correlation, two-dimensional tensor-product natural cubic 
       smoothing spline (TPNCSS), or a two-dimensional tensor-product penalized P-spline 
       model (TPPS) to account for the local spatial variation exhibited by a response variable 
       measured on a potentially irregular grid of rows and columns of the units. The data may 
       be arranged in sections, for each of which there is a grid and for which the model is to 
       be fitted separately. Also, the rows and columns of a grid are not necessarily one 
       observational unit wide. 
       
       No hypothesis testing or comparison of information criteria is made. To only change 
       the terms based on a comparison  of information criteria use 
       \code{\link{chooseSpatialModelOnIC.asrtests}}.
       
       The model fit supplied in the \code{asrtests.obj} should not include terms that will 
       be included in the local spatial model. All spatial model terms are fitted as fixed or 
       random. Consequently, the residual model does not have to be iid.

       One or more rows is added for each \code{section} to the \code{test.summary} 
       \code{data.frame}. Convergence and the occurrence of fixed correlations in fitting the 
       model is checked and a note included in the \code{action} if there was not. 
       All components of the \code{\link{asrtests.object}} are updated for the new model.
}
\usage{
\method{addSpatialModel}{asrtests}(asrtests.obj, spatial.model = "TPPS", 
                sections = NULL, 
                row.covar = "cRow", col.covar = "cCol", 
                row.factor = "Row", col.factor = "Col", 
                corr.funcs = c("ar1", "ar1"), 
                row.corrFitfirst = TRUE, 
                dropRowterm = NULL, dropColterm = NULL, 
                nsegs = NULL, nestorder = c(1,1), 
                degree = c(3,3), difforder = c(2,2),
                asreml.option = "mbf",  tpps4mbf.obj = NULL, 
                allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                checkboundaryonly = FALSE, update = FALSE, 
                IClikelihood = "full", ...)
}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components 
       (i) \code{asreml.obj}, (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{spatial.model}{A single \code{character} string nominating the type of spatial 
               model to fit. Possible values are \code{corr}, \code{TPNCSS} and 
               \code{TPPS}.}
 \item{sections}{A single \code{character} string that species the name of the column 
               in the \code{\link{data.frame}} that contains the \code{\link{factor}} 
               that identifies different sections of the data 
               to which separate spatial models are to be fitted.}
 \item{row.covar}{A single \code{character} string nominating a \code{\link{numeric}} 
               that contains the values of a centred covariate indexing the rows of a 
               grid. The \code{\link{numeric}} must be a column in the 
               \code{\link{data.frame}} stored in the \code{asreml.obj} that is a 
               component of the supplied \code{asrtests.obj}.}
 \item{col.covar}{A single \code{character} string nominating a \code{\link{numeric}} 
               that contains the values of a centred covariate indexing the columns 
               of a grid. The \code{\link{numeric}} must be a column in the 
               \code{\link{data.frame}} stored in the \code{asreml.obj} that is a 
               component of the supplied \code{asrtests.obj}.}
 \item{row.factor}{A single \code{character} string nominating a \code{\link{factor}} 
               that indexes the rows of a grid that are to be one dimension of a 
               spatial correlation model. The \code{\link{factor}} must a column in 
               the \code{\link{data.frame}} stored in the \code{asreml.obj} that is 
               a component of the supplied \code{asrtests.obj}.}
 \item{col.factor}{A single \code{character} string nominating a \code{\link{factor}} 
               that indexes the columns of a grid that are to be one dimension of a 
               spatial correlation model. The \code{\link{factor}} must a column in 
               the \code{\link{data.frame}} stored in the \code{asreml.obj} that is 
               a component of the supplied \code{asrtests.obj}.}
 \item{corr.funcs}{A single \code{character} string of length two that specifies the 
               \code{asreml} one-dimensional correlation or variance model function 
               for the row and column dimensions of a two-dimensional separable spatial 
               correlation model; the two-dimensional model is fitted as a random term. 
               If a correlation or variance model is not to be investigated for one 
               of the dimensions, specify \code{""} for that dimension.}
 \item{row.corrFitfirst}{If \code{TRUE} then the row correlation or variance 
               function is fitted first, followed by the addition of the column 
               correlation or variance function. If \code{FALSE}, the order of fitting 
               is reversed.}
 \item{dropRowterm}{A single \code{character} string nominating a \code{\link{factor}} 
               in the \code{\link{data.frame}} that has as many \code{levels} as there 
               are unique values in \code{row.covar}. This argument is required for 
               \code{spatial.model} set to \code{TPNCSS} or \code{TPPS}. It is used 
               to remove a term corresponding to the \code{dropRowterm} and a random 
               row deviations term based on \code{row.covar} will be included in the 
               model. If the argument is \code{NULL}, it is assumed that such a term 
               is not included in the fitted model stored in \code{asrtests.obj}.}
 \item{dropColterm}{A single \code{character} string nominating a \code{\link{factor}} 
               in the \code{\link{data.frame}} that has as many \code{levels} as there 
               are unique values in \code{col.covar}. This argument is required for 
               \code{spatial.model} set to \code{TPNCSS} or \code{TPPS}. It is used 
               to remove a term corresponding to the \code{dropColterm} and a random 
               column deviations term based on \code{col.covar} will be included in the 
               model. If the argument is \code{NULL}, it is assumed that such a term 
               is not included in the fitted model stored in \code{asrtests.obj}.}
 \item{nsegs}{A pair of \code{\link{numeric}} values giving the number of segments into
              which the column and row ranges are to be split, respectively (each value 
              specifies the number of internal knots + 1). If only one number is 
              specified, that value is used in both dimensions. If not specified, then 
              (number of unique values - 1) is used in each dimension; for a grid layout 
              with equal spacing, this gives a knot at each data value.}
 \item{nestorder}{A \code{\link{character}} of length 2. The order of nesting for column 
              and row dimensions, respectively; default=1 (no nesting). A value of 2 
              generates a spline with half the number of segments in that dimension, etc. 
              The number of segments in each direction must be a multiple of the order 
              of nesting.}
 \item{degree}{A \code{\link{character}}  of length 2. The degree of polynomial spline to 
              be used for column and row dimensions respectively; default=3.}
 \item{difforder}{A \code{\link{character}} of length 2. The order of differencing for 
              column and row dimensions, respectively; default=2.}
 \item{asreml.option}{A single \code{character} string specifying whether the \code{grp} or 
              \code{mbf} methods are to be used to supply externally formed covariate 
              matrices to \code{asreml}. If the \code{mbf} methods is to be used, then 
              \code{\link{makeTPPSplineMats.data.frame}} must be used before calling 
              \code{addSpatialModel.asrtests}. Compared to the \code{mbf} method, 
              the \code{grp} method creates large \code{\link{asrtests.object}}s, but is faster. 
              The \code{grp} method adds columns to the 
              \code{\link{data.frame}} containing the \code{data}; the \code{mbf} method 
              adds only fixed covariate to \code{data} and stores the random covariates 
              externally.}
 \item{tpps4mbf.obj}{An object made with \code{\link{makeTPPSplineMats.data.frame}} and 
              which contains the spline basis information, extra to that created by 
              \code{\link{makeTPPSplineMats.data.frame}}, that is needed to fit a 
              \code{TPPS} model using the \code{mbf} method of \code{asreml}.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
              even when it does not converge. If \code{FALSE} and the fit of the new 
              model does not converge, the supplied \code{asrtests.obj} is returned.  
              Also, if \code{FALSE} and the fit of the new model has converged, but that 
              of the old model has not, the new model will be accepted.}
 \item{allow.fixedcorrelation}{A \code{logical} indicating whether to accept a new model 
              even when it contains correlations in the model whose values have been 
              designated as fixed, bound or singular. If \code{FALSE} and the new model 
              contains correlations whose values have not been able to be estimated, 
              the supplied \code{asrtests.obj} is returned. The fit in the 
              \code{asreml.obj} component of the supplied \code{asrtests.obj} will 
              also be tested and a warning issued if both fixed correlations are found 
              in it and \code{allow.fixedcorrelation} is \code{FALSE}.}
 \item{checkboundaryonly}{If \code{TRUE} then boundary and singular terms are not removed by 
              \code{\link{rmboundary.asrtests}}; a warning is issued instead. Note that, for 
              correlation models, the fitting of each dimension and the test for a nugget 
              term are performed with \code{checkboundaryonly}  set to \code{TRUE} and its 
              supplied setting only honoured using a call to \code{\link{rmboundary.asrtests}} 
              immediately prior to returning the final result of the fitting.}
 \item{update}{If \code{TRUE} then \code{update.asreml} is called to fit the model 
               to be tested. In doing this the arguments \code{R.param} and 
              \code{G.param} are set to those in the \code{asreml} 
              object stored in \code{asrtests.obj} so that the values from the previous 
              model are used as starting values. If \code{FALSE} then a call is made to 
              \code{asreml} in which the only changes from the previous call are that 
              (i) models are modifed as specified and 
              (ii) modifications specified via \code{\dots} are made.}
 \item{IClikelihood}{A \code{character} that controls both the occurence and the type 
                    of likelihood for information criterion in the \code{test.summary} 
                    of the new \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML}, then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full}, 
                    then the AIC and BIC based on the full likelihood are included. 
                    (See also \code{\link{infoCriteria.asreml}}.)}
\item{\dots}{Further arguments passed to \code{\link{changeModelOnIC.asrtests}}, \code{asreml} and 
           \code{tpsmmb}.}
}
\details{
The model to which the spatial models is to be added is supplied in the \code{asrtests.obj}. It should not include terms that will be included in the local spatial model. All spatial model terms are fitted as fixed or random. Consequently, the residual model does not have to be iid. The improvement in the fit resulting from the addition of a spatial model to the supplied model is evaluated. Note that the data must be in the order that correponds to the \code{residual} argument with a variable to the right of another variable changes levels in the data frame faster than those of the other variable e.g. \code{Row:Column} implies that all levels for \code{Column} in consecutive rows of the \code{data.frame} with a single \code{Row} level. 

For the \code{corr} spatial model, the default model is an autocorrelation model of order one (\code{ar1}) for each dimension. However, any of the single dimension correlation/variance models from \code{asreml} can be specified for each dimension, as can no correlation model for a dimension; the models for the two dimensions can differ. A series of models are tried, without removing boundary or singular terms, beginning with the addition of row correlation and followed by the addition of column correlation or, if the \code{row.corrFitfirst} is set to \code{FALSE}, the reverse order. If the fitting of the first-fitted correlation did not result in a model change because the fitting did not converge or correlations were fixed, but the fit of the second correlation was successful, then adding the first correlation will be retried. If one of the metric correlation functions is specified (e.g. \code{exp}), then the \code{row.covar} or \code{col.covar} will be used in the spatial model. However, because the correlations are fitted separately for the two dimensions, the \code{row.factor} and \code{col.factor} are needed for all models and is used for a dimension that does not involve a correlation/variance function for the fit being performed. Also, the correlation models are fitted as \code{random} terms and so the correlation model will include a variance parameter for the grid even when \code{ar1} is used to specify the correlation model, i.e. the model fitted is a variance model and there is no difference between \code{ar1} and \code{ar1v} in fitting the model. The variance parameter for this term represents the spatial variance and the fit necessarily includes a nugget term, this being the residual variance. The nuggest term can be removed by fixing the residual variance to one using \code{\link{setvarianceterms.call}}, provided heterogeneous residual variances have not been specified. Once the fitting of the correlation model has been completed, the \code{rmboundary} function will be executed with the \code{checkboundaryonly} value suppied in the \code{addSpatialModel.asrtests} call.

The \code{TPNCSS} spatial model is as decribed by Verbyla et al. (2018) and the \code{TPPS} model is as described by Rodriguez-Alvarez et al. (2018). However, for the \code{TPPS} model, the degree of the polynomial  and the order of differencing can be varied. The defaults of 3 and 2, respectively, fit a cubic spline with second order differencing, while setting both the degree and order of differencing to 1 will fit a type of linear variance model (Piepho, Boer and Williams, 2022) The fixed terms for these models are \code{row.covar + col.covar + row.covar:col.covar} and the random terms are \code{spl(row.covar) + spl(col.covar) + dev(row.covar) + dev(col.covar) + spl(row.covar):col.covar +}\cr \code{row.covar:spl(col.covar) + spl(row.covar):spl(col.covar)}. The supplied model should not include any of these terms. However, any fixed or random main-effect term for either \code{dropRowterm} or \code{dropColterm} will be removed from the fit.

The \code{TPPS} model is fitted using the function \code{tpsmmb} from the \code{R} package \code{TPSbits} authored by Sue Welham (2022). There are two methods for supplying the spline basis information produced by \code{tpsmmb} to \code{asreml}. The \code{grp} method adds the it to the \code{data.frame} holding the information for the analysis. The \code{mbf} method requires the spline basis information to be in the same enviroment as the function that is called to make a fit using \code{asreml}. To this end and prior to invoking the calling function, \code{\link{makeTPPSplineMats.data.frame}} must be used produce the \code{\link{data.frame}s}.

All models utlize the function \code{\link{changeTerms.asrtests}} to fit the spatial model. Arguments from \code{tpsmmb} and \code{\link{changeTerms.asrtests}} can be supplied in calls to \code{addSpatialModel.asrtests} and will be passed on to the relevant function through the ellipses argument (\ldots).

The data for experiment can be divided \code{sections} and the same spatial model fitted separately to each. The fit over all of the sections is assessed.

Each combination of a row.coords and a col.coords does not have to specify a single observation; 
for example, to fit a local spatial model to the main units of a split-unit design, 
each combination would correspond to a main unit and all subunits of the main unit would would 
have the same combination.
}

\value{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary} for the model that includes the 
       spatial model, unless the spatial model fails to be fitted when \code{allow.unconverged} 
       and/or \code{allow.fixedcorrelation} is set to \code{FALSE}.}

\references{Piepho, H.-P., Boer, M. P., & Williams, E. R. (2022). Two-dimensional P-spline smoothing for spatial analysis of plant breeding trials. \emph{Biometrical Journal}, \bold{64}, 835-857.

Rodriguez-Alvarez, M. X., Boer, M. P., van Eeuwijk, F. A., & Eilers, P. H. C. (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. \emph{Spatial Statistics}, \bold{23}, 52-71. 

Verbyla, A. P., De Faveri, J., Wilkie, J. D., & Lewis, T. (2018). Tensor Cubic Smoothing Splines in Designed Experiments Requiring Residual Modelling. \emph{Journal of Agricultural, Biological and Environmental Statistics}, \bold{23(4)}, 478-508.

Welham, S. J. (2022) \code{TPSbits}: \emph{Creates Structures to Enable Fitting and Examination of 2D Tensor-Product Splines using ASReml-R.} Version 1.0.0 \url{https://mmade.org/tpsbits/}} 

\author{Chris Brien}

\seealso{\code{\link{as.asrtests}}, \code{\link{rmboundary.asrtests}}, 
         \code{\link{addSpatialModelOnIC.asrtests}}, \code{\link{testranfix.asrtests}}, \cr
         \code{\link{testresidual.asrtests}}, 
         \code{\link{newfit.asreml}}, \code{\link{reparamSigDevn.asrtests}}, 
         \code{\link{changeTerms.asrtests}}, \cr 
         \code{\link{chooseSpatialModelOnIC.asrtests}}, \code{\link{changeTerms.asrtests}}, 
         \code{\link{infoCriteria.asreml}}}
\examples{\dontrun{

data(Wheat.dat)

#Add row and column covariates
Wheat.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })

#Fit initial model
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column,
                      data=Wheat.dat)

#Create an asrtests object, removing boundary terms
current.asrt <- as.asrtests(current.asr, NULL, NULL, 
                            label = "Random Row and Column effects")
current.asrt <- rmboundary(current.asrt)

current.asrt <- addSpatialModel(current.asrt, spatial.model = "TPPS", 
                                row.covar = "cRow", col.covar = "cColumn",
                                dropRowterm = "Row", dropColterm = "Column",
                                asreml.option = "grp")
infoCriteria(current.asrt$asreml.obj)
}}
\keyword{asreml}
