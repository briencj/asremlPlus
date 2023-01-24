\name{chooseSpatialModelOnIC.asrtests}
\alias{chooseSpatialModelOnIC.asrtests}
\alias{chooseSpatialModelOnIC}
\title{Uses information criteria to choose the best fitting spatial model for accounting for local spatial variation.}
\description{For a response variable measured on a potentially irregular grid of rows and 
       columns of the units, uses information criteria to decide whether to add to the fitted 
       model stored in the supplied \code{\link{asrtests.object}} either a two-dimensional 
       exponential correlation model, a two-dimensional tensor-product natural cubic smoothing 
       spline model (TPNCSS), a two-dimensional tensor-product penalized P-spline model (TPPCS) 
       model, or a two-dimensional tensor-product penalized linear spline model with 
       first-difference penalties (TPP1LS) to account for the local spatial variation. The models 
       from which to select can be reduced to a subset of these four models. The data can be 
       arranged in sections, for each of which there is a grid and for which the model is to 
       be fitted separately. Also, the rows and columns of a grid are not necessarily one 
       observational unit wide. The spatial model is only added if the information criterion 
       of the supplied model is decreased with the addition of the local spatial model.

       One or more rows is added to the \code{test.summary} \code{data.frame} of the 
       \code{\link{asrtests.object}}, for each \code{section} and each spatial model, 
       stating whether or not the new model has been swapped for a model in which the 
       spatial model has been added to the  supplied model. Convergence in fitting the 
       model is checked and a note included in the \code{action} if there was not. 
       All components of the \code{\link{asrtests.object}} are updated to exhibit the 
       differences between the supplied and any new model.
}
\usage{
\method{chooseSpatialModelOnIC}{asrtests}(asrtests.obj, trySpatial = "all", 
                       sections = NULL, 
                       row.covar = "cRow", col.covar = "cCol", 
                       row.factor = NULL, col.factor = NULL, 
                       nsegs = NULL, nestorder = c(1,1), 
                       asreml.option = "mbf",  
                       tpps4mbf.obj = NULL, 
                       allow.unconverged = FALSE, 
                       allow.fixedcorrelation = FALSE,
                       checkboundaryonly = FALSE, update = FALSE, 
                       IClikelihood = "full", which.IC = "AIC", 
                       return.asrts = "best", ...)}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components 
       (i) \code{asreml.obj}, (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{trySpatial}{A \code{character} string nominating the types of spatial 
               model whose fits are to be assessed. Possible values are \code{corr}, 
               \code{TPNCSS}, \code{TPPCS}, and \code{TPP1LS}.}
 \item{sections}{A single \code{character} string that species the name of the column 
               in the \code{\link{data.frame}} that contains the \code{\link{factor}} 
               that identifies different sections of the data 
               to which separate spatial models are to be fitted.}
 \item{row.covar}{A single \code{character} string nominating a \code{\link{numeric}} 
               column in the \code{\link{data.frame}} that contains the values of a 
               centred covariate indexing the rows of the grid.}
 \item{col.covar}{A single \code{character} string nominating a \code{\link{numeric}} 
               column in the \code{\link{data.frame}} that contains the values of a 
               centred covariate indexing the columns of the grid.}
 \item{row.factor}{A single \code{character} string nominating a \code{\link{factor}} 
               in the \code{\link{data.frame}} that has as many \code{levels} as there 
               are unique values in \code{row.covar}. This argument is required for 
               \code{spatial.model} set to \code{TPNCSS} or \code{TPPS}. It is used 
               to remove a term corresponding to the \code{row.factor} and a random 
               row deviations term based on \code{row.covar} will be included in the 
               model. If the argument is \code{NULL}, it is assumed that such a term 
               is not included in the fitted model stored in \code{asrtests.obj}.}
 \item{col.factor}{A single \code{character} string nominating a \code{\link{factor}} 
               in the \code{\link{data.frame}} that has as many \code{levels} as there 
               are unique values in \code{col.covar}. This argument is required for 
               \code{spatial.model} set to \code{TPNCSS} or \code{TPPS}. It is used 
               to remove a term corresponding to the \code{col.factor} and a random 
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
 \item{asreml.option}{A single \code{character} string specifying whether the \code{grp} or 
              \code{mbf} methods are to be used to supply externally formed covariate 
              matrices to \code{asreml}. If the \code{mbf} methods is to be used, then 
              \code{\link{makeTPSPlineXZMats.data.frame}} must be used before calling 
              \code{addSpatialModelOnIC.asrtests}. Compared to the \code{mbf} method, 
              the \code{grp} method creates large \code{asreml} objects, but is faster. 
              The \code{grp} method adds columns to the 
              \code{\link{data.frame}} containing the \code{data}; the \code{mbf} method 
              adds only fixed covariate to \code{data} and stores the random covariates 
              externally.}
 \item{tpps4mbf.obj}{An object made with \code{\link{makeTPSPlineXZMats.data.frame}} and 
              which contains the spline basis information, that is extra to the 
              \code{\link{data.frame}s} created by \cr\code{\link{makeTPSPlineXZMats.data.frame}} 
              in the environment in which it is called and that is needed to fit a 
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
              \code{\link{rmboundary.asrtests}}; a warning is issued instead.}
 \item{update}{If \code{TRUE} then \code{update.asreml} is called to fit the model 
              to be tested. In doing this the arguments \code{R.param} and 
              \code{G.param} are set to those in the \code{asreml} 
              object stored in \code{asrtests.obj} so that the values from the previous 
              model are used as starting values. If \code{FALSE} then a call is made to 
              \code{asreml} in which the only changes from the previous call are that 
              (i) models are modifed as specified and 
              (ii) modifications specified via \code{\dots} are made.}
 \item{which.IC}{A \code{character} specifying the information criterion to be used in 
              selecting the best model. Possible values are \code{AIC} and 
              \code{BIC}. The values of the criterion for supplied model must 
              exceed that for changed model for the changed model to be returned.}
 \item{IClikelihood}{A \code{character} specifying whether Restricted Maximum Likelihood 
               (\code{REML}) or the full likelihood (\code{full}) are to be used in 
               calculating the information criteria.}
 \item{return.asrts}{A \code{character} string specifying whether the \code{\link{asrtests.object}} 
               for the \code{best} fitting model (smallest AIC or BIC) is returned or the 
               \code{\link{asrtests.object}}s resulting from the attempted fits of \code{all} of 
               the models specifed using \code{trySpatial} are returned.}
\item{\dots}{Further arguments passed to \code{\link{changeModelOnIC.asrtests}}, \code{asreml} and 
           \code{tpsmmb}.}
}

\details{
A fitted spatial model is only returned if it improves the fit over an above that achieved with the model fit supplied in the \code{asrtests.obj}. If \code{return.asrts} is \code{all}, then this applies to each spatial model specified by \code{trySpatial}. The model fit supplied in the \code{asrtests.obj} should not include terms that will be included in any local spatial model. All spatial model terms are fitted as fixed or random. Consequently, the residual model does not have to be iid. The improvement in the fit resulting from the addition of a spatial model to the supplied model is evaluated. 

For the \code{corr} spatial model, an exponential model (\code{exp}) is used for each dimension to model the spatial correlation. A series of models are tried, beginning with the addition of row correlation and followed by the addition of column correlation. Only if the model fit is improved is a correlation retained. Finally, if any correlation is retained, the improvment to the fit of a nuggest term is assessed. In this model, the correlation between observations from different rows is the correlation between observations in adjacent rows raised to the power equal to the absolute value of the difference in their \code{row.covar} values; similarly for the correlation in the column dimension. 

The tensor-product natural-cubic-smoothing-spline \code{TPNCSS} spatial model is as decribed by Verbyla et al. (2018), the tensor-product penalized-cubic-spline \code{TPPCS} model is as described by Rodriguez-Alvarez et al. (2018), and the tensor-product, first-difference-penalty, linear spline \code{TPP1LS} model that is amongst those described by Piepho, Boer and Williams (2022).  The fixed terms for these models are \code{row.covar + col.covar + row.covar:col.covar} and the random terms \code{spl(row.covar) + spl(col.covar) + dev(row.covar) + dev(col.covar) + spl(row.covar):col.covar + row.covar:spl(col.covar) + spl(row.covar):spl(col.covar)}. The supplied model should not include any of these terms. However, any fixed or random main-effect term for either \code{row.factor} or \code{col.factor} will be removed from the fit.

The \code{TPPCS} and \code{TPP1LS} models are fitted using the function \code{tpsmmb} from the \code{R} package \code{TPSbits} authored by Sue Welham (2022). There are two methods for supplying the spline basis information produced by \code{tpsmmb} to \code{asreml}. The \code{grp} method adds the it to the \code{data.frame} holding the information for the analysis. The \code{mbf} method requires the spline basis information to be in the same environment as the function that is called to make a fit using \code{asreml}. To this end, and prior to invoking the calling function, \code{\link{makeTPSPlineXZMats.data.frame}} must be used produce the \code{\link{data.frame}s}.

All models utlize the function \code{\link{changeModelOnIC.asrtests}} to assess the model fit, the information critera used in assessing the fit being calculated using \code{\link{infoCriteria}}. Arguments from \code{tpsmmb} and \code{\link{changeModelOnIC.asrtests}} can be supplied in calls to \code{chooseSpatialModelOnIC.asrtests} and will be passed on to the relevant function throught the ellipses argument (\ldots).

The data for experiment can be divided into \code{sections} and an attempt to fit the same spatial model to each is made. The fit may differ for each of the \code{sections}, but the fit over all of the \code{sections} is assessed.

Each combination of a row.coords and a col.coords does not have to specify a single observation; 
for example, to fit a local spatial model to the main units of a split-unit design, 
each combination would correspond to a main unit and all subunits of the main unit would would 
have the same combination.
}

\value{A \code{\link{list}} containing four components: (i) \code{asrts}, (ii) \code{spatial.IC}, 
       (iii) \code{best.spatial.mod}, and (iv) \code{best.spatial.IC}. 
       
       The component \code{asrts} itself holds a  \code{\link{list}} of one or more 
       \code{\link{asrtests.object}}s, either the best overall out of the supplied model and 
       the spatial models, or, for each spatial model, the best out of the supplied model 
       and that spatial model. Each \code{\link{asrtests.object}} contains the components:  
       (i) \code{asreml.obj}, (ii) \code{wald.tab}, and (iii) \code{test.summary}. 
       
       The \code{spatial.IC} component holds a \code{\link{data.frame}} with summary of the 
       values of the information criteria for the supplied model and those resulting from  
       adding the spatial model to the supplied model (if the spatial model did not fit, then 
       all values will be \code{NA}). 
       
       The \code{best.spatial} component is a character giving the name of the best spatial 
       model, and \code{best.spatial.AIC} gives the value of its \code{AIC}.
}

\references{Piepho, H.-P., Boer, M. P., & Williams, E. R. (2022). Two-dimensional P-spline smoothing for spatial analysis of plant breeding trials. \emph{Biometrical Journal}, \bold{64}, 835-857.

Rodriguez-Alvarez, M. X., Boer, M. P., van Eeuwijk, F. A., & Eilers, P. H. C. (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. \emph{Spatial Statistics}, \bold{23}, 52-71. 

Verbyla, A. P., De Faveri, J., Wilkie, J. D., & Lewis, T. (2018). Tensor Cubic Smoothing Splines in Designed Experiments Requiring Residual Modelling. \emph{Journal of Agricultural, Biological and Environmental Statistics}, \bold{23(4)}, 478-508.

Welham, S. J. (2022) \code{TPSbits}: \emph{Creates Structures to Enable Fitting and Examination of 2D Tensor-Product Splines using ASReml-R.} Version 1.0.0 \url{https://mmade.org/tpsbits/}} 

\author{Chris Brien}

\seealso{\code{\link{as.asrtests}}, \code{\link{rmboundary.asrtests}}, 
         \code{\link{addSpatialModelOnIC.asrtests}}, \code{\link{addSpatialModel.asrtests}}, ,
         \code{\link{testranfix.asrtests}}, \cr 
         \code{\link{testresidual.asrtests}}, 
         \code{\link{newfit.asreml}}, \code{\link{reparamSigDevn.asrtests}}, 
         \code{\link{chooseModel.asrtests}}, \cr 
         \code{\link{changeModelOnIC.asrtests}}, \code{\link{infoCriteria.asreml}}}
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

# Choose the best of three models the local spatial variation
current.asrt <- chooseSpatialModelOnIC(current.asrt, 
                                       row.covar = "cRow", col.covar = "cColumn",
                                       row.factor = "Row", col.factor = "Column",
                                       asreml.option = "grp")
}}
\keyword{asreml}
