\name{addSpatialModelOnIC.asrtests}
\alias{addSpatialModelOnIC.asrtests}
\alias{addSpatialModelOnIC}
\title{Uses information criteria to decide whether to add a spatial model to account for local spatial variation.}
\description{Adds either a correlation, two-dimensional tensor-product natural cubic 
       smoothing spline (TPNCSS), or a two-dimensional tensor-product penalized P-spline 
       model (TPPS) to account for the the local spatial variation exhibited by a response variable 
       measured on a potentially irregular grid of rows and columns of the units. The data may 
       be arranged in sections for each of which there is a grid and for which the model is to 
       be fitted separately. Also, the rows and columns of a grid are not necessarily one 
       observational unit wide. The spatial model is only added if the information criterion of 
       the supplied model is decreased with the addition of the local spatial model.  
       
       A row is added for each \code{section} to the \code{test.summary} \code{data.frame} 
       of the \code{\link{asrtests.object}} stating whether or not the new model has been 
       swapped for a model in which the spatial model has been add to the  supplied model. 
       Convergence and the occurrence of fixed correlations in fitting the 
       model is checked and a note included in the \code{action} if there was not. 
       All components of the \code{\link{asrtests.object}} are updated to exhibit the 
       differences between the supplied and the new model, if a spatial model is added.
}
\usage{
\method{addSpatialModelOnIC}{asrtests}(asrtests.obj, spatial.model = "TPPS", 
                    sections = NULL, 
                    row.covar = "cRow", col.covar = "cCol", 
                    row.factor = "Row", col.factor = "Col", 
                    corr.funcs = c("ar1", "ar1"), 
                    row.corrFitfirst = TRUE, 
                    dropRowterm = NULL, dropColterm = NULL, 
                    nsegs = NULL, nestorder = c(1,1), 
                    degree = c(3,3), difforder = c(2,2), 
                    rotateX = FALSE, ngridangles = c(18, 18), 
                    which.rotacriterion = "AIC", nrotacores = 1, 
                    asreml.option = "mbf", tpps4mbf.obj = NULL, 
                    allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                    checkboundaryonly = FALSE, update = FALSE, 
                    IClikelihood = "full", which.IC = "AIC", ...)
}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components 
       (i) \code{asreml.obj}, (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{spatial.model}{A single \code{character} string nominating the type of spatial 
               model to fit. Possible values are \code{corr}, \code{TPNCSS} and 
               \code{TPPS}.}
 \item{sections}{A single \code{character} string that specifies the name of the column 
               in the \code{\link{data.frame}} that contains the \code{\link{factor}} 
               that identifies different sections of the data to which separate spatial 
               models are to be fitted. Note that, if there are other terms that 
               involve \code{sections} in the random formula, there should be separate 
               terms for each level of \code{sections}. For example, in a blocked 
               experiment involving multiple sites, there should be a separate term for 
               the Blocks at each Site i.e. terms of the form at(Site, i):Block, or 
               equivalent, for each site. Otherwise, the combined term (e.g. Site:Block) 
               will affect the fitting of the local spatial model for a site.}
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
               correlation model to be fitted when \code{spatial.model} is \code{corr}; 
               the two-dimensional model is fitted as a random term. If a correlation 
               or variance model is not to be investigated for one of the dimensions, 
               specify \code{""} for that dimension.}
 \item{row.corrFitfirst}{If \code{TRUE} then, in fitting the model for 
                \code{spatial.model} set to \code{corr}, the row correlation or 
               variance function is fitted first, followed by the addition of the 
               column correlation or variance function. If \code{FALSE}, the order 
               of fitting is reversed.}
 \item{dropRowterm}{A single \code{character} string nominating a \code{\link{factor}} 
               in the \code{\link{data.frame}} that has as many \code{levels} as there 
               are unique values in \code{row.covar}. This argument is applicable for 
               \code{spatial.model} set to \code{TPNCSS} or \code{TPPS}. It is used 
               to remove a term corresponding to the \code{dropRowterm} and a random 
               row deviations term based on \code{row.covar} will be included in the 
               model. If the argument is \code{NULL}, it is assumed that such a term 
               is not included in the fitted model stored in \code{asrtests.obj}.}
 \item{dropColterm}{A single \code{character} string nominating a \code{\link{factor}} 
               in the \code{\link{data.frame}} that has as many \code{levels} as there 
               are unique values in \code{col.covar}. This argument is applicable for 
               \code{spatial.model} set to \code{TPNCSS} or \code{TPPS}. It is used 
               to remove a term corresponding to the \code{dropColterm} and a random 
               column deviations term based on \code{col.covar} will be included in the 
               model. If the argument is \code{NULL}, it is assumed that such a term 
               is not included in the fitted model stored in \code{asrtests.obj}.}
 \item{nsegs}{A pair of \code{\link{numeric}} values giving the number of segments into
              which the column and row ranges are to be split, respectively, for fitting 
              a P-spline model (\code{TPPS}) (each value specifies the number of internal 
              knots + 1). If not specified, then (number of unique values - 1) is used 
              in each dimension; for a grid layout with equal spacing, this gives a knot 
              at each data value. If \code{sections} is not \code{NULL} and the grid 
              differs between the \code{sections}, then \code{nsegs} will differ between 
              the \code{sections}.}
 \item{nestorder}{A \code{\link{numeric}} of length 2. The order of nesting for column 
              and row dimensions, respectively, in fitting a P-spline model (\code{TPPS}). 
              A value of 1 specifies no nesting, a value of 2 generates a spline with 
              half the number of segments in that 
              dimension, etc. The number of segments in each direction must be a 
              multiple of the order of nesting.}
 \item{degree}{A \code{\link{numeric}}  of length 2. The degree of polynomial spline to 
              be used for column and row dimensions respectively, in fitting a P-spline 
              (\code{TPSS}).}
 \item{difforder}{A \code{\link{numeric}} of length 2. The order of differencing for 
              column and row dimensions, respectively, in fitting a P-spline 
              (\code{TPSS}).}
 \item{rotateX}{A \code{\link{logical}} indicating whether to rotate the eigenvectors of the 
              penalty matrix, as described by Piepho, Boer and Williams (2022), when 
              fitting a P-spline (\code{TPSS}). Setting \code{rotateX} to \code{TRUE} 
              results in an optimized rotation produced by exploring 
              a two-dimensional grid of rotation angle pairs, and for each pair 
              analyzing the data under a model that omits the random interaction terms. 
              The angle pair with the minimum deviance is used to apply an optimized 
              rotation. Rotation of the eigenvectors is only relevant for 
              \code{difforder} values greater than 1 and has only been implemented for 
              \code{difforder} equal to 2.}
 \item{ngridangles}{A \code{\link{numeric}} of length 2. The numbers of angles between 0 
              and 90 degrees for row and column dimensions to be used in determining the 
              optimal pair of angles for rotating the eignevectors of the penalty matrix 
              of a P-spline (\code{TPSS}). Specifying factors of 90 will result in 
              integer-valued angles. The number of grid points, and hence 
              re-analyses will be the product of the two sums. a sum being one plus each 
              of the values for \code{ngridangles}.}
 \item{which.rotacriterion}{A single \code{character} string nominating which of the 
              criteria out of the \code{deviance}, the \code{likelihood}, the 
              \code{AIC} and the \code{BIC} in determining the optimal rotation of the 
              eigenvectors of the penalty matrix. The \code{deviance} uses the 
              REML value computed by \code{asreml}; the other criteria use the 
              full likelihood, evaluated using the REML estimates, that is computed 
              by \code{\link{infoCriteria.asreml}}.}
 \item{nrotacores}{A \code{\link{numeric}} specifying the number of cores to deploy for 
              running the analyses required to search the two-diemsional grid of rotation 
              angles when \code{rotateX} is \code{TRUE}. Parallel processing has been 
              implemented for analyzing, for each column angle, the set of angles to be 
              investigated for the row dimension.  The default value of one means that 
              parallel processing will not be used. The value chosen for \code{nrotacores} 
              needs to balanced against the other processes that are using parallel 
              processing at the same time.}
 \item{asreml.option}{A single \code{character} string specifying whether the \code{grp} or 
              \code{mbf} methods are to be used to supply externally formed covariate 
              matrices to \code{asreml} when fitting a P-spline (\code{TPSS}). Compared 
              to the \code{mbf} method, the \code{grp} method is somewhat faster, but 
              creates large \code{\link{asrtests.object}}s for which the time it takes 
              to save them can exceed any gains in execution speed. The \code{grp} 
              method adds columns to the \code{\link{data.frame}} containing the \code{data}. 
              On the other hand, the \code{mbf} method adds only the fixed covariates 
              to \code{data} and stores the random covariates in the environment of the 
              internal function that calls the spline-fitting function; there are three smaller 
              \code{\link{data.frame}s} for each \code{section} that are not stored in the 
              \code{asreml.object} resulting from the fitted model.}
 \item{tpps4mbf.obj}{An object made with \code{\link{makeTPPSplineMats.data.frame}} 
              that contains the spline basis information for fitting P-splines. 
              The argument \code{tpps4mbf.obj} only needs to be set when the \code{mbf} 
              option of \code{asreml.option} is being used and it is desired to use 
              \code{mbf} \code{\link{data.frame}s} that have been created and stored 
              prior to calling \cr 
              \code{addSpatialModelOnIC.asrtests}. If \code{tpps4mbf.obj} 
              is \code{NULL}, \code{\link{makeTPPSplineMats.data.frame}} will be called 
              internally to produce the required \code{mbf} \code{\link{data.frame}s}.}
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
 \item{update}{If \code{TRUE}, then \code{\link{newfit.asreml}} is called to fit the 
              model to be tested, using the values of the variance parameters stored in 
              the \code{asreml.object}, that is stored in \code{asrtests.obj}, as starting values. 
              If \code{FALSE}, then \code{\link{newfit.asreml}} will not use the stored 
              variance parameter values as starting values when fitting the new model, 
              the only modifications being (i) to add the terms for the spatial models 
              and (ii) those specified via \code{\dots}.}
 \item{which.IC}{A \code{character} specifying the information criterion to be used in 
                   selecting the best model. Possible values are \code{AIC} and 
                   \code{BIC}. The value of the criterion for supplied model must 
                   exceed that for changed model for the changed model to be returned. 
                   (For choosing the rotation angle of the eigenvectors of the penalty 
                    matrix, see \code{which.rotacriterion}.}
 \item{IClikelihood}{A \code{character} specifying whether Restricted Maximum Likelihood 
                   (\code{REML}) or the full likelihood, evaluated using REML estimates, 
                   (\code{full}) are to be used in calculating the information criteria 
                   to be included in the \code{test.summary} of an 
                   \code{\link{asrtests.object}} or to be used in choosing the best model.}
\item{\dots}{Further arguments passed to \code{\link{changeModelOnIC.asrtests}}, \code{asreml} and 
           \code{tpsmmb}.}
}
\details{
A fitted spatial model is only returned if it improves the fit over and above that of achieved with the model fit supplied in the \code{asrtests.obj}. To fit the spatial model without any hypothoses testing or comparison of informtion criteria use \code{\link{addSpatialModel.asrtests}}. The model fit supplied in the \code{asrtests.obj} should not include terms that will be included in the local spatial model. All spatial model terms are fitted as fixed or random. Consequently, the residual model does not have to be iid. Note that the data must be in the order that correponds to the \code{residual} argument with a variable to the right of another variable changes levels in the data frame faster than those of the other variable e.g. \code{Row:Column} implies that all levels for \code{Column} in consecutive rows of the \code{data.frame} with a single \code{Row} level.

For the \code{corr} spatial model, the default model is an autocorrelation model of order one (\code{ar1}) for each dimension. However, any of the single dimension correlation/variance models from \code{asreml} can be specified for each dimension, as can no correlation model for a dimension; the models for the two dimensions can differ. Using a forward selection procedure, a series of models are tried, without removing boundary or singular terms, beginning with the addition of row correlation and followed by the addition of column correlation or, if the \code{row.corrFitfirst} is set to \code{FALSE}, the reverse order. If the fitting of the first-fitted correlation did not result in a model change because the fitting did not converge or correlations were fixed, but the fit of the second correlation was successful, then adding the first correlation will be retried. If one of the metric correlation functions is specified (e.g. \code{exp}), then the \code{row.covar} or \code{col.covar} will be used in the spatial model. However, because the correlations are fitted separately for the two dimensions, the \code{row.factor} and \code{col.factor} are needed for all models and is used for a dimension that does not involve a correlation/variance function for the fit being performed. Also, the correlation models are fitted as \code{random} terms and so the correlation model will include a variance parameter for the grid even when \code{ar1} is used to specify the correlation model, i.e. the model fitted is a variance model and there is no difference between \code{ar1} and \code{ar1v} in fitting the model. The variance parameter for this term represents the spatial variance and the fit necessarily includes a nugget term, this being the residual variance. If any correlation is retained, the need for a nuggest term is assessed by fixing the residual variance to one, which will have no effect if heterogeneous residual variances have been specified. Once the fitting of the correlation model has been completed, the \code{rmboundary} function will be executed with the \code{checkboundaryonly} value suppied in the \code{addSpatialModelOnIC.asrtests} call.

The tensor-product natural-cubic-smoothing-spline (\code{TPNCSS}) spatial models are as decribed by Verbyla et al. (2018). The tensor-product penalized-spline (\code{TPPS}) models are as described by (Piepho, Boer and Williams, 2022). For the \code{TPPS} model, the degree of the polynomial and the order of differencing can be varied. The defaults of 3 and 2, respectively, fit a cubic spline with second order differencing, which is similar to those of Rodriguez-Alvarez et al. (2018). Setting both the degree and order of differencing to 1 will fit a type of linear variance model. The fixed terms for the spline models are \code{row.covar + col.covar + row.covar:col.covar} and the random terms are \code{spl(row.covar) + spl(col.covar) + dev(row.covar) + dev(col.covar) + spl(row.covar):col.covar +}\cr \code{row.covar:spl(col.covar) + spl(row.covar):spl(col.covar)}, except that \code{spl(row.covar) + spl(col.covar)} is not included in \code{TPPS} models. The supplied model should not include any of these terms. However, any fixed or random main-effect Row or Column term that has been included as an initial model for comparison with a spatial model can be removed prior to fitting the spatial model using \code{dropRowterm} or \code{dropColterm}.  For \code{TPPS} with second-order differencing, the model matrices used to fit the random terms \code{spl(row.covar):col.covar} and \code{row.covar:spl(col.covar)} are transformed using the spectral decomposition of their penalty matrices, and unstructured variance models across the columns of the linear component of each of them will be tried.

The \code{TPPCS} and \code{TPP1LS} models are fitted using functions from the \code{R} package \code{TPSbits} authored by Sue Welham (2022). There are two methods for supplying the spline basis information produced by \code{tpsmmb} to \code{asreml}. The \code{grp} method adds it to the \code{data.frame} supplied in the \code{data} argument of the \code{asreml} call. The \code{mbf} method creates smaller \code{\link{data.frame}s} with the spline basis information in the same environment as the internal function that calls the spline-fitting function. For \code{TPPS} with second-order differencing, it is also possible to investigate the rotation of the penalty matrix eigenvectors for the random terms \code{spl(row.covar):col.covar} and \code{row.covar:spl(col.covar)} (for more information see Piepho, Boer and Williams, 2022).  If it is desired to use in a later session, an \code{asreml} function, or \code{asrtests} function that calls \code{asreml}, (e.g. \code{predict.asreml}, \code{\link{predictPlus.asreml}}, or \code{\link{changeTerms.asrtests}}) on an \code{asreml.object} created using \code{mbf} terms, then the \code{mbf} \code{\link{data.frame}s} will need to be recreated using \code{\link{makeTPPSplineMats.data.frame}} in the new session, supplying, if there has been rotation of the penalty matrix eigenvectors, the \code{theta} values that are returned as the attribute \code{theta.opt} of the \code{asreml.obj}. 

All models utlize the function \code{\link{changeModelOnIC.asrtests}} to assess the model fit, the information critera used in assessing the fit being calculated using \code{\link{infoCriteria}}. Any \code{bound} terms are remopved from the model. Arguments from \code{tpsmmb} and \code{\link{changeModelOnIC.asrtests}} can be supplied in calls to \code{addSpatialModelOnIC.asrtests} and will be passed on to the relevant function through the ellipses argument (\ldots).

The data for experiment can be divided \code{sections} and the same spatial model fitted separately to each. The fit over all of the sections is assessed.

Each combination of a row.coords and a col.coords does not have to specify a single observation; 
for example, to fit a local spatial model to the main units of a split-unit design, 
each combination would correspond to a main unit and all subunits of the main unit would would 
have the same combination.
}

\value{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj}, 
       possibly with attribute \code{theta.opt},
       (ii) \code{wald.tab}, and (iii) \code{test.summary} for the model whose fit has 
       the smallest information criterion between the supplied and spatial model. The values 
       of the degrees of freedom and the information criteria in the \code{test.summary} are 
       differences between those of the changed model and those of the model supplied to 
       \code{addSpatialModelOnIC}. If the 
       \code{\link{asrtests.object}} is the result of fitting a \code{TPPCS} model with 
       an exploration of the rotation of the eigenvectors of the penalty matrix for the linear 
       components, then the \code{asreml.obj} will have an attribute \code{theta.opt} that contains 
       the optimal rotation angles of the eigenvectors.}
\references{Piepho, H.-P., Boer, M. P., & Williams, E. R. (2022). Two-dimensional P-spline smoothing for spatial analysis of plant breeding trials. \emph{Biometrical Journal}, \bold{64}, 835-857.

Rodriguez-Alvarez, M. X., Boer, M. P., van Eeuwijk, F. A., & Eilers, P. H. C. (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. \emph{Spatial Statistics}, \bold{23}, 52-71. 

Verbyla, A. P., De Faveri, J., Wilkie, J. D., & Lewis, T. (2018). Tensor Cubic Smoothing Splines in Designed Experiments Requiring Residual Modelling. \emph{Journal of Agricultural, Biological and Environmental Statistics}, \bold{23(4)}, 478-508.

Welham, S. J. (2022) \code{TPSbits}: \emph{Creates Structures to Enable Fitting and Examination of 2D Tensor-Product Splines using ASReml-R.} Version 1.0.0 \url{https://mmade.org/tpsbits/}} 

\author{Chris Brien}

\seealso{\code{\link{as.asrtests}}, 
         \code{\link{makeTPPSplineMats.data.frame}}, 
         \code{\link{addSpatialModel.asrtests}}, \cr
         \code{\link{chooseSpatialModelOnIC.asrtests}}, 
         \code{\link{changeModelOnIC.asrtests}}, 
         \code{\link{changeTerms.asrtests}},  \cr
         \code{\link{rmboundary.asrtests}},  
         \code{\link{testranfix.asrtests}}, 
         \code{\link{testresidual.asrtests}}, 
         \code{\link{newfit.asreml}}, \cr
         \code{\link{reparamSigDevn.asrtests}}, 
         \code{\link{changeTerms.asrtests}}, 
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

current.asrt <- addSpatialModelOnIC(current.asrt, spatial.model = "TPPS", 
                                    row.covar = "cRow", col.covar = "cColumn",
                                    dropRowterm = "Row", dropColterm = "Column",
                                    asreml.option = "grp")
infoCriteria(current.asrt$asreml.obj)
}}
\keyword{asreml}
