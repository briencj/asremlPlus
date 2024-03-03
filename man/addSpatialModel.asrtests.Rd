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
       observational unit wide. For \code{TPPS} models for which the order of differencing the 
       penalty matrix is two, the an optimal rotation of the null-space eigenvectors of the 
       penalty matrix can be investigated.
       
       No hypothesis testing or comparison of information criteria is made. To use information 
       criteria to decide whether to change the model use \code{\link{chooseSpatialModelOnIC.asrtests}}.
       
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
                corr.funcs = c("ar1", "ar1"), corr.orders = c(0, 0), 
                row.corrFitfirst = TRUE, allow.corrsJointFit = TRUE, 
                dropFixed = NULL, dropRandom = NULL, 
                nsegs = NULL, nestorder = c(1,1), 
                degree = c(3,3), difforder = c(2,2), 
                usRandLinCoeffs = TRUE, 
                rotateX = FALSE, ngridangles = NULL, 
                which.rotacriterion = "AIC", nrotacores = 1, 
                asreml.option = "grp", tpps4mbf.obj = NULL, 
                allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                checkboundaryonly = FALSE, update = FALSE, 
                maxit = 30, IClikelihood = "full", which.IC = "AIC", ...)
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
               models are to be fitted. Note that, for other terms that involve 
               \code{sections} in the random formula, there should be separate terms for 
               each level of \code{sections}. For example, in a blocked experiment 
               involving multiple sites, there should be the sum of separate terms for 
               the Blocks at each Site i.e. a formula that contains terms like 
               at(Site, i):Block for each site and these are separated by \code{'+'}. 
               Otherwise, the combined term (e.g. Site:Block) will impact on the 
               fitting of the local spatial models for the different Sites. Similarly, 
               a separate residual variance for each of the \code{sections} should be 
               fitted, unless there is a need to fit a different variance structure to 
               the residual, e.g. heterogeneous residual variances depending on 
               treatments. Separate residual variances for \code{sections} can be 
               achieved using the \code{asreml} functions \code{dsum} or \code{idh}. 
               Because, unlike random terms, terms for residual variances are not 
               removed from the model, compound residual terms can be used to include 
               them in the model, e.g. terms with \code{idh} or \code{dsum} with multiple
               \code{levels} in the \code{list} or leaving \code{levels} out altogether. 
               In addition to allowing the independent fitting of models to the 
               \code{sections}, separate residual variance terms allows a nugget variance 
               to be fitted in a correlation model for each of the \code{sections}.}
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
               specify \code{""} for that dimension. If the correlation model is  
               \code{corb}, the values of \code{corr.orders} are used for its order 
               argument (\code{b}).}
 \item{corr.orders}{A \code{numeric} of length two that specifies the order argument 
               (\code{b}) values for the row and column dimensions of a two-dimensional 
               separable spatial correlation model when \code{spatial.model} is 
               \code{corr} and the \code{corr.funcs} for a dimension is \code{corb}, 
               the \code{asreml} banded correlation model. If one of the dimensions 
               does not involve an order argument, set the value of \code{corr.orders} 
               for that dimension to zero. For a dimension for which the 
               \code{corr.funcs} is \code{corb} and \code{corr.orders} is zero, a 
               model with a single band, the correlation between immediate neighbours, 
               will be fitted and then further bands, up to a maximum of 10 bands, will 
               be added until the addition of an extra band does not reduce the 
               information criterion nominated using \code{which.IC}. Note that the 
               two-dimensional spatial model is fitted as a random term.}
 \item{row.corrFitfirst}{A \code{\link{logical}}. If \code{TRUE} then, in fitting the 
               model for \code{spatial.model} set to \code{corr}, the row correlation  
               or variance function is fitted first, followed by the addition of the 
               column correlation or variance function. If \code{FALSE}, the order 
               of fitting is reversed.}
 \item{allow.corrsJointFit}{A \code{\link{logical}} which, if \code{TRUE}, will allow 
               the simultaneous fitting of correlation functions for the two dimensions 
               of the grid when separate fits have failed to fit any correlation 
               functions. This argument is available for when a joint fit 
               hangs the system.}
 \item{dropFixed}{A single \code{character} string or a \code{character} vector of strings 
               with an element for each level of \code{sections} in the same order as the 
               \code{sections} levels. Each string, which if it is not \code{NA} and after 
               the addition of \code{". ~ . -"} and conversion to a formula that is then 
               expanded, specifies the  sum of a set of terms to be dropped from the fixed 
               formula in fitting splines (\code{TPPS} and \code{TPNCSS}). The result is 
               that the fitted model supplied in the \code{asrtests.obj}, that includes these 
               terms, will be compared with the fitted model that has had them removed and 
               a spatial model added.

               An element that is \code{NA} indicates that no term pertaining to the 
               corresponding \code{sections} level is to be removed. If \code{sections} 
               is not \code{NULL} and a single \code{character} string has been supplied, 
               the terms specified by the string are taken to be terms that are 
               independent of the \code{sections} and will be removed when adding the 
               spatial model for the first \code{sections}. 
 
               The terms must match those in the \code{wald.tab} component of the 
               \code{asrtests.obj}. The fixed terms will be reordered so that 
               single-variable terms come first, followed by two-variable terms and 
               so on. Note also that multiple terms specified using a single 
               \code{asreml::at} function can only be dropped as a whole. If the term 
               was specified using an \code{asreml::at} function with a single level, 
               then it can be removed and either the level itself or its 
               \code{\link{numeric}} position in the levels returned by the 
               \code{\link{levels}} function can be specified.}
 \item{dropRandom}{A single \code{character} string or a \code{character} vector of strings 
               with an element for each level of \code{sections} in the same order as the 
               \code{sections} levels. Each string, which if it is not \code{NA} and after 
               the addition of \code{" ~ . -"} and conversion to a formula that is then 
               expanded, specifies the  sum of a set of terms to be dropped from the random 
               formula in fitting splines (\code{TPPS} and \code{TPNCSS}). The result is 
               that the fitted model supplied in the \code{asrtests.obj}, that includes 
               these terms, will be compared with the fitted model that has had them removed 
               and a spatial model added.

               An element that is \code{NA} indicates that no term pertaining to the 
               corresponding \code{sections} level is to be removed. If \code{sections} 
               is not \code{NULL} and a single \code{character} string has been supplied, 
               the terms specified by the string are taken to be terms that are 
               independent of the \code{sections} and will be removed when adding the 
               spatial model for the first \code{sections}. 
 
               The terms must match those in the \code{vparameters} component of the 
               \code{asreml.obj} component in the \code{asrtests.obj}. Note also that 
               multiple terms specified using a single \code{asreml::at} function 
               can only be dropped as a whole. If the term was specified using 
               an \code{asreml::at} function with a single level, then it can be 
               removed and either the level itself or its \code{\link{numeric}} 
               position in the levels returned by the \code{\link{levels}} function 
               can be specified.}
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
              (\code{TPPS}).}
 \item{difforder}{A \code{\link{numeric}} of length 2. The order of differencing for 
              column and row dimensions, respectively, in fitting a P-spline 
              (\code{TPPS}).}
 \item{usRandLinCoeffs}{A \code{\link{logical}} which, if \code{TRUE}, will attempt to  
               fit an unstructured variance model to the constant and linear terms in 
               the interactions for constant and linear terms in one grid dimension 
               interacting with smooth terms in the second grid dimension. The 
               unstructured variance model can only be fitted if both the constant and 
               linear interaction terms have been retained in the fitted model. 
               This argument can be used to omit the attempt to fit an unstructured 
               variance model when the attempt results in a system error.}
 \item{rotateX}{A \code{\link{logical}} indicating whether to rotate the eigenvectors of the 
              penalty matrix, as described by Piepho, Boer and Williams (2022), when 
              fitting a P-spline (\code{TPPS}). Setting \code{rotateX} to \code{TRUE} 
              results in a search for an optimized rotation under a model that omits 
              the random spline interaction terms. If \code{ngridangles} is set to 
              \code{NULL}, the optimal rotation us found using an optimizer 
              (\code{nloptr::bobyqa}). Otherwise, the optimal rotation is found by 
              exploring the fit over a  two-dimensional grid of rotation angle pairs. 
              The optimization seeks to optimize the criterion nominated in 
              \code{which.rotacriterion}. Rotation of the eigenvectors is only relevant 
              for \code{difforder} values greater than 1 and has only been implemented 
              for \code{difforder} equal to 2.}
 \item{ngridangles}{A \code{\link{numeric}} of length 2. If \code{NULL} (the default), 
              the optimal pair of angles for rotating the eigenvectors of the penalty 
              matrix of a P-spline (\code{TPPS}) will be determined using a nonlinear 
              optimizer (\code{nloptr::bobyqa}). Otherwise, its two values specify the 
              numbers of angles between 0 and 90 degrees for each of the row and column 
              dimensions to be used in determining the optimal pair of angles. Specifying 
              factors of 90 will result in integer-valued angles. The number of grid 
              points, and hence re-analyses will be the product of the values of 
              (\code{ngridangles} + 1).}
 \item{which.rotacriterion}{A single \code{character} string nominating which of the 
              criteria, out of the \code{deviance}, the \code{likelihood}, the 
              \code{AIC} and the \code{BIC}, is to be used in determining the optimal 
              rotation of the eigenvectors of the penalty matrix. The \code{deviance} 
              uses the REML value computed by \code{asreml}; the other criteria use 
              the full likelihood, evaluated using the REML estimates, that is computed 
              by \code{\link{infoCriteria.asreml}}.}
 \item{nrotacores}{A \code{\link{numeric}} specifying the number of cores to deploy for 
              running the analyses required to search the two-dimensional grid of rotation 
              angles when \code{rotateX} is \code{TRUE}. Parallel processing has been 
              implemented for analyzing, for each column angle, the set of angles to be 
              investigated for the row dimension.  The default value of one means that 
              parallel processing will not be used. The value chosen for \code{nrotacores} 
              needs to balanced against the other processes that are using parallel 
              processing at the same time.}
 \item{asreml.option}{A single \code{character} string specifying whether the \code{grp} or 
              \code{mbf} methods are to be used to supply externally formed covariate 
              matrices to \code{asreml} when fitting a P-spline (\code{TPPS}). Compared 
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
              prior to calling \code{addSpatialModel.asrtests}. If \code{tpps4mbf.obj} 
              is \code{NULL}, \cr
              \code{\link{makeTPPSplineMats.data.frame}} will be called 
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
 \item{maxit}{A \code{\link{numeric}} specifying the maximum number of iterations that 
              \code{asreml} should perform in fitting a model.}
 \item{IClikelihood}{A \code{character} that controls both the occurrence and the type 
                    of likelihood for information criterion in the \code{test.summary} 
                    of the new \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML}, then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full}, 
                    then the AIC and BIC based on the full likelihood, evaluated 
                    using REML estimates, are included. 
                    (See also \code{\link{infoCriteria.asreml}}.)}
 \item{which.IC}{A \code{character} specifying the information criterion to be used in 
              selecting the best model. Possible values are \code{AIC} and \code{BIC}. 
              The value of the criterion for supplied model must exceed that for changed 
              model for the changed model to be returned. (For choosing the rotation angle 
              of the eigenvectors of the penalty matrix, see \code{which.rotacriterion}.}
\item{\dots}{Further arguments passed to \code{\link{changeModelOnIC.asrtests}}, \code{\link{newfit.asreml}}, \code{asreml} and 
           \code{tpsmmb}.}
}
\details{
The model to which the spatial models is to be added is supplied in the \code{asrtests.obj}. It should not include terms that will be included in the local spatial model. All spatial model terms are fitted as fixed or random. Consequently, the residual model does not have to be iid. The improvement in the fit resulting from the addition of a spatial model to the supplied model is evaluated. Note that the data must be in the order that corresponds to the \code{residual} argument with a variable to the right of another variable changes levels in the data frame faster than those of the other variable e.g. \code{Row:Column} implies that all levels for \code{Column} in consecutive rows of the \code{data.frame} with a single \code{Row} level. 

For the \code{corr} spatial model, the default model is an autocorrelation model of order one (\code{ar1}) for each dimension. However, any of the single dimension correlation/variance models from \code{asreml} can be specified for each dimension, as can no correlation model for a dimension; the models for the two dimensions can differ. Using a forward selection procedure, a series of models are tried, without removing boundary or singular terms, beginning with the addition of row correlation and followed by the addition of column correlation or, if the \code{row.corrFitfirst} is set to \code{FALSE}, the reverse order. If the fitting of the first-fitted correlation did not result in a model change because the fitting did not converge or correlations were fixed, but the fit of the second correlation was successful, then adding the first correlation will be retried. If one of the metric correlation functions is specified (e.g. \code{exp}), then the \code{row.covar} or \code{col.covar} will be used in the spatial model. However, because the correlations are fitted separately for the two dimensions, the \code{row.factor} and \code{col.factor} are needed for all models and is used for a dimension that does not involve a correlation/variance function for the fit being performed. Also, the correlation models are fitted as \code{random} terms and so the correlation model will include a variance parameter for the grid even when \code{ar1} is used to specify the correlation model, i.e. the model fitted is a variance model and there is no difference between \code{ar1} and \code{ar1v} in fitting the model. The variance parameter for this term represents the spatial variance and the fit necessarily includes a nugget term, this being the residual variance. If any correlation is retained in the model, for a section if \code{sections} is not \code{NULL}, then the need for a nugget term is assessed by fixing the corresponding residual variance to one, unless there are multiple residual variances and these are not related to the \code{sections}. Once the fitting of the correlation model has been completed, the \code{rmboundary} function will be executed with the \code{checkboundaryonly} value supplied in the \code{addSpatialModel.asrtests} call. Finally, checking for bound and singular random terms associated with the correlation model and residual terms will be carried out when there are correlation terms in the model and \code{checkboundaryonly} has been set to \code{FALSE}; as many as possible will be removed from the fitted model, in some cases by fixing variance terms to one.

%This paragraph differs from the others in not mentioning that changes will be made provided they improve the fit.
The tensor-product natural-cubic-smoothing-spline (\code{TPNCSS}) spatial model is as described by Verbyla et al. (2018), the tensor-product penalized-cubic-spline (\code{TPPSC2}) model with second-order differencing of the penalty is similar to that described by Rodriguez-Alvarez et al. (2018), and the tensor-product, first-difference-penalty, linear spline (\code{TPPSL1}) model is amongst those described by Piepho, Boer and Williams (2022).  The fixed terms for the spline models are \code{row.covar + col.covar + row.covar:col.covar} and the random terms are \code{spl(row.covar) + spl(col.covar) + dev(row.covar) + dev(col.covar) + spl(row.covar):col.covar + row.covar:spl(col.covar) + spl(row.covar):spl(col.covar)}, except that \code{spl(row.covar) + spl(col.covar)} is replaced with \code{spl(row.covar):int(col.covar) + int(row.covar):spl(col.covar)} in the \code{TPPSC2} model, where \code{int(.)} indicates an intercept or constant value specific to its argument. For \code{TPPSL1} models, the terms \code{spl(row.covar):col.covar + row.covar:spl(col.covar)} are omitted,  The supplied model should not include any of these terms. However, any fixed or random main-effect Row or Column term that has been included as an initial model for comparison with a spatial model can be removed prior to fitting the spatial model using \code{dropFixed} or \code{dropRandom}.  For the P-spline models with second-order differencing, the model matrices used to fit the pairs of random terms (i) \code{spl(row.covar):int(col.covar)} and  \code{spl(row.covar):col.covar} and (ii) \code{int(row.covar):spl(col.covar)} and \code{row.covar:spl(col.covar)} are transformed using the spectral decomposition of their penalty matrices. An unstructured variance model is tried for each of these pairs. For \code{TPPSC2}, it is also possible to optimize the rotation of the null-space eigenvectors of the penalty matrix for each of these random-term pairs (for more information see Piepho, Boer and Williams, 2022). The optimization is achieved either using an optimizer or takes the form of a search over a grid of rotation angles for a reduced model; the fit of the full model with rotation using the optimal rotation angles will be returned.

The \code{TPPCS} and \code{TPP1LS} models are fitted using functions from the \code{R} package \code{TPSbits} authored by Sue Welham (2022). There are two methods for supplying the spline basis information produced by \code{tpsmmb} to \code{asreml}. The \code{grp} method adds it to the \code{data.frame} supplied in the \code{data} argument of the \code{asreml} call. The \code{mbf} method creates smaller \code{\link{data.frame}s} with the spline basis information in the same environment as the internal function that calls the spline-fitting function. If it is desired to use in a later session, an \code{asreml} function, or \code{asrtests} function that calls \code{asreml}, (e.g. \code{predict.asreml}, \code{\link{predictPlus.asreml}}, or \code{\link{changeTerms.asrtests}}) on an \code{asreml.object} created using \code{mbf} terms, then the \code{mbf} \code{\link{data.frame}s} will need to be recreated using \code{\link{makeTPPSplineMats.data.frame}} in the new session, supplying, if there has been rotation of the penalty matrix eigenvectors, the \code{theta} values that are returned as the attribute \code{theta.opt} of the \code{asreml.obj}. 

All models utlize the function \code{\link{changeTerms.asrtests}} to fit the spatial model. Arguments from \code{tpsmmb} and \code{\link{changeTerms.asrtests}} can be supplied in calls to \code{addSpatialModel.asrtests} and will be passed on to the relevant function through the ellipses argument (\ldots).

The data for experiment can be divided \code{sections} and the same spatial model fitted separately to each. The fit over all of the sections is assessed. For more detail see \code{sections} above.

Each combination of a row.coords and a col.coords does not have to specify a single observation; 
for example, to fit a local spatial model to the main units of a split-unit design, 
each combination would correspond to a main unit and all subunits of the main unit would 
have the same combination.
}

\value{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj}, 
       possibly with attribute \code{theta.opt}, 
       (ii) \code{wald.tab}, and (iii) \code{test.summary} for the model that includes the 
       spatial model, unless the spatial model fails to be fitted when \code{allow.unconverged} 
       and/or \code{allow.fixedcorrelation} is set to \code{FALSE}. If the 
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
         \code{\link{addSpatialModelOnIC.asrtests}}, \cr
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

#Create an asrtests object with a P-spline spatial variation model 
spatial.asrt <- addSpatialModel(current.asrt, spatial.model = "TPPS", 
                                row.covar = "cRow", col.covar = "cColumn",
                                dropRowterm = "Row", dropColterm = "Column",
                                asreml.option = "grp")
infoCriteria(current.asrt$asreml.obj)

#Create an asrtests object with a P-spline spatial variation model 
#that includes rotation of the eigenvectors of the penalty matrix
spatial.asrt <- addSpatialModel(current.asrt, spatial.model = "TPPS", 
                                row.covar = "cRow", col.covar = "cColumn",
                                dropRowterm = "Row", dropColterm = "Column",
                                rotateX = TRUE, 
                                which.rotacriterion = "dev", 
                                nrotacores = parallel::detectCores(), 
                                asreml.option = "mbf")
infoCriteria(current.asrt$asreml.obj)
}}
\keyword{asreml}
