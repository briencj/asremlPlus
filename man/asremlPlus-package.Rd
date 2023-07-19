\newcommand{\packageVer}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Version")}}
\newcommand{\packageDate}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Date")}}
\name{asremlPlus-package}
\alias{asremlPlus-package}
\alias{asremlPlus}
\docType{package}
\title{
\packageTitle{asremlPlus}
}
\description{\packageDescription{asremlPlus}

\bold{Version:} {\packageVer{asremlPlus}}

\bold{Date:} \packageDate{asremlPlus}
}

\section{Index}{
Note that many of the function below are S3 methods so that the suffix can be omitted.Of course, whether ot not the suffix is omitted, the object supplied to the first argument must be of the class specified by the suffix.  For example \code{\link{getFormulae.asreml}} is a \code{\link{getFormulae}} method for an \code{asreml.object} and so  \code{.asreml} can be omitted and the object supplied to the first argument must be of class \code{asreml}.
\tabular{ll}{
(i) Data \cr\cr

\code{\link{Oats.dat}} 
      \tab  Data for an experiment to investigate nitrogen response of \cr
      \tab 3 oats varieties. \cr
\code{\link{Wheat.dat}} 
      \tab  Data for an experiment to investigate 25 varieties of \cr
      \tab  wheat. \cr
\code{\link{WaterRunoff.dat}}
      \tab Data for an experiment to investigate the quality of \cr
      \tab water runoff over time\cr
\cr
(ii) Object manipulation \cr\cr

\code{\link{as.alldiffs}}
      \tab Forms an \code{\link{alldiffs.object}} from the supplied \cr
      \tab predictions, along with those statistics, associated with the \cr 
      \tab predictions and their pairwise differnces, that have been supplied. \cr
\code{\link{asrtests}}
      \tab Pseudonym for \code{\link{as.asrtests}}. \cr
\code{\link{as.asrtests}}
      \tab Forms an \code{\link{asrtests.object}} that stores (i) a fitted asreml object, \cr
      \tab (ii) a pseudo-anova table for the fixed terms and \cr
      \tab (iii) a history of changes and hypthesis testing \cr
      \tab used in obtaining the model. \cr
\code{\link{as.predictions.frame}}
      \tab Forms a \code{\link{predictions.frame}} from a data.frame, ensuring that \cr
      \tab the correct columns are present. \cr
\code{\link{facCombine.alldiffs}}
      \tab Combines several factors into one in the components of \cr
      \tab an \code{\link{alldiffs.object}}.  \cr
\code{\link{facRecast.alldiffs}}
      \tab Reorders and/or revises the factor levels using the order of old levels \cr
      \tab in \code{levels.order} and the new labels for the levels given in \code{newlabels}. \cr
\code{\link{facRename.alldiffs}}
      \tab Renames \code{\link{factor}}s in the \code{prediction} component of \cr
      \tab an \code{\link{alldiffs.object}}.    \cr
\code{\link{getFormulae.asreml}}
      \tab Gets the formulae from an asreml object. \cr
\code{\link{is.alldiffs}}
      \tab A single-line function that tests whether an object is \cr 
      \tab of class alldiffs.   \cr
\code{\link{is.asrtests}}
      \tab A single-line function that tests whether an object is \cr 
      \tab of class asrtests.   \cr
\code{\link{is.predictions.frame}}
      \tab A single-line function that tests whether an object is \cr 
      \tab of classes \code{\link{predictions.frame}} and \code{\link{data.frame}}.   \cr
\code{\link{makeTPPSplineMats.data.frame}}
      \tab Make the spline basis matrices and data needed to fit \cr 
      \tab Tensor Product P-Splines. \cr
\code{\link{print.alldiffs}}
      \tab Prints the values in an \code{\link{alldiffs.object}} in a nice format. \cr
\code{\link{print.asrtests}}
      \tab Prints the values in an \code{\link{asrtests.object}}. \cr
\code{\link{print.LSDdata}}
      \tab Prints the components of a list containing data on the \cr
      \tab LSDs for all pairwise differences of predictions.      \cr
\code{\link{print.predictions.frame}}
      \tab Prints the values in a \code{\link{predictions.frame}}, with or without title and heading. \cr
\code{\link{print.test.summary}}
      \tab Prints a data.frame containing a test.summary. \cr
\code{\link{print.wald.tab}}
      \tab Prints a data.frame containing a Wald or pseudoanova table. \cr
\code{\link{printFormulae.asreml}}
      \tab Prints the formulae from an asreml object. \cr
\code{\link{sort.alldiffs}}
      \tab Sorts the components of an \code{\link{alldiffs.object}} according to \cr
      \tab the predicted values associated with a factor. \cr
\code{\link{subset.alldiffs}}
       \tab Subsets the components in an \code{\link{alldiffs.object}} according \cr 
       \tab to the supplied condition. \cr
\code{\link{subset.list}}
       \tab Forms a \code{list} that contains a subset of the components of \cr
       \tab the supplied \code{\link{list}}.     \cr
\code{\link{validAlldiffs}}
       \tab Checks that an object is a valid \code{\link{alldiffs.object}}. \cr
\code{\link{validAsrtests}}
       \tab Checks that an object is a valid \code{\link{asrtests.object}}. \cr
\code{\link{validPredictionsFrame}}
       \tab Checks that an object is a valid \code{\link{predictions.frame}}. \cr
\cr
(iii) Model modification \cr\cr

\code{\link{addSpatialModel.asrtests}}
      \tab Adds, to a supplied model, a spatial model that accounts for \cr
      \tab local spatial variation.      \cr
\code{\link{addSpatialModelOnIC.asrtests}}
      \tab Uses information criteria to decide whether to add a spatial \cr
      \tab model to account for local spatial variation. \cr
\code{\link{changeTerms.asrtests}}
      \tab Adds and drops terms from one or both of the fixed or random \cr 
      \tab model, replaces the residual (rcov) model with a new model \cr 
      \tab and changes bounds or intial values of terms.  \cr
\code{\link{iterate.asrtests}} 
      \tab Subject the fitted asreml.obj stored in an asrtests.object \cr
      \tab  to further iterations of the fitting process. \cr
\code{\link{newfit.asreml}} 
      \tab  Refits an asreml model with modified model formula \cr
      \tab  using either a call to 'update.asreml' or a direct \cr
      \tab  call to 'asreml'. \cr
\code{\link{reparamSigDevn.asrtests}}
       \tab Reparamterizes each random (deviations) term \cr
       \tab involving 'devn.fac' to a fixed term and ensures \cr
       \tab that the same term, with 'trend.num' replacing \cr
       \tab 'devn.fac', is included if any other term with \cr
       \tab 'trend.num' is included in 'terms'. \cr
\code{\link{rmboundary.asrtests}} 
       \tab Removes any boundary or singular variance components \cr
       \tab from the fit stored in 'asreml.obj' and records their \cr
       \tab removal in an \code{\link{asrtests.object}}. \cr
\code{\link{setvarianceterms.call}} 
       \tab Allows the setting of bounds and initial values \cr
       \tab for terms in the 'random' and 'residual' arguments of an \cr
       \tab 'asreml' call. \cr
\cr
(iv) Model selection \cr\cr

\code{\link{addto.test.summary}}
       \tab Adds a row to a \code{test.summary} \code{data.frame}. \cr
\code{\link{changeModelOnIC.asrtests}}
       \tab Uses information criteria to decide whether to change an \cr
       \tab already fitted model. \cr
\code{\link{chooseModel.asrtests}}
       \tab Determines and records the set of significant terms using an \cr
       \tab \code{\link{asrtests.object}}, taking into account the hierarchy \cr 
       \tab or marginality relations of the terms..\cr
\code{\link{chooseModel.data.frame}}
       \tab Determines the set of significant terms from results stored \cr
       \tab in a \code{data.frame}, taking into account the marginality  \cr 
       \tab relations of terms and recording the tests used in a \cr
       \tab \code{data.frame}. \cr
\code{\link{chooseSpatialModelOnIC.asrtests}}
       \tab Uses information criteria to choose the best fitting \cr 
       \tab spatial model for accounting for local spatial variation. \cr
\code{\link{getTestPvalue.asrtests}}
       \tab Gets the p-value for a test recorded in the test.summary \cr
       \tab data.frame of an \code{\link{asrtests.object}}.  \cr
\code{\link{infoCriteria.asreml}}
       \tab Computes AIC and BIC for models. \cr
\code{\link{infoCriteria.list}}
       \tab Computes AIC and BIC for models. \cr
\code{\link{recalcWaldTab.asrtests}}
       \tab Recalculates the denDF, F.inc and P values for a table \cr
       \tab of Wald test statistics obtained using 'wald.asreml'. \cr
\code{\link{REMLRT.asreml}} 
       \tab Performs a REML ratio test. \cr
\code{\link{bootREMLRT.asreml}} 
       \tab Performs a REML ratio test using the parametric \cr
       \tab bootstrap. \cr
\code{\link{testranfix.asrtests}}
       \tab Tests for a single fixed or random term in model \cr
       \tab fitted using 'asreml' and records the result in an \cr
      \tab  \code{\link{asrtests.object}}. \cr
\code{\link{testresidual.asrtests}}
       \tab Fits a new residual formula using 'asreml', tests \cr
       \tab whether the change is significant and records the \cr
       \tab result in an \code{\link{asrtests.object}}. \cr
\code{\link{testswapran.asrtests}}
       \tab Tests, using a REMLRT, the significance of the difference \cr
       \tab between the current random model and one in which oldterms \cr
       \tab are dropped and newterms are added. The result is recorded \cr
       \tab in an \code{\link{asrtests.object}}.\cr
\cr
(v) Model diagnostics and simulation\cr\cr

\code{\link{plotVariofaces}}
       \tab Plots empirical variogram faces, including envelopes, \cr
       \tab from supplied residuals as described by Stefanova, Smith \cr
       \tab & Cullis (2009). \cr
\code{\link{variofaces.asreml}}
       \tab Calculates and plots empirical variogram faces, including \cr
       \tab envelopes, as described by Stefanova, Smith & Cullis (2009). \cr
\code{\link{estimateV.asreml}}
       \tab Forms the estimated variance, random or residual matrix for \cr
       \tab the observations from the variance parameter estimates. \cr
\code{\link{simulate.asreml}}
       \tab Produce sets of simulated data from a multivariate normal \cr
       \tab distribtion and save quantites related to the simulated data. \cr
\cr
(vi) Prediction production and presentation \cr\cr

\code{\link{addBacktransforms.alldiffs}}
       \tab Adds or recalculates the backtransforms component of an \cr 
       \tab \code{\link{alldiffs.object}}. \cr
\code{\link{allDifferences.data.frame}}
       \tab Using supplied predictions and standard errors of pairwise \cr
       \tab differences or the variance matrix of predictions, forms \cr
       \tab all pairwise differences between the set of predictions, and \cr
       \tab p-values for the differences. \cr
\code{\link{exploreLSDs}}
       \tab Explores the computed LSD values for pairwise differences \cr
       \tab between predictions.    \cr
\code{\link{linTransform.alldiffs}} 
       \tab Calculates a linear transformation of the \cr
       \tab predictions stored in an \code{\link{alldiffs.object}}. \cr
\code{\link{pairdiffsTransform.alldiffs}} 
       \tab Calculates the differences between nominated pairs of \cr 
       \tab predictions stored in an \code{\link{alldiffs.object}} \cr
\code{\link{pickLSDstatistics}}
       \tab Pick LSDstatistics whose values minimize the number of \cr 
       \tab errors in pairwise comparisons of predictions. \cr
\code{\link{plotLSDerrors.data.frame}}
       \tab Plots a map of the supplied errors that occur in using the \cr
       \tab computedLSD values for pairwise differences between \cr
       \tab predictions. \cr
\code{\link{plotLSDerrors.alldiffs}}
       \tab Plots a map of the errors that occur in using the computed \cr
       \tab LSD values for pairwise differences between predictions. \cr
\code{\link{plotLSDs.data.frame}}
       \tab Plots a heat map of computed LSD values for pairwise \cr
       \tab differences between predictions. \cr
\code{\link{plotLSDs.alldiffs}}
       \tab Plots a heat map of computed LSD values for pairwise \cr
       \tab differences between predictions. \cr
\code{\link{plotPredictions.data.frame}}
       \tab Plots the predictions for a term, possibly with \cr
       \tab error bars. \cr
\code{\link{plotPvalues.alldiffs}}
       \tab Plots the p-values in the p.differences components \cr
       \tab of an \code{\link{alldiffs.object}} as a heat map. \cr
\code{\link{plotPvalues.data.frame}}
       \tab Plots the p-values in data.frame as a heat map.\cr
\code{\link{predictPlus.asreml}}
       \tab Forms the predictions and associated statistics for \cr
       \tab a term, using an asreml object and a wald.tab and \cr
       \tab taking into account that a numeric vector \cr
       \tab and a factor having parallel values may occur in the \cr
       \tab model. It stores the results in an object of class \cr
       \tab 'alldifffs' and may print the results. It can be \cr
       \tab when there are not parallel values.  \cr
\code{\link{predictPresent.asreml}}
       \tab Forms the predictions for each of one or more terms \cr
       \tab and presents them in tables and/or graphs. \cr
\code{\link{ratioTransform.alldiffs}} 
       \tab Calculates the ratios of noinated pairs of predictions \cr 
       \tab stored in an \code{\link{alldiffs.object}}.   \cr
\code{\link{recalcLSD.alldiffs}} 
       \tab Adds or recalculates the \code{\link{LSD.frame}} that is a \cr 
       \tab component of an \code{\link{alldiffs.object}}.           \cr
\code{\link{redoErrorIntervals.alldiffs}} 
       \tab Adds or replaces the error intervals stored in the \cr
       \tab prediction component of an \code{\link{alldiffs.object}}. \cr
\code{\link{renewClassify.alldiffs}}
       \tab Renews the components in an \code{\link{alldiffs.object}} \cr 
       \tab according to a new classify. \cr
\code{\link{sort.alldiffs}}
       \tab Sorts the components in an \code{\link{alldiffs.object}} \cr 
       \tab according to the predicted values associated with a factor. \cr
\code{\link{subset.alldiffs}}
       \tab Subsets the components in an \code{\link{alldiffs.object}} according \cr 
       \tab to the supplied condition. \cr
\code{\link{sort.predictions.frame}}
       \tab Sorts a \code{\link{predictions.frame}} according to the \cr
       \tab predicted values. associated with a factor.\cr

(vii) Response transformation \cr\cr

\code{\link{angular}}
       \tab Applies the angular transformation to proportions. \cr
\code{\link{angular.mod}}
       \tab Applies the modified angular transformation to a \cr
       \tab vector of counts. \cr
\code{\link{powerTransform}}
       \tab Performs a combination of a linear and a power \cr
       \tab transformation on a variable. The transformed \cr
       \tab variable is stored in the 'data.frame data'. \cr
\cr
(viii) Miscellaneous \cr\cr

\code{\link{getASRemlVersionLoaded}}
       \tab Finds the version of asreml that is loaded and \cr
       \tab returns the initial characters in version. \cr
\code{\link{loadASRemlVersion}}
       \tab Ensures that a specific version of asreml is loaded.\cr
\code{\link{num.recode}}
       \tab Recodes the unique values of a vector using the values  \cr
       \tab in a new vector. \cr
\code{\link{permute.square}}  
       \tab          Permutes the rows and columns of a square matrix. \cr
\code{\link{permute.to.zero.lowertri}}
       \tab Permutes a square matrix until all the lower \cr
       \tab triangular elements are zero. \cr
}
The functions whose names end in 'alldiffs" utilize an \code{\link{alldiffs.object}} that stores:
   (i) a \code{\link{predictions.frame}}, being a data frame containing predicted values, variables indexing them and their standard errores and estimability status; 
       the lower and upper limits of error intervals will be included when these are requested,
   (ii) optionally, square matrices containing all pairwise differences, the standard errors and p-values of the differences, 
        and a \code{data.frame} containing LSD values and their summary statistics, 
   (iii) optionally, the variance matrix of the predictions, and 
   (iv) if the response was trandformed for analysis, a data frame with backtransforms of the predicted values.

The functions whose names end in 'asrtests', which are most of the model functions, utilize an \code{\link{asrtests.object}} that stores:  
  (i) the currently fitted model in \code{asreml.obj}, 
  (ii) the table of test statistics for the fixed effects in \code{wald.tab}, and 
  (iii) a data frame that contains a history of the changes made to the model in \code{test.summary}. 
}
\author{
\packageAuthor{asremlPlus}

Maintainer: \packageMaintainer{asremlPlus}
}
\references{Butler, D. G., Cullis, B. R., Gilmour, A. R., Gogel, B. J. and 
            Thompson, R. (2018). \emph{{ASReml-R} Reference Manual Version 4.}  
            VSN International Ltd, \url{https://asreml.kb.vsni.co.uk/}.
            
            Stefanova, K. T., Smith, A. B. & Cullis, B. R. (2009) Enhanced diagnostics for the 
            spatial analysis of field trials. \emph{Journal of Agricultural, Biological, 
            and Environmental Statistics}, \bold{14}, 392--410.}
\keyword{ package }
\seealso{\code{asreml}}
\examples{\dontrun{
## Analyse wheat dat using asreml and asremlPlus (see also the Wheat Vignette)
## Set up for analysis
library(dae)
library(asreml)
library(asremlPlus)
## use ?Wheat.dat for data set details
data(Wheat.dat)

# Fit initial model
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
summary(current.asr)

# Intialize a testing sequence by loading the current fit into an asrtests object
current.asrt <- as.asrtests(current.asr, NULL, NULL)

# Check for and remove any boundary terms
current.asrt <- rmboundary(current.asrt)

#Unbind Rep, Row and Column components and reload into an asrtests object
current.asr <- setvarianceterms(current.asr$call, 
                                terms = c("Rep", "Rep:Row", "Rep:Column"), 
                                bounds = "U")
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary(current.asrt)
summary(current.asrt$asreml.obj)$varcomp
print(current.asrt, which = "testsummary")
print(current.asrt, which = "pseudoanova")

# Check term for within Column pairs (a post hoc covariate)
current.asrt <- testranfix(current.asrt, "WithinColPairs", drop.fix.ns=TRUE)
  
# Test nugget term
current.asrt <- testranfix(current.asrt, "units", positive=TRUE)
  
# Test Row autocorrelation
current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
                             label="Row autocorrelation", simpler=TRUE)
  
# Test Col autocorrelation (depends on whether Row autocorrelation retained)
(p <- getTestPvalue(current.asrt, label = "Row autocorrelation"))
{ if (p <= 0.05)
    current.asrt <- testresidual(current.asrt, "~ ar1(Row):Column", 
                                 label="Col autocorrelation", 
                                 simpler=TRUE, update=FALSE)
  else
    current.asrt <- testresidual(current.asrt, "~ Row:Column", 
                                 label="Col autocorrelation", 
                                 simpler=TRUE, update=FALSE)
}

# Output the results
print(current.asrt, which = "test")
info <- infoCriteria(current.asrt$asreml.obj)
summary(current.asrt$asreml.obj)$varcomp

# Get current fitted asreml object and update to include standardized residuals
current.asr <- current.asrt$asreml.obj
current.asr <- update(current.asr, aom=TRUE)
Wheat.dat$res <- residuals(current.asr, type = "stdCond")
Wheat.dat$fit <- fitted(current.asr)

#### Do diagnostic checking  
# Do residuals-versus-fitted values plot
with(Wheat.dat, plot(fit, res))
  
#Produce variogram and variogram faces plot (Stefanaova et al, 2009)
plot.varioGram(varioGram(current.asr))
faces <- variofaces(current.asr, V=NULL, units="addtores", 
                    maxiter=50, update = FALSE)

#Get Variety predictions, sorted in increasing order for the predicted values, 
#and all pairwise prediction differences and p-values
Var.diffs <- predictPlus(classify = "Variety", 
                         asreml.obj=current.asr, 
                         error.intervals="halfLeast",
                         wald.tab=current.asrt$wald.tab, 
                         sortFactor = "Variety",
                         tables = "predictions")
print(Var.diffs, which = c("differences", "p.differences"))

# Plot the Variety predictions, with halfLSD intervals, and the p-values
plotPredictions(Var.diffs$predictions, 
                classify = "Variety", y = "predicted.value", 
                error.intervals = "half")
plotPvalues(Var.diffs)
}}
