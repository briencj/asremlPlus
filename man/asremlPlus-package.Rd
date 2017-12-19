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
\tabular{ll}{
(i) Data \cr\cr

\code{\link{Wheat.dat}} 
      \tab  Data for an experiment to investigate 25 varieties of \cr
      \tab  wheat. \cr
\code{\link{WaterRunoff.dat}}
      \tab Data for an experiment to investigate the quality of \cr
      \tab water runoff over time\cr
\cr
(ii) Object manipulation \cr\cr

\code{\link{as.alldiffs}}
      \tab Forms an object of S3-class 'alldiffs' that stores the \cr
      \tab predictions for a model fitted using asreml, along \cr
      \tab with statistics for all pairwise differences.\cr
\code{\link{asrtests}} 
      \tab Forms an object of S3-class 'asrtests' that stores (i) a \cr
      \tab fitted asreml object, (ii) a pseudo-anova table for the \cr
      \tab fixed and (iii) a history of changes and hypthesis testing \cr
      \tab used in obtaining the model. \cr
\code{\link{print.alldiffs}}
      \tab Prints the values in an 'alldiffs' object in a nice format. \cr
\code{\link{print.asrtests}}
      \tab Prints the values in an 'asrtests' object. \cr
\code{\link{sort.alldiffs}}
      \tab Sorts the components of an 'alldiffs' object according to \cr
      \tab the predicted values associated with a factor. \cr
\code{\link{subset.alldiffs}}
       \tab Subsets the components in an \code{alldiffs} object according \cr 
       \tab to the supplied condition. \cr
\cr
(iii) Model modification \cr\cr

\code{\link{changeTerms.asrtests}}
      \tab  Adds and drops the specified sets of terms from one \cr
      \tab or both of the fixed  or random model and/or replaces the \cr
      \tab residual (rcov) model with a new model. \cr
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
       \tab removal in a data.frame. \cr
\code{\link{setvarianceterms.call}} 
       \tab Allows the seting of bounds and initial values \cr
       \tab for terms in the 'random' and 'residual' arguments of an \cr
       \tab 'asreml' call. \cr
\cr
(iv) Model testing \cr\cr

\code{\link{chooseModel.asrtests}}
       \tab Determines the set of significant terms taking into \cr
       \tab account hierarchy or marginality relations and records \cr 
       \tab the tests performed in a data.frame. \cr
\code{\link{infoCriteria.asreml}}
       \tab        Computes AIC and BIC for a model. \cr
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
       \tab fitted using 'asreml' and records the result in a \cr
      \tab  data.frame. \cr
\code{\link{testresidual.asrtests}}
       \tab Fits a new residual formula using 'asreml', tests \cr
       \tab whether the change is significant and records the \cr
       \tab result in a data.frame. \cr
\code{\link{testswapran.asrtests}}
       \tab Tests, using a REMLRT, the significance of the difference \cr
       \tab between the current random model and one in which oldterms \cr
       \tab are dropped and newterms are added. The result is recorded \cr
       \tab in a data.frame.\cr
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

\code{\link{allDifferences.data.frame}}
       \tab Using supplied predictions, forms all pairwise \cr
       \tab differences between a set of predictions, their \cr
       \tab standard errors and p-values for the differences. \cr
\code{\link{plotPredictions.data.frame}}
       \tab Plots the predictions for a term, possibly with \cr
       \tab error bars. \cr
\code{\link{plotPvalues.alldiffs}}
       \tab Plots the p-values in the p.differences components \cr
       \tab of an \code{alldiffs} object as a heat map. \cr
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
\code{\link{recalcLSD.alldiffs}} 
       \tab Adds or recalculates the LSD component of an \cr 
       \tab \code{alldiffs} object.                           \cr
\code{\link{redoErrorIntervals.alldiffs}} 
       \tab Adds or replaces the error intervals stored in the \cr
       \tab prediction component of an alldiffs.object. \cr
\code{\link{sort.alldiffs}}
       \tab Sorts the components in an \code{alldiffs} object \cr 
       \tab according to the predicted values associated with a factor. \cr
\code{\link{subset.alldiffs}}
       \tab Subsets the components in an \code{alldiffs} object according \cr 
       \tab to the supplied condition. \cr
\cr
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

\code{\link{num.recode}}
       \tab Recodes the unique values of a vector using the values  \cr
       \tab in a new vector. \cr
\code{\link{permute.square}}  
       \tab          Permutes the rows and columns of a square matrix. \cr
\code{\link{permute.to.zero.lowertri}}
       \tab Permutes a square matrix until all the lower \cr
       \tab triangular elements are zero. \cr
}
The functions whose names end in 'asrtests', which are most of the model functions, utilize an \code{asrtests} object that stores:  
  (i) the currently fitted model in \code{asreml.obj}, 
  (ii) the table of test statistics for the fixed effects in \code{wald.tab}, and 
  (iii) a data frame that contains a history of the changes made to the model in \code{test.summary}. 
}
\author{
\packageAuthor{asremlPlus}

Maintainer: \packageMaintainer{asremlPlus}
}
\references{Butler, D. G., Cullis, B. R., Gilmour, A. R., Gogel, B. J. and 
            Thompson, R. (2017). \emph{{ASReml-R} User Guide Version 4.}  
            VSN International Ltd, \url{http://www.vsni.co.uk/software/asreml/}.
            
            Stefanova, K. T., Smith, A. B. & Cullis, B. R. (2009) Enhanced diagnostics for the 
            spatial analysis of field trials. \emph{Journal of Agricultural, Biological, 
            and Environmental Statistics}, \bold{14}, 392--410.}
\keyword{ package }
\seealso{\code{asreml}}
\examples{\dontrun{
## Analyse wheat dat using asreml4 and asremlPlus
## Set up for analysis
library(dae)
library(asreml4)
library(asremlPlus)
## use ?Wheat.dat for data set details
data(Wheat.dat)

# Fit initial model
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
summary(current.asr)

# Load current fit into an asrtests object
current.asrt <- asrtests(current.asr, NULL, NULL)

# Check for and remove any boundary terms
current.asrt <- rmboundary(current.asrt)

# Check term for within Column pairs
current.asrt <- testranfix(current.asrt, "WithinColPairs", drop.fix.ns=TRUE)
  
# Test nugget term
current.asrt <- testranfix(current.asrt, "units", positive=TRUE)
  
# Test Row autocorrelation
current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
                             label="Row autocorrelation", simpler=TRUE)
  
# Test Col autocorrelation (depends on whether Row autocorrelation retained)
k <- match("Row autocorrelation", current.asrt$test.summary$terms)
p <- current.asrt$test.summary$p
{ if (p[k] <= 0.05)
    current.asrt <- testresidual(current.asrt, "~ ar1(Row):Column", 
                                 label="Col autocorrelation", simpler=TRUE,
                                 update=FALSE)
    else
      current.asrt <- testresidual(current.asrt, "~ Row:Column", 
                                   label="Col autocorrelation", simpler=TRUE,
                                   update=FALSE)
}
print(current.asrt)
info <- infoCriteria(current.asrt$asreml.obj)

# Get current fitted asreml object
current.asr <- current.asrt$asreml.obj
current.asr <- update(current.asr, aom=TRUE)
  
  
# Do residuals-versus-fitted values plot
plot(fitted(current.asr),residuals(current.asr))
  
#Produce variogram and variogram faces plot (Stefanaova et al, 2009)
plot.varioGram(varioGram(current.asr))
faces <- variofaces(current.asr, V=V, maxiter=50, units="addtores")

#Get Variety predictions, sorted in increasing order for the predicted values, 
#and all pairwise prediction differences and p-values. Plot the latter.
Var.diffs <- predictPlus(classify = "Variety", 
                         asreml.obj=current.asr, 
                         error.intervals="halfLeast",
                         wald.tab=current.asrt$wald.tab,
                         tables = "predictions", 
                         sortFactor = "Variety")
print(Var.diffs, which = c("differences", "p.differences"))
plotPvalues(Var.diffs)
}}
