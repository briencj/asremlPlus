\name{addBacktransforms.alldiffs}
\alias{addBacktransforms.alldiffs}
\alias{addBacktransforms}
\title{Adds or recalculates the backtransforms component of an \code{\link{alldiffs.object}}.}
\description{Given an \code{\link{alldiffs.object}}, adds or recalculate its \code{backtransforms} component. 
             The values of \code{transform.power}, \code{offset}, \code{scale} and \code{transform.function} 
             from the \code{backtransforms} component will be used, unless this component is \code{NULL} 
             when the values supplied in the call will be used.}
\usage{
\method{addBacktransforms}{alldiffs}(alldiffs.obj, 
                  transform.power = 1, offset = 0, scale = 1, 
                  transform.function =  "identity", ...)}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
 \item{transform.power}{A \code{\link{numeric}} specifying the power of a transformation, if 
              one has been applied to the response variable. Unless it is equal 
              to 1, the default, back-transforms of the predictions will be 
              obtained and presented in tables or graphs as appropriate. 
              The back-transformation raises the predictions to the power equal 
              to the reciprocal of \code{transform.power}, unless it equals 0 in 
              which case the exponential of the predictions is taken.} 
 \item{offset}{A \code{\link{numeric}} that has been added to each value of the 
              response after any scaling and before applying any power transformation.}
 \item{scale}{A \code{\link{numeric}} by which each value of the response has been multiplied 
              before adding any offset and applying any power transformation.}
 \item{transform.function}{A \code{\link{character}} giving the name of a function that 
              specifies the scale on which the predicted values are defined. This may be the 
              result of a transformation of the data using the function or the use of the 
              function as a link function in the fitting of a generalized linear (mixed) 
              model (GL(M)M). The possible \code{transform.function}s are 
              \code{identity}, \code{log}, \code{inverse}, \code{sqrt}, \code{logit}, 
              \code{probit}, and \code{cloglog}.  The \code{predicted.values} and 
              \code{error.intervals}, if not \code{StandardError} intervals, will be 
              back-transformed using the inverse function of the \code{transform.function}. 
              The \code{standard.error} column will be set to \code{NA}, unless (i) 
              \code{asreml} returns columns named \code{transformed.value} and 
              \code{approx.se}, as well as those called \code{predicted.values} and 
              \code{standard.error} (such as when a GLM is fitted) and 
              (ii) the values in \code{transformed.value} are equal to those obtained by 
              backtransforming the \code{predicted.value}s using the inverse function 
              of the \code{transform.function}. Then, the \code{approx.se} values will be 
              saved in the \code{standard.error} column of the \code{backtransforms} 
              component of the returned \code{alldiffs.obj}. Also, the 
              \code{transformed.value} and \code{approx.se} columns are removed from both 
              the \code{predictions} and \code{backtransforms} components of the 
              \code{alldiffs.obj}. Note that the values that end up in the \code{standard errors} 
              column are approximate for the backtransformed values and are not used in 
              calculating \code{error.intervals}.}
 \item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
}
\value{An \code{\link{alldiffs.object}} with components 
       \code{predictions}, \code{vcov}, \code{differences}, \code{p.differences}, 
       \code{sed}, \code{LSD} and \code{backtransforms}.

       The \code{backtransforms} component will have the attributes (i) \code{LSDtype}, 
       \code{LSDby} and \code{LSDstatistic} added from the \code{predictions} component and 
       (ii) \code{transform.power}, \code{offset}, \code{scale}, and \code{link}.
}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{as.alldiffs}}, \code{\link{sort.alldiffs}}, 
         \code{\link{subset.alldiffs}}, \code{\link{print.alldiffs}},  \cr 
         \code{\link{renewClassify.alldiffs}}, \code{\link{redoErrorIntervals.alldiffs}},  
         \code{\link{plotPredictions.data.frame}}, \cr 
         \code{\link{predictPlus.asreml}}, \code{\link{predictPresent.asreml}}}
\examples{
\donttest{
##Subset WaterRunoff data to reduce time to execute
data(WaterRunoff.dat)
tmp <- subset(WaterRunoff.dat, Date == "05-18" & Benches != "3")

##Use asreml to get predictions and associated statistics

\dontrun{
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = log.Turbidity ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= tmp)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
TS.diffs <- predictPlus(classify = "Sources:Type", 
                        asreml.obj = current.asr, 
                        wald.tab = current.asrt$wald.tab, 
                        present = c("Sources", "Type", "Species"))
}

##Use lmeTest and emmmeans to get predictions and associated statistics

if (requireNamespace("lmerTest", quietly = TRUE) && 
    requireNamespace("emmeans", quietly = TRUE))
{
  m1.lmer <- lmerTest::lmer(log.Turbidity ~ Benches + (Sources * (Type + Species)) + 
                              (1|Benches:MainPlots),
                            data=tmp)
  TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Species)
  TS.preds <- summary(TS.emm)
  den.df <- min(TS.preds$df, na.rm = TRUE)
  ## Modify TS.preds to be compatible with a predictions.frame
  TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
  
  ## Form an all.diffs object and check its validity
  TS.vcov <- vcov(TS.emm)
  TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Species", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

## Recalculate the back-transforms of the predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  TS.diffs <- addBacktransforms.alldiffs(TS.diffs, transform.power = 0)
}
}}
\keyword{asreml}