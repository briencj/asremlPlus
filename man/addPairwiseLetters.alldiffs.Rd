\name{addPairwiseLetters.alldiffs}
\alias{addPairwiseLetters.alldiffs}
\alias{addPairwiseLetters}
\title{Adds a column of letters to a predictions.frame, there being a unique letter for 
       each set of predictions that are not significantly different.}
\description{Given an \code{\link{alldiffs.object}} with a \code{p.differences} component, 
             identify the sets of predictions such that the predictions within a set are 
             not significantly different and assign a lowercase letter to these, and only 
             these, predictions.}
\usage{
\method{addPairwiseLetters}{alldiffs}(alldiffs.obj, within = NULL, 
                   alpha = 0.05, ...)}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
\item{within}{A \code{\link{character}} (vector) of variables names, being the names of the 
             \code{\link{factor}s} or \code{\link{numeric}s} in the \code{classify}; for each 
             combination of the values of the \code{\link{factor}s} and \code{\link{numeric}s},  
             the letters indicating the differences between the predictions within the 
             combination are added to the \code{\link{predictions.frame}} in the 
             \code{\link{alldiffs.object}}.}
\item{alpha}{A \code{\link{numeric}} specifying the significance level for t-test for judging 
             whether a pair of predictions is different.}
\item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
}
\value{A \code{\link{alldiffs.object}} whose \code{\link{predictions.frame}} has a column 
       named \code{pairwiseSignif} that contains the letters that specify the sets of 
       predictions that are not significantly different}

\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}}
\examples{
data(WaterRunoff.dat)

##Use asreml to get predictions and associated statistics

\dontrun{
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
TS.diffs <- predictPlus(classify = "Sources:Type", 
                        asreml.obj = current.asr, 
                        wald.tab = current.asrt$wald.tab, 
                        present = c("Sources", "Type", "Species"))
}

## Use lmeTest and emmmeans to get predictions and associated statistics

if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * (Type + Species)) + 
                              (1|Benches:MainPlots),
                            data=na.omit(WaterRunoff.dat))
  TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Type)
  TS.preds <- summary(TS.emm)
  den.df <- min(TS.preds$df, na.rm = TRUE)
  ## Modify TS.preds to be compatible with a predictions.frame
  TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
   
  ## Form an all.diffs object and check its validity
  TS.vcov <- vcov(TS.emm)
  TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Type", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

## Add significance letters for predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  ## Obtain the letters that indicate the significances with each Type
  TS.diffs <- addPairwiseLetters(TS.diffs, within = "Type")
  
  ## Plot the predictions with the letters
  library(ggplot2)
  ggplot(data = TS.diffs$predictions, 
         mapping = aes(x = Sources, y = predicted.value, 
                       colour = Sources, fill = Sources)) + 
    facet_grid(cols = vars(Type)) + 
    geom_bar(stat = "identity") +
    labs(y = "pH") + 
    geom_errorbar(aes(ymin=lower.Confidence.limit, 
                      ymax=upper.Confidence.limit,
                      width=0.5), colour = "black") +
    geom_text(aes(y = (upper.Confidence.limit + 0.025*upper.Confidence.limit), 
                  label = pairwiseSignif), colour = "black") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle=90))
}
}
\keyword{asreml}