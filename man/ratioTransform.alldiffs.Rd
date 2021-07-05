\name{ratioTransform.alldiffs}
\alias{ratioTransform.alldiffs}
\alias{ratioTransform}
\title{Calculates the ratios of nominated pairs of predictions stored in an \code{\link{alldiffs.object}}.}
\description{Ratio predictions and error intervals are formed for two levels of a factor, 
             the \code{ratio.factor}. For each pair of a level of the \code{ratio.factor} 
             in \code{numerator.levels} with a level in \code{denominator.levels}, 
             the ratio predictions are formed from all combinations of the other factors 
             as the ratio of the two predictions for each combination, along with confidence 
             intervals for the ratio predictions computed using the Fieller (1954) method.

            The printing of the components produced is controlled by the 
             \code{tables} argument.
}
\usage{
\method{ratioTransform}{alldiffs}(alldiffs.obj, ratio.factor, 
               numerator.levels, denominator.levels, 
               method = "Fieller", alpha = 0.05,  
               response = NULL, response.title = NULL, 
               tables = "predictions", ...)
}
\arguments{
 \item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
 \item{ratio.factor}{A \code{\link{character}} string giving the name of the factor 
                 for whose levels the ratios are to be calculated.}
 \item{numerator.levels}{A \code{\link{character}} string containing the levels 
               of \code{ratio.factor} to be used as numerators of the ratio.}
 \item{denominator.levels}{A \code{\link{character}} string containing the levels 
               of \code{ratio.factor} to be used as denominators of the ratio.}
 \item{method}{A \code{\link{character}} string specifying the method to use in calculating 
               the ratios and their \code{error.intervals}. At present only \code{Fieller} 
               is available. For the \code{Fieller} \code{method}, ratios of predictions 
               are formed and confidence intervals formed for them using Fieller's (1954) theorem.}
 \item{alpha}{A \code{\link{numeric}} giving the significance level for LSDs or one minus 
              the confidence level for confidence intervals.}
 \item{response}{A \code{character} specifying the response variable for the 
                    predictions. It is stored as an attribute to the \code{\link{alldiffs.object}} .}
 \item{response.title}{A \code{character} specifying the title for the response variable 
                    for the predictions. It is stored as an attribute to the 
                    \code{\link{alldiffs.object}}.}
 \item{tables}{A \code{\link{character}} vector containing either \code{none} or \code{predictions}}
 \item{\dots}{further arguments passed to \code{\link{linTransform.alldiffs}}.}
}
\value{A \code{list} of \code{\link{predictions.frame}}s, each containing the ratio predictions 
       and their confidence limits for a combination of the \code{numerator.levels} with the 
       \code{denominator.levels}.

       If \code{sortFactor} attribute of the \code{\link{alldiffs.object}} is set and is not the 
       \code{ratio.factor}, the predictions and their backtransforms will be sorted using 
       the \code{sortOrder} attribute of the  \code{\link{alldiffs.object}}.
}
\author{Chris Brien}
\references{Fieller, E. C. (1954). Some Problems in Interval Estimation. \emph{Journal of the 
               Royal Statistical Society.Series B (Methodological)}, \bold{16}, 175-185.}


\seealso{\code{\link{pairdiffsTransform}}, \code{\link{linTransform}}, \code{\link{predictPlus.asreml}}, 
         \code{\link{as.alldiffs}}, \cr 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}},  \code{\link{subset.alldiffs}}, 
         \code{\link{allDifferences.data.frame}}, \cr 
         \code{\link{redoErrorIntervals.alldiffs}},  
         \code{\link{recalcLSD.alldiffs}},  \code{\link{predictPresent.asreml}}, \cr
         \code{\link{plotPredictions.data.frame}}, \cr
         \code{as.Date}, \code{predict.asreml}}
\examples{
#### Form the ratios and Fieller CIs for RGR Salinity
load(system.file("extdata", "testDiffs.rda", package = "asremlPlus", mustWork = TRUE))
Preds.ratio.RGR <- ratioTransform(diffs.RGR,
                                  ratio.factor = "Salinity", 
                                  numerator.levels = "Salt",
                                  denominator.levels = "Control")

#### Form the ratios and Fieller CIs for Nitrogen compared to no Nitrogen                                  
data("Oats.dat")
\dontrun{
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                 random=~Blocks/Wplots,
                 data=Oats.dat)
current.asrt <- as.asrtests(m1.asr)
wald.tab <-  current.asrt$wald.tab
Var.diffs <- predictPlus(m1.asr, classify="Nitrogen:Variety", pairwise = TRUE,
                         Vmatrix = TRUE, error.intervals = "halfLeast",
                         meanLSD.type = "factor", LSDby = "Variety",
                         wald.tab = wald.tab)
}

 ## Use lme4 and emmmeans to get predictions and associated statistics
if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
  ## Set up a wald.tab
  int <- as.data.frame(rbind(rep(NA,4)))
  rownames(int) <- "(Intercept)"
  wald.tab <- anova(m1.lmer, ddf = "Kenward", type = 1)[,3:6]
  names(wald.tab) <- names(int) <- c("Df", "denDF", "F.inc", "Pr")
  wald.tab <- rbind(int, wald.tab)
  #Get predictions
  Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
  Var.preds <- summary(Var.emm)
  ## Modify Var.preds to be compatible with a predictions.frame
  Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                    se = "SE", interval.type = "CI", 
                                    interval.names = c("lower.CL", "upper.CL"))
  Var.vcov <- vcov(Var.emm)
  Var.sed <- NULL
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  
  #Create alldiffs object
  Var.diffs <- as.alldiffs(predictions = Var.preds, 
                           sed = Var.sed, vcov = Var.vcov, 
                           classify = "Nitrogen:Variety", response = "Yield", tdf = den.df)
} 

if (exists("Var.diffs"))
  Preds.ratio.OatsN <- ratioTransform(alldiffs.obj = Var.diffs,
                                      ratio.factor = "Nitrogen", 
                                      numerator.levels = c("0.2","0.4","0.6"),
                                      denominator.levels = "0.2")
}  
 
\keyword{asreml}