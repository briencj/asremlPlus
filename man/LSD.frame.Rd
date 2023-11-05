\name{LSD.frame}
\alias{LSD.frame}
%\alias{LSD.frame-class}
\title{Description of an LSD frame}
\description{A \code{data.frame} that stores  
             Least Significant differences (LSDs) for \code{predictions} for a fitted model.}

%             \code{\link{as.LSD.frame}} is function that converts a 
%             \code{\link{data.frame}} to an object of this class.
%
%             \code{\link{is.LSD.frame}} is the membership function for this class; it tests 
%             that an object has class \code{LSD.frame}.

%             \code{\link{validLSDFrame}} can be used to test the validity of a 
%             \code{LSD.frame}.
%}
\value{A \code{data.frame} that can be a component of an \code{\link{alldiffs.object}} and that 
       contains LSD values and statistics to be used in determining the significance of the 
       pairwise diffrences. In particular, they are used in calculating 
       \code{halfLeastSignificant} limits to be included in a \code{predictions.frame}.

       Exactly what an \code{LSD.frame} contains is 
       determined by the following arguments to functions that return an 
       \code{\link{alldiffs.object}}: \code{LSDtype}, \code{LSDby}, \code{LSDstatistic}, 
       \code{LSDaccuracy} and \code{LSDsupplied}. The \code{rownames} of the \code{LSD.frame} 
       indicate, for each of its rows, for what group of predictions the entries in the row were calculated, 
       this being controlled by the \code{LSDtype} and \code{LSDby} arguments. The values for 
       all of the LSD arguments are stored as attributes to the \code{\link{alldiffs.object}} and the 
       \code{predictions} and, if present \code{backtransforms}, components of the 
       \code{\link{alldiffs.object}}.
       
       An \code{LSD.frame} always has the eight columns \code{c}, \code{minimumLSD}, \code{meanLSD}, 
       \code{maximumLSD}, \code{assignedLSD}, \code{accuracyLSD}, \code{falsePos} and 
       \code{falseNeg}.
       \enumerate{
        \item \code{c}: This gives the number of pairwise coomparison of predictions for the combinations of 
             the factor levels given by the row name. If the row name is \code{overall} then it is for 
             all predictions.
        \item \code{minimumLSD, meanLSD, maximumLSD}: These are computed for either \code{overall}, 
             \code{factor.combinations}, \code{per.prediction} or \code{supplied} LSD values, as specified by the 
             \code{LSDtype} argument.  The \code{meanLSD} is calculated using the square root of the mean of 
             the variances of set of pairwise differences appropriate to the specific \code{LSDtype} argument.
             
             For \code{overall}, the mean, minimum and maximum of the LSDs for all pairwise 
             comparisons are computed.  
             
             If \code{factor.combinations} was specified for 
             \code{LSDtype} when the LSDs were being calculated, then the \code{\link{LSD.frame}} 
             contains a row for each combination of the values of the \code{\link{factor}s} and 
             \code{\link{numeric}s} specified by \code{LSDby}. The values in a row  are calculated 
             from the LSD values for the pairwise differences for each combination of the 
             \code{\link{factor}s} and \code{\link{numeric}s} values, unless there is only one 
             prediction for a combination, when notional LSDs are calculated that are based on the 
             standard error of the prediction multiplied by the square root of two. 
             
             For \code{per.prediction}, the minimum, mean and maximum LSD, based, for each prediction, 
             on the LSD values for all pairwise differences involving that prediction are computed.
             
             For \code{supplied}, the \code{LSD.frame} is set up based on the setting of \code{LSDby}: 
             a single row with name \code{overall} if \code{LSDby} is \code{NULL} or, if \code{LSDby} 
             is a vector of \code{\link{factor}} and \code{\link{numeric}} names, rows for each observed 
             combinations of the values of the named \code{\link{factor}s} and \code{\link{numeric}s}. 
             The \code{LSDsupplied} argument is used to provide the values to be stored in the column 
             \code{assignedLSD}.
         \item \code{assignedLSD}: The \code{assignedLSD} column contains the values that are assigned for 
             use in calculating \code{halfLeastSignificant} \code{error.intervals}. Its contents are 
             determined by \code{LSDstatistic} and \code{LSDsupplied} arguments. The 
             \code{LSDsupplied} argument allows the direct specification of values to be placed 
             in the \code{assignedLSD} column of the \code{\link{LSD.frame}}. The default is to use the 
             values in the \code{meanLSD} column.
         \item \code{LSDaccuracy}: The \code{LSDaccuracy} gives an indication of the proportion that the 
              correct LSD for a single \code{predicted.value} might deviate from its \code{assignedLSD} value. 
              The contents of the \code{accuracyLSD} column is controlled by the \code{LSDaccuracy} 
              argument.
         \item \code{falsePos} and \code{falseNeg}: These columns contain the number of false positives 
              and negatives if the \code{assignedLSD} value(s) is(are) used to determine the significance 
              of the pairwise predictions differences. Each LSD value in the \code{assignedLSD} 
              column is used to determine the significance of pairwise differences that involve 
              predictions for the combination of values given by the row name for the LSD value.
      }
       
       See \code{\link{recalcLSD.alldiffs}} for more information.

}
\examples{
  data(Oats.dat)
  
  ## Use asreml to get predictions and associated statistics

  \dontrun{
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  Var.diffs <- predictPlus(m1.asr, classify="Nitrogen:Variety", 
                          wald.tab = current.asrt$wald.tab, 
                          tables = "none")
  }
  
  ## Use lmerTest and emmmeans to get predictions and associated statistics
   if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    #Get predictions
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    ## Modify Var.preds to be compatible with a predictions.frame
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    Var.vcov <- vcov(Var.emm)
    Var.sed <- NULL

    #Set up an alldiffs object, which includes overall LSDs
    Var.diffs <- allDifferences(predictions = Var.preds, classify = "Variety:Nitrogen", 
                                     sed = Var.sed, vcov = Var.vcov, tdf = 45)
  }

  if (exists("Var.diffs"))
  {
    ## Use recalcLSD to get LSDs for within Variety differences
    Var.LSD.diffs <- recalcLSD(Var.diffs, 
                               LSDtype = "factor.combinations", LSDby = "Variety")
    print(Var.LSD.diffs$LSD)
  }
}

\author{Chris Brien}
\seealso{\code{\link{recalcLSD.alldiffs}}, \code{\link{redoErrorIntervals.alldiffs}}, 
          \code{\link{predictPresent.asreml}}, \cr 
          \code{\link{predictPlus.asreml}}}
\keyword{asreml}
\keyword{htest}