\name{as.predictions.frame}
\alias{as.predictions.frame}
\title{Forms a \code{\link{predictions.frame}} from a \code{data.frame}, ensuring that the correct 
       columns are present.}
\description{Creates a \code{\link{predictions.frame}} from a \code{\link{data.frame}} by adding the 
             \code{class} \code{\link{predictions.frame}} to it, and renaming the columns containing 
             the  \code{predictions}, \code{se}, \code{est.status} and \code{error.intervals}.}

\usage{as.predictions.frame(data, classify = NULL, 
                     predictions = NULL, se = NULL, est.status = NULL, 
                     interval.type = NULL, interval.names = NULL)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing columns giving the variables that uniquely index 
                    the predicted values and columns with the predicted values, their standard errors 
                    and, optionally, their estimation status (\code{est.status}).}
 \item{classify}{A \code{\link{character}} string giving the variables that 
                    define the margins of the multiway table that was predicted. 
                    Multiway tables are specified by forming an interaction type 
                    term from the classifying variables, that is, separating the 
                    variable names with the \code{:} operator. For predicting the overall 
                    mean, the \code{classify} is set to "\code{(Intercept)}".}
 \item{predictions}{A \code{\link{character}} giving the name of the column in \code{data} that 
                    contains the predicted values. This column will be renamed to 
                    \code{predicted.value}.}
 \item{se}{A \code{\link{character}} giving the name of the column in \code{data} that 
                    contains the standard errors of the predicted values. This column will 
                    be renamed to \code{standard.error}.}
  \item{est.status}{A \code{\link{character}} giving the name of the column in \code{data} that 
                    contains the estimation status of the predicted values. It will have a value 
                    \code{Estimable} for predicted values that have been estimated and a value 
                    \code{Aliased} for predicted values that are \code{NA}. If a column named 
                    \code{est.status} is not present in \code{data} and \code{est.status} is 
                    \code{NULL}, a column \code{est.status} will be generated.}
 \item{interval.type}{A \code{\link{character}} specifying the type of \code{error.intervals} 
                    stored in \code{data} that require renaming. If \code{NULL}, \code{error.intervals} 
                    will not be renamed, even if they are present. Otherwise, \code{interval.type} 
                    should be set to one of "CI", "SE" or "halfLSD".}
 \item{interval.names}{A \code{\link{character}} specifying the column names of the \code{lower} and 
                    \code{upper} limits stored in \code{data} that are to be renamed. The 
                    \code{character} must be of length two, with the first element being the name of 
                    the 'lower' limit and the second element being the name of the 'upper' limit.}
}
\value{An S3-class \code{\link{predictions.frame}}.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{predictions.frame}},  
          \code{\link{is.predictions.frame}}, \code{\link{predictions.frame}}, \cr 
          \code{\link{validPredictionsFrame}}}
\examples{
  data(Oats.dat)
  
  ## Use asreml to get predictions and associated statistics

  \dontrun{
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  Var.pred <- asreml::predict.asreml(m1.asr, classify="Nitrogen:Variety", 
                                      sed=TRUE)
  if (getASRemlVersionLoaded(nchar = 1) == "3")
    Var.pred <- Var.pred$predictions
 #Form predictions.frame changing asreml-R4 names to the standard names, if these are present
 Var.preds <- as.predictions.frame(Var.pred$pvals, se = "std.error", 
                                    est.status = "status")
  }
  
  ## Use lmerTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
  }

  ## Check the class and validity of the alldiffs object
  if (exists("Var.preds"))
  {
    is.predictions.frame(Var.preds)
    validPredictionsFrame(Var.preds)
  }
}

\keyword{asreml}
\keyword{htest}