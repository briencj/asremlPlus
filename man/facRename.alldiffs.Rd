\name{facRename.alldiffs}
\alias{facRename.alldiffs}
\alias{facRename}
\title{Renames \code{\link{factor}}s in the \code{prediction} component of an 
       \code{\link{alldiffs.object}}.}
\description{Renames \code{\link{factor}}s in the \code{prediction} component of an 
       \code{\link{alldiffs.object}}. These changes are propogated to the other 
        components and attributes of the \code{\link{alldiffs.object}}.}
\usage{\method{facRename}{alldiffs}(object, factor.names, newnames,  ...)}
\arguments{
 \item{object}{An \code{\link{alldiffs.object}}.}
 \item{factor.names}{A \code{\link{character}} containing the names of the 
                     \code{\link{factor}}s in the \code{prediction} component 
                     of \code{object} that are to be renamed.}
 \item{newnames}{A \code{\link{character}} containing the new names of the 
                     \code{\link{factor}}s in the \code{prediction} component 
                     of \code{object}.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A modified \code{\link{alldiffs.object}}.}
\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, \cr
         \code{\link{facCombine.alldiffs}}, \code{\link{facRecast.alldiffs}}, 
         \code{\link{renewClassify.alldiffs}}; 
         \code{fac.recode} in package \pkg{dae}.}
\examples{
  data("Ladybird.dat")
  
  ## Use asreml to get predictions and associated statistics

  \dontrun{
  m1.asr <- asreml(logitP ~ Host*Cadavers*Ladybird, 
                   random = ~ Run,
                   data = Ladybird.dat)
  current.asrt <- as.asrtests(m1.asr)
  HCL.pred <- asreml::predict.asreml(m1.asr, classify="Host:Cadavers:Ladybird", 
                                     sed=TRUE)
  HCL.preds <- HCL.pred$pvals
  HCL.sed <- HCL.pred$sed
  HCL.vcov <- NULL
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Host:Cadavers:Ladybird", rownames(wald.tab)), "denDF"]
  }
  
  ## Use lmeTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(logitP ~ Host*Cadavers*Ladybird + (1|Run),
                              data=Ladybird.dat)
    HCL.emm <- emmeans::emmeans(m1.lmer, specs = ~ Host:Cadavers:Ladybird)
    HCL.preds <- summary(HCL.emm)
    den.df <- min(HCL.preds$df)
    ## Modify HCL.preds to be compatible with a predictions.frame
    HCL.preds <- as.predictions.frame(HCL.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    HCL.vcov <- vcov(HCL.emm)
    HCL.sed <- NULL
  }
  
  ## Use the predictions obtained with either asreml or lmerTest
  if (exists("HCL.preds"))
  {
    ## Form an all.diffs object
    HCL.diffs <- allDifferences(predictions = HCL.preds, 
                                classify = "Host:Cadavers:Ladybird", 
                                sed = HCL.sed, vcov = HCL.vcov, tdf = den.df)
    
    ## Check the class and validity of the alldiffs object
    is.alldiffs(HCL.diffs)
    validAlldiffs(HCL.diffs)

    ## Rename Cadavers 
    HCL.diffs <- facRename(HCL.diffs, factor.names = "Cadavers", newnames = "Cadaver.nos")
    
    ## Check the validity of HCL.diffs
    validAlldiffs(HCL.diffs)
  }
}
\keyword{asreml}
\keyword{factor}
\keyword{manip}