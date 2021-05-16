\name{asremlPlus-deprecated}
\alias{asremlPlus-deprecated}
\alias{addrm.terms.asreml}
\alias{addrm.terms.asrtests}
\alias{alldiffs}
\alias{choose.model.asreml}
\alias{choose.model.asrtests}
\alias{facRecode}
\alias{facRecode.alldiffs}
\alias{info.crit}
\alias{info.crit.asreml}
\alias{newrcov.asrtests}
\alias{plotvariofaces.asreml}
\alias{predictiondiffs.asreml}
\alias{predictionplot.asreml}
\alias{predictparallel.asreml}
\alias{pred.present.asreml}
\alias{power.transform}
\alias{recalc.wald.tab.asreml}
\alias{recalc.wald.tab.asrtests}
\alias{reorderClassify}
\alias{reorderClassify.alldiffs}
\alias{reml.lrt}
\alias{reml.lrt.asreml}
\alias{rmboundary.asreml}
\alias{setvarianceterms.asreml}
\alias{sig.devn.reparam.asreml}
\alias{sig.devn.reparam.asrtests}
\alias{testranfix.asreml}
\alias{testrcov.asreml}
\alias{testrcov.asrtests}
\alias{testswapran.asreml}
\title{Deprecated Functions in the Package asremlPlus}
\description{These functions have been renamed and deprecated in \code{asremlPlus}: 
       \enumerate{
         \item addrm.terms.asreml and addrm.terms.asrtests -> \code{\link{changeTerms.asrtests}}, 
         \item alldiffs -> \code{\link{as.alldiffs}}, 
         \item choose.model.asreml and choose.model.asrtests -> \code{\link{chooseModel.asrtests}}, 
         \item facRecode and facRecode.alldiffs -> \code{\link{facRecast.alldiffs}},
         \item info.crit and info.crit.asreml -> \code{\link{infoCriteria.asreml}}, 
         \item newrcov.asrtests -> \code{\link{changeTerms.asrtests}}, 
         \item plotvariofaces.asreml -> \code{\link{plotVariofaces.data.frame}}, 
         \item power.transform -> \code{\link{powerTransform}}, 
         \item predictiondiffs.asreml -> \code{\link{allDifferences.data.frame}}, 
         \item predictionplot.asreml -> \code{\link{plotPredictions.data.frame}}, 
         \item predictparallel.asreml -> \code{\link{predictPlus.asreml}}, 
         \item pred.present.asreml -> \code{\link{predictPresent.asreml}}, 
         \item recalc.wald.tab.asreml and recalc.wald.tab.asrtests -> \code{\link{recalcWaldTab.asrtests}}, 
         \item reorderClassify and reorderClassify.alldiffs -> \code{\link{renewClassify.alldiffs}}, 
         \item reml.lrt and reml.lrt.asreml -> \code{\link{REMLRT.asreml}},
         \item rmboundary.asreml -> \code{\link{rmboundary.asrtests}}, 
         \item setvarianceterms.asreml  -> \code{\link{setvarianceterms.call}}, 
         \item sig.devn.reparam.asreml and sig.devn.reparam.asrtests -> \code{\link{reparamSigDevn.asrtests}}, 
         \item testranfix.asreml -> \code{\link{testranfix.asrtests}}, 
         \item testrcov.asreml and testrcov.asrtests -> \code{\link{testresidual.asrtests}}, 
         \item testswapran.asreml -> \code{\link{testswapran.asrtests}}}}
\usage{addrm.terms.asreml(...)
addrm.terms.asrtests(...)
alldiffs(...)
choose.model.asreml(...)
choose.model.asrtests(...)
facRecode(...)
facRecode.alldiffs(...)
info.crit(...)
info.crit.asreml(...)
newrcov.asrtests(...)
plotvariofaces.asreml(...)
power.transform(...)
predictiondiffs.asreml(...)
predictionplot.asreml(...)
predictparallel.asreml(...)
pred.present.asreml(...)
recalc.wald.tab.asreml(...)
recalc.wald.tab.asrtests(...)
reml.lrt(...)
reml.lrt.asreml(...)
\method{reorderClassify}{alldiffs}(...)
\method{rmboundary}{asreml}(...)
setvarianceterms.asreml(...)
sig.devn.reparam.asreml(...)
sig.devn.reparam.asrtests(...)
testranfix.asreml(...)
testrcov.asreml(...)
testrcov.asrtests(...)
\method{testswapran}{asreml}(...)
}
\arguments{
  \item{...}{absorbs arguments passed from the old functions of the style foo.bar().}
}
\author{Chris Brien}
