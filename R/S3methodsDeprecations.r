#S3 method generics

addSpatialModel <- function(asrtests.obj, ...) UseMethod("addSpatialModel")
addSpatialModelOnIC <- function(asrtests.obj, ...) UseMethod("addSpatialModelOnIC")
changeModelOnIC <- function(asrtests.obj, ...) UseMethod("changeModelOnIC")
changeTerms <- function(asrtests.obj, ...) UseMethod("changeTerms")
chooseModel <- function(object, ...) UseMethod("chooseModel")
chooseSpatialModelOnIC <- function(asrtests.obj, ...) UseMethod("chooseSpatialModelOnIC")
getTestEntry <- function(asrtests.obj, ...) UseMethod("getTestEntry")
getTestPvalue <- function(asrtests.obj, ...) UseMethod("getTestPvalue")
iterate <- function(asrtests.obj, ...) UseMethod("iterate")
recalcWaldTab <- function(asrtests.obj, ...) UseMethod("recalcWaldTab")
rmboundary <- function(asrtests.obj, ...) UseMethod("rmboundary")
reparamSigDevn <- function(asrtests.obj, ...) UseMethod("reparamSigDevn")
testranfix <- function(asrtests.obj, ...) UseMethod("testranfix")
testresidual <- function(asrtests.obj, ...) UseMethod("testresidual")
testswapran <- function(asrtests.obj, ...) UseMethod("testswapran")

bootREMLRT <- function(h0.asreml.obj, h1.asreml.obj, ...) UseMethod("bootREMLRT")
convAsremlobj <- function(asreml.obj, ...) UseMethod("convAsremlobj")
convEffectNames2DataFrame <- function(asreml.obj, ...) UseMethod("convEffectNames2DataFrame")
estimateV <- function(asreml.obj, ...) UseMethod("estimateV")
getFormulae <- function(asreml.obj, ...) UseMethod("getFormulae")
infoCriteria <- function(object, ...) UseMethod("infoCriteria")
isFixedCorrelOK <- function(asreml.obj, ...) UseMethod("isFixedCorrelOK")
newfit <- function(asreml.obj, ...) UseMethod("newfit")
predictPlus <- function(asreml.obj, ...) UseMethod("predictPlus")
predictPresent <- function(asreml.obj, ...) UseMethod("predictPresent")
printFormulae <- function(asreml.obj, ...) UseMethod("printFormulae")
R2adj <- function(asreml.obj, ...) UseMethod("R2adj")
REMLRT <- function(h0.asreml.obj, h1.asreml.obj, ...) UseMethod("REMLRT")
setmbfenv <- function(asreml.obj, ...) UseMethod("setmbfenv")
variofaces <- function(asreml.obj, ...) UseMethod("variofaces")

addBacktransforms <- function(alldiffs.obj, ...) UseMethod("addBacktransforms")
exploreLSDs <- function(alldiffs.obj, ...) UseMethod("exploreLSDs")
findLSDminerrors <- function(alldiffs.obj, ...) UseMethod("findLSDminerrors")
linTransform <- function(alldiffs.obj, ...) UseMethod("linTransform")
pairdiffsTransform <- function(alldiffs.obj, ...) UseMethod("pairdiffsTransform")
pickLSDstatistics <- function(alldiffs.obj, ...) UseMethod("pickLSDstatistics")
ratioTransform <- function(alldiffs.obj, ...) UseMethod("ratioTransform")
recalcLSD <- function(alldiffs.obj, ...) UseMethod("recalcLSD")
redoErrorIntervals <- function(alldiffs.obj, ...) UseMethod("redoErrorIntervals")
renewClassify <- function(alldiffs.obj, ...) UseMethod("renewClassify")

setvarianceterms <- function(call, ...) UseMethod("setvarianceterms")

allDifferences <- function(predictions, ...) UseMethod("allDifferences")
facCombine <- function(object, ...) UseMethod("facCombine")
facRecast <- function(object, ...) UseMethod("facRecast")
facRename <- function(object, ...) UseMethod("facRename")
makeTPPSplineMats <- function(data, ...) UseMethod("makeTPPSplineMats")
plotLSDerrors <- function(object, ...) UseMethod("plotLSDerrors")
plotLSDs <- function(object, ...) UseMethod("plotLSDs")
plotPvalues <- function(object, ...) UseMethod("plotPvalues")
plotVariofaces <- function(data, ...) UseMethod("plotVariofaces")
plotPredictions <- function(data, ...) UseMethod("plotPredictions")

isCompoundSymmetric <- function(object, ...) UseMethod("isCompoundSymmetric")

#Deprecations

addrm.terms.asreml <- function(...)
{ .Deprecated(new = "changeTerms.asrtests", package = "asremlPlus")
  invisible()
}

addrm.terms.asrtests <- function(...)
{ .Deprecated(new = "changeTerms.asrtests", package = "asremlPlus")
  invisible()
}

alldiffs <- function(...)
{ .Deprecated(new = "as.alldiffs", package = "asremlPlus")
  invisible()
}

asrtests <- function(...)
{ .Deprecated(new = "as.asrtests", package = "asremlPlus")
  invisible()
}

choose.model.asreml <- function(...)
{ .Deprecated(new = "chooseModel.asrtests", package = "asremlPlus")
  invisible()
}

choose.model.asrtests <- function(...)
{ .Deprecated(new = "chooseModel.asrtests", package = "asremlPlus")
  invisible()
}


facRecode <- function(...)
{ .Deprecated(new = "facRecast.alldiffs", package = "asremlPlus")
  invisible()
}

facRecode.alldiffs <- function(...)
{ .Deprecated(new = "facRecast.alldiffs", package = "asremlPlus")
  invisible()
}

info.crit <- function(...)
{ .Deprecated(new = "infoCriteria.asreml", package = "asremlPlus")
  invisible()
}

info.crit.asreml <- function(...)
{ .Deprecated(new = "infoCriteria.asreml", package = "asremlPlus")
  invisible()
}

newrcov.asrtests <- function(...)
{ .Deprecated(new = "changeTerms.asrtests", package = "asremlPlus")
  invisible()
}


plotvariofaces.asreml <- function(...)
{ .Deprecated(new = "plotVariofaces.data.frame", package = "asremlPlus")
  invisible()
}


power.transform <- function(...)
{ .Deprecated(new = "powerTransform", package = "asremlPlus")
  invisible()
}

predictiondiffs.asreml <- function(...)
{ .Deprecated(new = "allDifferences.data.frame", package = "asremlPlus")
  invisible()
}

predictionplot.asreml <- function(...)
{ .Deprecated(new = "plotPredictions.data.frame", package = "asremlPlus")
  invisible()
}

predictparallel.asreml <- function(...)
{ .Deprecated(new = "predictPlus.asreml", package = "asremlPlus")
  invisible()
}

pred.present.asreml <- function(...)
{ .Deprecated(new = "predictPresent.asreml", package = "asremlPlus")
  invisible()
}

recalc.wald.tab.asreml <- function(...)
{ .Deprecated(new = "recalcWaldTab.asrtests", package = "asremlPlus")
  invisible()
}

recalc.wald.tab.asrtests <- function(...)
{ .Deprecated(new = "recalcWaldTab.asrtests", package = "asremlPlus")
  invisible()
}

reml.lrt <- function(...)
{ .Deprecated(new = "REMLRT.asreml", package = "asremlPlus")
  invisible()
}

reml.lrt.asreml <- function(...)
{ .Deprecated(new = "REMLRT.asreml", package = "asremlPlus")
  invisible()
}

reorderClassify  <- function(...)
{ .Deprecated(new = "renewClassify.alldiffs", package = "asremlPlus")
  invisible()
}

reorderClassify.alldiffs  <- function(...)
{ .Deprecated(new = "renewClassify.alldiffs", package = "asremlPlus")
  invisible()
}

rmboundary.asreml <- function( ...)
{ .Deprecated(new = "rmboundary.asrtests", package = "asremlPlus")
  invisible()
}

setvarianceterms.asreml <- function(...)
{ .Deprecated(new = "setvarianceterms.call", package = "asremlPlus")
  invisible()
}


sig.devn.reparam.asreml <- function(...)
{ .Deprecated(new = "reparamSigDevn.asrtests", package = "asremlPlus")
  invisible()
}

sig.devn.reparam.asrtests <- function(...)
{ .Deprecated(new = "reparamSigDevn.asrtests", package = "asremlPlus")
  invisible()
}

testranfix.asreml <- function(...)
{ .Deprecated(new = "testranfix.asrtests", package = "asremlPlus")
  invisible()
}
testrcov.asreml <- function(...)
{ .Deprecated(new = "testresidual.asrtests", package = "asremlPlus")
  invisible()
}
testrcov.asrtests <- function(...)
{ .Deprecated(new = "testresidual.asrtests", package = "asremlPlus")
  invisible()
}
testswapran.asreml <- function(...)
{ .Deprecated(new = "testswapran.asrtests", package = "asremlPlus")
  invisible()
}

# addrm.terms.asrtests,
# choose.model.asrtests,
# info.crit.asreml,
# newfit.asreml,
# newresidual.asrtests,
# plotvariofaces.asreml,
# pred.present.asreml,
# predictiondiffs.asreml,
# predictionplot.asreml,
# predictparallel.asreml,
# recalc.wald.tab.asrtests, 
# reml.lrt.asreml,
# rmboundary.asrtests,
# setvarianceterms.asreml,
# sig.devn.reparam.asrtests,
# simulate.asreml,
# testresidual.asrtests,
# testranfix.asrtests,
# testswapran.asrtests