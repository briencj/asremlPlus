\name{news}
\title{News for Package asremlPlus}
\section{asremlPlus Version 4.3.44 (2022-12-31)}{
  \itemize{
  \item{Add arguments \code{nestorder}, \code{degree} and \code{difforder} to functions addSpatialModelOnIC.asrtests and  makeTPSPlineXZMats.data.frame, as well as argument \code{nestorder} to chooseSpatialModelOnIC.asrtests. Add the option \code{TPLV} to the argument \code{trySpatial} from chooseSpatialModelOnIC.asrtests to allow a linear spatial model with first-order penalty to be fitted.}
  \item{Fix bug in testing for fixed correlations that meant they were not always detected.}
  \item{Fix error in changeModelsonIC.asrtests in the action when the model is \code{Unswapped} - it was being reported as  \code{Changed random}.}  
  \item{Stop crashing of changeTerms.asrtests and ChangeModelOnIC.asrtests when rmboundary.asrtests unable to remove singular terms.}
}}
\section{asremlPlus Version 4.3.43 (2022-11-30)}{
  \itemize{
  \item{Fix bugs in choosing the best spatial model in chooseSpatialModelOnIC.asrtests. Add components \code{best.spatial} and \code{best.spatial.AIC} to the value returned by chooseSpatialModelOnIC.asrtests.}
  \item{Add linear-row-by-linear-column-covariates interaction to the model fitted for TPNCSS option of addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests.}
}}
\section{asremlPlus Version 4.3.42 (2022-11-23)}{
  \itemize{
  \item{Add functions makeTPSPlineZMats.data.frame, addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests.}
}}
\section{asremlPlus Version 4.3.41 (2022-11-16)}{
  \itemize{
  \item{Revise test for fixed correlation so that allow.fixedcorrelation recognises any correlation term, rather than being restricted to ar1.}
}}
\section{asremlPlus Version 4.3.40 (2022-11-16)}{
  \itemize{
  \item{Compress pdf files and resubmit to CRAN.}
}}
\section{asremlPlus Version 4.3.39 (2022-11-13)}{
  \itemize{
  \item{Rebuild for R 4.2.2 and submission to CRAN.}
}}
\section{asremlPlus Version 4.3.38 (2022-10-24)}{
  \itemize{
  \item{Added function pickLSDstatistics.alldiffs.Rd.}
  \item{Revise documentation for the arguments \code{LSDtype}, \code{LSDsupplied}, \code{LSDby} and \code{LSDstatistic }.}
}}
\section{asremlPlus Version 4.3.37 (2022-09-19)}{
  \itemize{
  \item{Fix bugs in function newfit.asreml that involves making changes to models that involve  the \code{at} function. This bug affected the functions changeTerms.asrtests, changeModelOnIC.asrtests, rmboundary.asrtests, testswapran.asrtests and testranfix.asrtests. In removing a term that involve a single lewvel of a factor in an at function, it is now possible to specify either the level or the level index.}
  \item{Modify rmboundary.asrtests to allow the removal of random term for a single level of an \code{at} function when the levels are specified singly in the random model.}
  \item{Add function addto.test.summary.}
}}
\section{asremlPlus Version 4.3.36 (2022-08-17)}{
  \itemize{
  \item{Minor correction for resubmitting to CRAN.}
}}
\section{asremlPlus Version 4.3.35 (2022-08-16)}{
  \itemize{
  \item{Rebuild under R 4.2.1 for new submission.}
}}
\section{asremlPlus Version 4.3-34 (2022-04-20)}{
  \itemize{
  \item{Add functions plotLSDerrors.alldiffs and plotLSDerrors.data.frame.}
  \item{Replace the use of the package \code{reshape} with the use of \code{reshape2}.}
  \item{Replace \code{X1} and \code{X2} with \code{Rows} and \code{Columns} as the default names for the the first two columns of the \code{data.frame} produced by plotPvalues.alldiffs, plotLSDs.alldiffs and plotLSDerrors.alldiffs and provided as input to plotPvalues.data.frame, plotLSDs.data.frame and plotLSDerrors.data.frame.}
  \item{Change plotPvalues.alldiffs and plotLSDs.alldiffs to return a two component list, one component of which is a \code{data.frame} that contains the statistics to be plotted and the other of which is a \code{list} containing the plots that have be produced as \code{ggplot} objects.}
}}
\section{asremlPlus Version 4.3-33 (2022-04-09)}{
  \itemize{
  \item{Add an error if all residual variances are bound.}
}}
\section{asremlPlus Version 4.3-32 (2021-12-21)}{
  \itemize{
  \item{Fix bug in chooseModelOnIC result from a change in as.vector for data.frames and lists in R 4.2.0.}
}}
\section{asremlPlus Version 4.3-31 (2021-12-06)}{
  \itemize{
  \item{Rebuild under R 4.1.2 for new submission.}
}}
\section{asremlPlus Version 4.3-30 (2021-11-22)}{
  \itemize{
  \item{Put in a trap for fixed-effects variances being inestimable when calculating the full likelihood.}
}}
\section{asremlPlus Version 4.3-29 (2021-11-06)}{
  \itemize{
  \item{Modify linTransform.alldiffs so that an intercept-only submodel can be fitted.}
  \item{Revise the manual entry for exploreLSDs.alldiffs and for pairdiffsTransform.alldiffs.}
}}
\section{asremlPlus Version 4.3-28 (2021-10-26)}{
  \itemize{
  \item{Fix bug in newfit.asreml that re-orders the fixed model because of R limitations.}
}}
\section{asremlPlus Version 4.3-27 (2021-10-24)}{
  \itemize{
  \item{Fix problem in linTransform.alldiffs when the product of the numbers of levels in the classiify factors is large and not all combinations of the levels are observed.}
  \item{Relax the requirement for \code{LSDby} variable to be in the \code{classify} to beng that they must be in the \code{predictions} component of the \code{alldiffs.obj}.}
}}
\section{asremlPlus Version 4.3-26 (2021-10-21)}{
  \itemize{
  \item{Fix bug in false.pos and false.neg for \code{LSDtype} set to \code{supplied} when \code{LSDby} is not \code{NULL}.}
  \item{Add quantiles 25 (\code{q25}) and 75 (\code{q75}) to the values that \code{LSDstatistic} can take and the exploreLSDs.alldiffs investigates.}
  \item{Allow for multiple \code{LSDstatistic} values when \code{LSDtype} is \code{factor.combination} and \code{LSDby} is set.}
}}
\section{asremlPlus Version 4.3-25 (2021-10-15)}{
  \itemize{
  \item{Add \code{sig.size}, \code{sig.face} and \code{sig.colour} arguments to plotPvalues.data.frame and plotPvalues.alldiffs.}
}}
\section{asremlPlus Version 4.3-24 (2021-10-14)}{
  \itemize{
  \item{Make sure that when \code{predicted.values} and the associated \code{vcov} and \code{sed} in an \code{\link{alldiffs.object}} have values that are close to zero as a result of the limitations of real number calculations, they are set to zero so as to avoid further calculation problems.}
  \item{Add \code{zero.tolerance} argument to linTransform.alldiffs.}
}}
\section{asremlPlus Version 4.3-23 (2021-10-03)}{
  \itemize{
  \item{Modify the test for fixed correlation from a test for \code{Fixed} (\code{F}) to cover also correlations that are \code{bound} (\code{B}) or \code{singular} (\code{S}).}
  \item{Add tests on incoming fits for fixed correlation.}
}}
\section{asremlPlus Version 4.3-22 (2021-10-01)}{
  \itemize{
  \item{Support is no longer provided in asremlPlus for versions of ASReml-R prior to version 4.0.}
  \item{Add a new argument \code{allow.fixedcorrelation} to changeModelOnIC.alldiffs, newfit.asreml, chooseModel.asrtests, changeTerms.asrtests, testranfix.asrtests, testresidual.asrtests, testswapran.asrtests and reparamSigDevn.asrtests so that a model with a fixed correlation parameter may be rejected.}
}}
\section{asremlPlus Version 4.3-21 (2021-09-23)}{
  \itemize{
  \item{Fix a bug in pairdiffsTransform.alldiffs where the factors indexing the predictions were not generated correctly.}
  \item{In changeModelOnIC.alldiffs, add detection of absent terms to be dropped and do not change the model if one or more are detected.}
}}
\section{asremlPlus Version 4.3-20 (2021-09-03)}{
  \itemize{
  \item{Fix bug in the calculation of false positives and negatives in exploreLSDs.alldiffs for \code{LSDtype} set to \code{overall}.}
  \item{Change default for \code{plotHistogram} in exploreLSDs.alldiffs to \code{FALSE}.}
}}
\section{asremlPlus Version 4.3-19 (2021-09-01)}{
  \itemize{
  \item{Rebuild under R 4.1.1 for new submission.}
}}
\section{asremlPlus Version 4.3-18 (2021-08-31)}{
  \itemize{
  \item{Change number of comparisons symbol in an \code{LSD.frame} from \code{n} to \code{c}.}
}}
\section{asremlPlus Version 4.3-17 (2021-08-30)}{
  \itemize{
  \item{Add False positives and negatives counts to LSD statistics.}
}}
\section{asremlPlus Version 4.3-16 (2021-08-24)}{
  \itemize{
  \item{Fix bug in redoErrorIntervals when only one factor in \code{LSDby}.}
  \item{Add function print.LSDdata.}
}}
\section{asremlPlus Version 4.3-15 (2021-08-19)}{
  \itemize{
  \item{Fixed bug in \code{halfLeastSignificant.limit} when \code{LSDtype} is \code{factor.combination} and the \code{assignedLSD} is zero to now ensure that the limit is \code{NA}.}
  \item{Reinstated passing of parameters from linTransform.alldiffs to redoErrorIntervals.alldiffs.}
}}
\section{asremlPlus Version 4.3-14 (2021-08-18)}{
  \itemize{
  \item{Fix newly introduced  bug for \code{LSDtype = "halfLeastSignificant"} when there is not a numeric covariate involved in the \code{LSDby} argument.}
  \item{Stop extra frequency plot produced by exploreLSDs.}
  \item{Minor internal modification.}
}}
\section{asremlPlus Version 4.3-13 (2021-08-16)}{
  \itemize{
  \item{Fix bug for \code{LSDtype = "halfLeastSignificant"} when there is a numeric covariate involved in the \code{LSDby} argument.}
  \item{Minor internal modification.}
}}
\section{asremlPlus Version 4.3-12 (2021-08-11)}{
  \itemize{
  \item{Fix bug in pairdiffsTransform when only two factors in \code{classify}.}
  \item{Add function subset.list.}
  \item{Add functions exploreLSDs, plotLSDs.alldiffs, plotLSDs.matrix and plotLSDs.data.frame.}
  \item{Incorporate arguments \code{retain.zerLSDs} and \code{zero.tolerance} into AllDifferences.data.frame, redoErrotIntervals.alldiffs amd exploreLSDs.alldiffs to allow control of the retention of the zero LSD values when calculating LSD summaries.}
}}
\section{asremlPlus Version 4.3-11 (2021-07-28)}{
  \itemize{
  \item{Deprecate the function asrtests, which has been replaced by the function as.asrtests.}
  \item{There has been a substantial revision of the calculation and saving of the LSD component of an \code{alldiffs} object. 
    \itemize{
  \item{The argument \code{meanLSD.type} has been changed to \code{LSDtype}. All functions that are S3 methods for \code{alldiffs} objects will check for the \code{meanLSD.type} attribute and, if found, change it to \code{LSDtype}.} 
  \item{The option \code{supplied} has been added to the argument \code{LSDtype} so that the user can supply LSD values to replace those in the \code{assignedLSD} column of the \code{LSD} component.}
  \item{The two new arguments \code{LSDstatistic} and \code{LSDaccuracy} have been added. The argument \code{LSDstatistic} allows a choice of \code{minimum}, \code{mean}, \code{median} or \code{maximum}. The values of the nominated LSD statistics are stored in a column of the LSD component named \code{assignedLSD}. The argument \code{LSDaccuracy} allows the specification of a statistic summarizing the differences between the the LSDs and the assigned LSD value for each row of the LSD component. The values of the chosen measure are stored in a column of the \code{LSD} component named \code{accuracyLSD}. These two arguments have been added to the functions allDifferences.data.frame, as.alldiffs.object, as.alldiffs, linTransform.alldiffs, pairdiffsTransform, predictPlus.asreml, predictPresent.asreml, ratioTransform.alldiffs, recalcLSD.alldiffs and redoErrorIntervals.alldiffs.}
  \item{The \code{LSD} component is always calculated and stored as a component when \code{pairwise} is \code{TRUE}.}}}
  \item{A new argument \code{accuracy.threshold} has been introduced to control the addition of a column names \code{LSDwarning} to the \code{predictions} of an \code{alldiffs} object. THis column indicates whether or not the an \code{halfLeastSignificant} \code{error.interval} is reliable as determined using the supplied value of  \code{accuracy.threshold}.}
  \item{The attributes of an alldiffs object  and of the \code{predictions} component have been modified. The arguments \code{LSDtype}, \code{LSDstatistic} and \code{LSDby}, if set, are attributes of an \code{alldiffs} object. They are also attributes of a \code{prediction.frame} and the \code{backtransforms} component of an \code{alldiffs} object, as is \code{LSDvalues}, when \code{halfLeastSignificant} \code{error.intervals} are specified. The \code{LSDvalues} attribute has the values used in calculating half-LSD intervals. The \code{meanLSD} attribute has been removed.}
  \item{A separate entry for \code{LSD.frame} has been introduced into the manual to dedscribe the \code{data.frame} that is stored in the \code{LSD} component of an \code{alldiffs} object.}
}}
\section{asremlPlus Version 4.3-10 (2021-07-5)}{
  \itemize{
  \item{Fix bug in facRecast that meant that the argument \code{levels.order} was being ignored.}
  \item{Fix a bug in facRecast where only the predictions component was appropriately re-ordered.}
  \item{Add new functions ratioTransform.alldiffs and pairdiffsTransform.}
  \item{Add new method to sort for a \code{predictions.frame}: function sort.predictions.frame.}
  \item{Change the argument \code{sortWithinVals} to \code{sortParallelToCombo} in sort.alldiffs.}
  \item{Add the new argument \code{sortNestinFactor} to sort.alldiffs.}
  \item{Add new function getTestEntry.asrtests.}
  \item{Fix bug in infoCriteria.asreml that resulted in values of \code{Inf} for the \code{full} likelihood.}
}}
\section{asremlPlus Version 4.2-36 (2021-05-16)}{
  \itemize{
  \item{Deprecate facRecode.alldiffs because it only make sense to redo the current levels or labels.}
  \item{Tidy up function facRecast.alldiffs.}
}}
\section{asremlPlus Version 4.2-35 (2021-05-16)}{
  \itemize{
  \item{Add new function facRecast.alldiffs that supercedes facRecode.alldiffs.}
  \item{Revise sort.alldiffs to use facRecast.alldiffs to reorder the alldiffs components. This fixed a bug when not all levels of the \code{saortFactor} are present in the subset on which the sort is based.}
  \item{Add argument \code{alpha} to plotPvalues.data.frame and plotPvalues.alldiffs.}
}}
\section{asremlPlus Version 4.2-34 (2021-05-12)}{
  \itemize{
    \item{Add argument \code{axis.text.size} to plotPvalues.data.frame and plotPvalues.alldiffs.}
}}
\section{asremlPlus Version 4.2-33 (2021-05-06)}{
  \itemize{
    \item{Modify testresidual to not change the variance model when there is no difference in the variance DF between the models being compared.}
    \item{In as.asrtests, add variance DF to \code{test.summary} when forming an \code{asrtests} object.}
}}
\section{asremlPlus Version 4.2-32 (2021-03-19)}{
  \itemize{
    \item{Rebuild vignerttes.}
    \item{Modify more links to asreml.}
}}
\section{asremlPlus Version 4.2-31 (2021-03-19)}{
  \itemize{
    \item{Modify link to asreml.}
}}
\section{asremlPlus Version 4.2-30 (2021-03-18)}{
  \itemize{
    \item{Revise doi specification in a link to a reference in the manual.}
    \item{Increment version number for resubmission.}
}}
\section{asremlPlus Version 4.2-29 (2021-03-18)}{
  \itemize{
    \item{Increment version number for resubmission.}
}}
\section{asremlPlus Version 4.2-28 (2021-03-17)}{
  \itemize{
    \item{Recompile for CRAN submission.}
}}
\section{asremlPlus Version 4.2-27 (2021-03-13)}{
  \itemize{
    \item{Fix bugs in redoErrorIntervals.alldiff concerning whether or not meanLSD.type and LSDby settings result in the LSD component being recalculated.}
    \item{Add rowname \code{overall} to the LSD component when \code{meanLSD.type} is set to \code{overall}}
}}
\section{asremlPlus Version 4.2-26 (2020-11-11)}{
  \itemize{
    \item{Modify URLs for CRAN submission.}
}}
\section{asremlPlus Version 4.2-25 (2020-11-11)}{
  \itemize{
    \item{Try to neutralize problem that occurs with \code{colourise} in print.predictions.frame.}
    \item{Rebuild for R 4.0.3 and recompile vignettes for CRAN submission.}
}}
\section{asremlPlus Version 4.2-24 (2020-09-25)}{
  \itemize{
    \item{Include a \code{key} option in the \code{which} that produces \code{vparametersummary} rather than the full \code{asremlsummary}.}
    \item{Remove attributes from formuale when executing printFormulae.}
    \item{Add the ChickpeaEnd data set.}
}}
\section{asremlPlus Version 4.2-23 (2020-09-07)}{
  \itemize{
    \item{Modify code so that error intervals are set to NA when the LSDs are zero.}
}}
\section{asremlPlus Version 4.2-22 (2020-09-04)}{
  \itemize{
    \item{Modify code for calculating LSDs so that all zero LSDs are allowed; this is needed when a linear transformation has resulted in predictions that are equal.}
}}
\section{asremlPlus Version 4.2-21 (2020-07-11)}{
  \itemize{
    \item{Rebuild the vignettes prior to CRAN submission.}
    \item{Add \code{envir} argument to getFormulae.asreml and printFormulae.asreml to make them more robust to where the objects in the \code{asreml.obj} have been defined.}
}}
\section{asremlPlus Version 4.2-20 (2020-06-29)}{
  \itemize{
    \item{Fix bug in getFormulae.asreml and printFormulae.asreml that occurs when formula(e) have been supplied to asreml in an object.}
}}
\section{asremlPlus Version 4.2-19 (2020-04-02)}{
  \itemize{
    \item{Change the Maintainer email address.}
    \item{Add ellipsis to as.asrtests call within an rmboundary.asrtests call so that arguments in the calling routine are passed through.}
    \item{Add a trap in recalcWaldTab.asrtests for a \code{NULL} \code{wald.tab} before dealing with \code{denDF} so that recalcWaldTab.asrtests throws a warning rather than an error.}
    \item{Modify recalcWaldTab.asrtests to deal with a Wald table that has Chi-square instead of F statistics (this happens when \code{denDF = "none"}.)}
}}
\section{asremlPlus Version 4.2-18 (2020-03-16)}{
  \itemize{
    \item{Increment version number for resubmission to CRAN.}
}}
\section{asremlPlus Version 4.2-17 (2020-03-16)}{
  \itemize{
    \item{Fix bug in changeTerms when both \code{addFixed} and \code{dropFixed} are \code{NULL}.}
}}
\section{asremlPlus Version 4.2-16 (2020-02-25)}{
  \itemize{
    \item{Revisions to the manual.}
    \item{Fix bug in forming test.summary.}
}}
\section{asremlPlus Version 4.2-15 (2020-02-17)}{
  \itemize{
    \item{Added the \code{IClikelihood} argument to chooseModel.asrtests, reparamSigDevn.asrtests, rmboundary.asrtests, testranfix.asrtests, testresidual.asrtests, and testswapran.asrtests.}
    \item{Fix bug in printFormulae when the formula is long.}
    \item{Various bug fixes associated with asreml-R version 3.}
}}
\section{asremlPlus Version 4.2-14 (2020-02-14)}{
  \itemize{
    \item{Fix bug in processing a formula that includes an \code{at} function.}
}}
\section{asremlPlus Version 4.2-13 (2020-02-11)}{
  \itemize{
    \item{Removed the \code{material.diff} argument and the \code{both} option from the \code{which.IC} argument of changeModelOnIC.}
}}
\section{asremlPlus Version 4.2-12 (2020-02-09)}{
  \itemize{
    \item{Change \code{likelihood} in infoCriteria to \code{IClikelihood} to make it consistent with other functions.}
    \item{Replace \code{tol.diff} with \code{material.diff} and set the default to 0.5 in changeModelOnIC.}
    \item{Revise all vignettes.}
}}
\section{asremlPlus Version 4.2-11 (2020-02-08)}{
  \itemize{
    \item{Add a new vignette for the use of information criteria with the wheat experiment.}
    \item{Fixed bugs in changeModelOnIC associated with \code{addFixed} and \code{dropFixed} and with dealing with unconverged models.}
}}
\section{asremlPlus Version 4.2-10 (2020-02-07)}{
  \itemize{
    \item{This version represents a major revision in that the facilities for employing information criteria (AIC and BIC) have been considerabley expanded. The infoCriteria method has been modified to add the possibility of calculating not only the the information criteria based on the REML likelihood, but those based on the full likelihood instead. There are now two infoCriteria methods, one for \code{asreml} objects and the other for \code{list}s of \code{asreml} objects. An alternative is to use changeTerms.asrtests wiith \code{IClikelihood} set to \code{REML} or \code{full} and the information criteria will be included in the The \code{test.summary} of the \code{asrtest} object. A new function changeModelOnIC has been included to allow model selection to be based on the information criteria, rather than on hypothesis tests. However, the two methods of selecting models can be interweaved.}
    \item{Add  argument \code{omit.columns} to print.test.summary that allows the specification of columns to be omitted from the printout.}
    \item{Various modifications to the internal code for the handling of the \code{test.summary} component of an \code{asrtest.object} to allow more flexibility in the make up of \code{test.summary}.}
    \item{Add function iterate.asrtests to execute additional iterations of the fit stored in an  \code{asreml.obj} within an \code{asrtests} object; the model itself cannot be changed.}
    \item{Add functions getFormulae.asreml and printFormulae.asreml to extract and print the \code{formula}e stored in the \code{call} of an \code{asreml} object.}
    \item{Add \code{trace} argument to newfit.asreml and activated it for asreml-R version 4.} 
    \item{Made R >= 3.5.0 a requirement for asremlPlus.}
}}
\section{asremlPlus Version 4.1-37 (2020-01-27)}{
  \itemize{
    \item{Add chooseModel.data.frame function.}
    \item{Modify linTransform.alldiffs so that submodels involving covariates are possible.}
    \item{Add print.test.summary function}
    \item{Split Ladybird vignette into two, one based on asreml and the other on lm, lmerTest and emmeans.}
    \item{Modified newfit to circumvent problems in the test functions in dealing with terms that include the \code{at} function and modified testranfix to improve its detection of multiple terms.}
    }
}
\section{asremlPlus Version 4.1-36 (2020-01-12)}{
  \itemize{
    \item{Correct date in DESCRIPTION for CRAN rsubmission.}
    }
}
\section{asremlPlus Version 4.1-35 (2020-01-11)}{
  \itemize{
    \item{Rebuild for CRAN submission.}
    }
}
\section{asremlPlus Version 4.1-34 (2019-12-14)}{
  \itemize{
    \item{Change default for \code{subset} argument in subset.alldiffs so that the complete set of predictions is returned. Also, change check on the resulting classify so that it checks whether the variables in the classify uniquely index the predictions.}
    \item{Include in renewClassify.alldiffs a check that the \code{newclassify} uniquely indexes the predictions.}
    }
}
\section{asremlPlus Version 4.1-33 (2019-11-08)}{
  \itemize{
    \item{Add \code{printPlot = FALSE} to both plotPvalues.data.frame and plotPvalues.alldiffs. Change return for plotPvalues.data.frame from nothing to the plot produced using ggplot.}
    }
}
\section{asremlPlus Version 4.1-32 (2019-09-30)}{
  \itemize{
    \item{Fix bug in check for predict.asreml returning sed and/or vcov when requested.}
    }
}
\section{asremlPlus Version 4.1-31 (2019-09-19)}{
  \itemize{
    \item{Add new functions facRecode.alldiffs and facRename.alldiffs.}
    }
}
\section{asremlPlus Version 4.1-30 (2019-09-05)}{
  \itemize{
    \item{Implement a check that predict.asreml has returned sed and/or vcov when requested.}
    }
}
\section{asremlPlus Version 4.1-29 (2019-07-08)}{
  \itemize{
    \item{Allow for calculation of notional LSDs when LSDs for a single observation are specified.}
    \item{Removed restriction that overall LSD is always used when the range of the SEDs divided by the average of the SEDs is less than the \code{avsed.tolerance}.}
    \item{A bug in detecting valid wald.tab objects was fixed.}
    }
}
\section{asremlPlus Version 4.1-28 (2019-07-04)}{
  \itemize{
    \item{Modify testranfix.asrtests, as.asrtests, recalcWaldTab.asrtests and chooseModel.asrtests to allow for p-values based on a conditional Wald F-statistic if one is present in the \code{wald.tab}.}
    \item{Allow for 6-column \code{wald.tab} such as is obtained for \code{ssType} set to \code{"conditional"}.}
    }
}
\section{asremlPlus Version 4.1-27 (2019-06-01)}{
  \itemize{
    \item{Revise Wheat vignette to use qqplotr for Normal Probability Plots.}
    \item{Revise test for tilde (~) in formula supplied to changeTerms so that only a leading tilde causes an error.}
    }
}
\section{asremlPlus Version 4.1-26 (2019-05-12)}{
  \itemize{
    \item{Rebuild for R 3.6.0.}
    }
}
\section{asremlPlus Version 4.1-25 (2019-04-18)}{
  \itemize{
    \item{Fix a bug in changeTerms.asrtests that reset the \code{residual} model to the default when changing  fixed or random terms.}
    \item{Check for \code{NULL} \code{wald.tab} when setting its \code{class} prior to storing it in an \code{asrtests.object}.}
    \item{Add function addBacktransforms.alldiffs.}
    }
}
\section{asremlPlus Version 4.1-24 (2019-04-11)}{
  \itemize{
    \item{Allow numerics, in addition to factors, in \code{LSDby} argument of predictPlus.}
    \item{Rename reorderClassify to renewClassify.}
    \item{Fix a bug in retention of attributes of the \code{alldiffs.object} in  renewClassify.}
    \item{Add \code{transform.power}, \code{offset} and \code{scale} as attributes of the \code{backtransforms} component of an \code{alldiffs.object}.}
    \item{Rework the storage of the \code{backtransforms} component of an \code{alldiffs.object} to ensure that it is retained when manipulating \code{alldiffs.object}.}
    }
}
\section{asremlPlus Version 4.1-23 (2019-04-03)}{
  \itemize{
    \item{Build in R version 3.5.3.}
    \item{Update the Author field of the DESCRIPTION to include an ORCID.}
    \item{Allow backtransformed halfLSD intervals in the \code{backtransforms} component of an \code{alldiffs.object}.}
    }
}
\section{asremlPlus Version 4.1-22 (2019-03-04)}{
  \itemize{
    \item{Build in R version 3.5.2.}
    \item{Update the license to remove blank line and resubmit to CRAN.}
    }
}
\section{asremlPlus Version 4.1-21 (2019-03-01)}{
  \itemize{
    \item{Build in R version 3.5.2.}
    \item{Update the license to include 2019.}
    }
}
\section{asremlPlus Version 4.1-20 (2019-02-23)}{
  \itemize{
    \item{Add print.predictions.frame and print.wald.tab as S3 methods to give more control over the printing of objects of class \code{wald.tab} and \code{predictions.frame}.}
    \item{Made sure that class \code{wald.tab} is associated with a \code{wald.tab} and classes asreml.predict and predictions.frame are associated with a \code{predictions.frame} when these are produced by asreml.}
    }
}
\section{asremlPlus Version 4.1-19 (2019-02-22)}{
  \itemize{
    \item{Add \code{colourise} argument to print.alldiffs so that the default is FALSE and it overrides the \code{asreml.options} setting.}
    }
}
\section{asremlPlus Version 4.1-18 (2019-02-21)}{
  \itemize{
    \item{Use R.rsp to install vignettes so that hyperlinks are active.}
    \item{Add \code{colourise} argument to print.asrtests so that the default is FALSE and it overrides the \code{asreml.options} setting.}
    }
}
\section{asremlPlus Version 4.1-17 (2019-02-16)}{
  \itemize{
    \item{Fix bug in infoCriteria.asreml when DF supplied. For both infoCriteria.asreml and REMLRT.asreml, calculate the bound terms and if appropriate give arning messages.}
    }
}
\section{asremlPlus Version 4.1-16 (2019-02-07)}{
  \itemize{
    \item{More vignettes changes, including adding manual to vignettes.}
    \item{Fix accidental reordering of \code{alldiffs} object when \code{meanLSD.type} is \code{factor combinatons}.}
    }
}
\section{asremlPlus Version 4.1-15 (2019-02-05)}{
  \itemize{
    \item{Make vignettes visible.}
    \item{Modify estimateV in response to correction of corb naming in latest build.}
    }
}
\section{asremlPlus Version 4.1-14 (2019-02-02)}{
  \itemize{
    \item{Revise Wheat vignette and entry in asremlPlus-package documentation.}
    \item{Add Ladybird predictions vignette.}
    \item{Modify printing of \code{wald.tab} by print.asrtests to remove colouring because it causes problems with knitr.}
    }
}
\section{asremlPlus Version 4.1-13 (2019-01-30)}{
  \itemize{
    \item{Add a check to plotPredictions.data.frame to warn when \code{linear.transformation} is supplied.}
    \item{Also remove zero LSDs when calculating overall LSD and sed.range.}
    \item{Add a check to linTransform.alldiffs that factors in LSDby are in the submodel.}
    \item{Add the Wheat example as a vignette.}
    }
}
\section{asremlPlus Version 4.1-12 (2019-01-26)}{
  \itemize{
    \item{Remove LSDs for which sed is effectively zero from the calculation of LSD statistics, provided not all are zero.}
    }
}
\section{asremlPlus Version 4.1-11 (2019-01-21)}{
  \itemize{
    \item{Fix bug in test of the range of SEDs against \code{avsed.tolerance} when caclulating LSDs.}
    \item{Add attributes \code{meanLSD.type} and \code{LSDby} to \code{alldiffs} objects so that redoErrorIntervals.alldiffs can determine these attributes when its \code{meanLSD.type} argument is not set.}
    \item{Add function getTestPvalue.asrtests.}
    \item{Minor documentation changes.}
  }
}
\section{asremlPlus Version 4.1-10 (2019-01-08)}{
  \itemize{
    \item{Add entries for \code{predictions.frame}, \code{alldiffs.object} and \code{asrtests.object} and made changes to the documentation to provide links to these pages.}
    \item{Defined classes alldiffs, asrtests and predictions.frame.}
    \item{Added the functions is.alldiffs, validAlldiffs, is.asrtests, validAsrtests, is.predictions.frame and validPredictionsFrame for validating the objects.}
    \item{Added tests for the validity of objects to methods functions.}
    \item{Provided as.asrtests function, while retaining asrtests function as a pseudonym for as.asrtests. Added an as.prediction.frame function.}
    \item{Add examples for \code{alldiffs} and \code{data.frame} methods functions that use lme4 to show how to use asremlPlus with packages othere than asreml.}
  }
}
\section{asremlPlus Version 4.1-08 (2018-12-10)}{
  \itemize{
    \item{Correct predictPlus, plotPredictions and allDifferences.dat.frame to deal with a \code{classify} equal to \code{"(Intercept)"}.}
    \item{Implement fix of \code{R.param} and \code{G.param} in supplied object for \code{corb} in estimateV function.}
    \item{Make documentation changes for asreml-R4.1.}
  }
}
\section{asremlPlus Version 4.1-07 (2018-09-05)}{
  \itemize{
    \item{Add Ladybird.dat data set.}
  }
}
\section{asremlPlus Version 4.1-06 (2018-09-03)}{
  \itemize{
    \item{Fix minor bugs in redoErrorIntervals.alldiffs.}
  }
}
\section{asremlPlus Version 4.1-05 (2018-08-29)}{
  \itemize{
    \item{Allow no \code{Vmatrix} in linTransform.alldiffs when \code{error.intervals} is \code{none} and allow transformation of \code{predicted.value} only when \code{vcov component} of \code{alldiffs} object is not available.}
    \item{Allow no \code{sed} in allDifferences.data.frame when not \code{pairwise} is \code{TRUE}.}
  }
}
\section{asremlPlus Version 4.1-04 (2018-08-16)}{
  \itemize{
    \item{Modify legend size and position in plotPvalues.alldiffs.}
  }
}
\section{asremlPlus Version 4.1-03 (2018-08-12)}{
  \itemize{
    \item{Add new function reorderClassify.alldiffs.}
    \item{Add \code{rmClassifyVars} to subset.alldiffs.}
  }
}
\section{asremlPlus Version 4.1-02 (2018-06-13)}{
  \itemize{
    \item{Rebuild under latest version of R 3.5.0 with requirement for R >= 3.1.0.}
    \item{Fix bug in bootREMLRT.}
  }
}
\section{asremlPlus Version 4.1-01 (2018-06-03)}{
  \itemize{
    \item{Changes to make asremlPlus with the release version of asreml version 4.1.}
    \item{Add new functions getASRemlVersionLoaded and loadASRemlVersion to facilitate handling of different versions of asreml.}
    \item{Fixed bug in predictPlus that meant the denominator degrees of freedom could not be obtained from \code{wald.tab} when \code{error.intervals} was set to \code{"none"} or \code{"StandardError"}.}
    \item{Fixed bug in \code{LSDby} option that meant the correct subsets were not being identified when the levels of one factor in the predictions are subset of another factor.}
  }
}
\section{asremlPlus Version 4.0-28 (2018-04-18)}{
  \itemize{
    \item{Add new function facCombine.alldiffs.}
    \item{Fix minor bug in plotPvalues.alldiffs associated with \code{factors.per.grid} when there ais more than one factor indexing the different grids.}
  }
}
\section{asremlPlus Version 4.0-27 (2018-03-23)}{
  \itemize{
    \item{Add new function facCombine.alldiffs.}
    \item{Add \code{factors.per.grid} plotPvalues.alldiffs to allow automatic calculation of grid spacing.}
  }
}
\section{asremlPlus Version 4.0-26 (2018-03-16)}{
  \itemize{
    \item{Add argument \code{linear.transformation} to predictPresent.}
  }
}
\section{asremlPlus Version 4.0-25 (2018-03-15)}{
  \itemize{
    \item{Minor bug fix in print.alldiffs.}
  }
}
\section{asremlPlus Version 4.0-24 (2018-03-15)}{
  \itemize{
    \item{Add argument \code{vcov} to allDifferences.data.frame and as.alldiffs and modify print.alldiffs, recalcLSD.alldiffs, sort.alldiffs and subset.alldiffs to process the \code{vcov} component of an \code{alldiffs} object.}
    \item{Add argument \code{linear.transformation} to predictPlus to enable a linear transformation of the predictiions as specified by a \code{matrix} or a model \code{formula}.}
    \item{Add function linTransform.alldiffs that effects a linear transformation of the predictions stored in an \code{alldiffs} object.}
    \item{Fixed bug in backtransforming predictions in predictPlus when offset is nonzero and scale is not one.}
  }
}
\section{asremlPlus Version 4.0-23 (2018-02-19)}{
  \itemize{
    \item{Fix bug in gridline in plotPvalues when\code{triangles} is set to \code{upper} or \code{lower}.}
  }
}
\section{asremlPlus Version 4.0-22 (2018-02-01)}{
  \itemize{
    \item{Fix bug in reordering components of an alldiffs object.}
  }
}
\section{asremlPlus Version 4.0-21 (2018-01-25)}{
  \itemize{
    \item{Change license to MIT.}
  }
}
\section{asremlPlus Version 4.0-20 (2017-12-28)}{
  \itemize{
    \item{Modifications to checking of options for function arguments.}
    \item{Add tip for installing asremlPlus from github.}
    \item{Modifications to internals.}
  }
}
\section{asremlPlus Version 4.0-19 (2017-12-18)}{
  \itemize{
    \item{Fix bug in redoErrorIntervals.alldiffs to make sure \code{status} or \code{est.status} is the last column in the \code{predictions} component of the \code{alldiffs} object.}
    \item{Modifications to internals.}
  }
}
\section{asremlPlus Version 4.0-18 (2017-12-11)}{
  \itemize{
    \item{Add \code{colours} argument to plotPvalues.data.frame and to plotPvalues.alldiffs and add \code{triangles} argument to plotPvalues.alldiffs. In plotPvalues.data.frame, move \code{legend.position} to bottom of plot and increase \code{legend.width} to 2 cm. Also, remove \code{panel.grid} from plot and enforce an \code{aspect.ratio} of one.}
    \item{Fix bug in subset.asreml that did not reduce the levels of the in the prediction component of an \code{alldiffs} object to those remaining after subsetting.}
  }
}
\section{asremlPlus Version 4.0-17 (2017-12-04)}{
  \itemize{
    \item{Rebuild the package under the latest version of R.}
  }
}
\section{asremlPlus Version 4.0-16 (2017-11-20)}{
  \itemize{
    \item{Add code to allDifferences.data.frame that ensures that the order of the predictions, and associated compoenents are in standard order for the classify on which the predictions are based.}
  }
}
\section{asremlPlus Version 4.0-15 (2017-11-17)}{
  \itemize{
    \item{Add subset.alldiffs function.}
    \item{Fix bug in REMLRT.asreml related to bound.exclusions with asreml3.}
    \item{Build under the latest version of R.}
  }
}
\section{asremlPlus Version 4.0-14 (2017-10-26)}{
  \itemize{
    \item{Have modified testresidual.asrtests, testranfix.asrtests and testswapran.asrtests to fit the new random/residual model when the old random/residual model has not converged and \code{allow.unconverged} is set to \code{FALSE}. In such cases the value of p will be \code{NA}. For these functions, when airthe the new or the old model has failed to converge, an indication of failure to converge will be added to the \code{action} reported in \code{test.summary}.}
  }
}
\section{asremlPlus Version 4.0-13 (2017-10-22)}{
  \itemize{
    \item{Modified \code{meanLSD.type} argument in allDifferences.data.frame, predictPlus.asreml and predictPresent.asreml to include a \code{factor.combinations} option that obtains the mean LSD for each combinaton of the \code{levels} of one or more \code{factors}.}
    \item{Added functions recalcLSD.alldiffs and redoErrorIntervals.alldiffs to enable recalculations on \code{alldiffs} objects.}
    \item(Added call to asreml.options with \code{keep.order = TRUE} in newfit.asreml when asreml4 is loaded, because asreml4 requires it.)
  }
}
\section{asremlPlus Version 4.0-12 (2017-10-09)}{
  \itemize{
    \item{Modified recalcWaldTab.asrtests to pass empty \code{wald.tab} when wald.asreml has failed.}
    \item{Made minor changes to testresidual.asrtests to make clear the still use \code{rcov} in asreml3.}
 }
}
\section{asremlPlus Version 4.0-11 (2017-07-13)}{
  \itemize{
    \item{Fixed bug in REMLRT.asreml and rmboundary.asrtests.}
 }
}
\section{asremlPlus Version 4.0-10 (2017-07-06)}{
  \itemize{
    \item{Fixed bugs in allDifferences.data.frame and plotPvalues.alldiffs.}
    \item{Compiled with the latest version of R}
    \item{Some internal recoding to replace calls to summary.asreml.}
 }
}
\section{asremlPlus Version 4.0-09 (2017-05-14)}{
  \itemize{
    \item{Added \code{gridspacing} option to plotPvalues.alldiffs and plotPvalues.data.frame so that a grid could be formed over groups of rows and columns``.}
    \item{Added function sort.alldiffs to sort the components of an \code{alldiffs} object according to the predicted values. Add arguments and a call to it to predictPlus.asreml, predictPresent.asreml and plotPvalues.alldiffs.}
    \item{Added attributes \code{sortFactor} and \code{sortOrder} to an \code{alldiffs} object so that this is accessible to users for use in other calls.}
    \item{Added the oats split-plot design from Yates (1937).}
 }
}
\section{asremlPlus Version 4.0-08 (2017-05-4)}{
  \itemize{
    \item{Minor change to as.alldiffs to make sure that using \code{standard.error} and \code{est.status} rather than \code{std.error} and \code{status} in the predictions data.frame.}
 }
}
\section{asremlPlus Version 4.0-07 (2017-4-12)}{
  \itemize{
    \item{Added \code{which.matrix} argument to estimateV so that could estimate the varaince, random or residual matrix.}

    \item{Incorporate newresidual.asrtests into changeTerms.asrtests and added check for  boundary and singular terms.}
 }
}
\section{asremlPlus Version 4.0-06 (2017-3-25)}{
  \itemize{
    \item{Added \code{checkboundaryonly} argument to rmboundary.asrtest and all functions that call it so that could control removal of boundary and singular terms.}

    \item{Incorporate newresidual.asrtests into changeTerms.asrtests and added check for  boundary and singular terms.}
 }
}
\section{asremlPlus Version 4.0-05 (2017-3-19)}{
  \itemize{
    \item{Change full.asreml.obj and reduced.asreml.obj to h1.asreml.obj and h0.asreml.obj in REMLRT.asreml and reordered them.}
 }
}
\section{asremlPlus Version 4.0-04 (2017-3-16)}{
  \itemize{
    \item{Add function bootREMLRT.asreml and replace addrm.terms.asrtests with changeTerms.asrtests.}
 }
}
\section{asremlPlus Version 4.0-03 (2017-2-26)}{
  \itemize{
    \item{Make the default in variofaces.asreml for \code{V} be \code{NULL}, in which case V is obtained using estimateV.asreml.}
 }
}
\section{asremlPlus Version 4.0-02 (2017-2-19)}{
  \itemize{
    \item{Fix bug in plotPvalues arising from asreml not always returning \code{factors} in the predicted values.}
 }
}
\section{asremlPlus Version 4.0-01 (2017-2-12)}{
  \itemize{
    \item{Modify varofaces.asreml and simulate.asreml to perform the simulations and analysis of data sets using parallel processng.}
    \item{Add recalcWaldTab to asrtests so that wald.tab options for denominator degrees of freedom are taken into account whenever wald.tab is recalculated.}
    \item{Add arguments \code{DF} and \code{bound.exclusions} to REMLRT.asreml and infoCriteria.asreml to allow more control over how the DF are determined. The number of bound parameters was added to the data.frame produced by these methods. Also, the arguments \code{bound.exclusions} and \code{REMLDF} were added to the methods testranfix.asrtests, testresidual.asrtests and testswapran.asrtests. }
    \item{Added \code{meanLSD.type} to allDifferences.data.frame, predictPlus.asreml and predictPresent.asreml to provide for an overall or per-prediction meanLSD. Allow \code{avsed.tolerance} in predictPlus to be set to \code{NA} which forces the calculation of LSDs for the overall mean LSD and set calculations of per-prediction LSDs whatever the value of code{avsed.tolerance}. Changed to calculate the mean SED to be used for an LSD as the square root of the mean(s) of the variances of pairwise differences.}
    \item{Added functions plotPvalues.data.frame and plotPvalues.alldiffs to plot heat maps of p-values for pairwise differences between predictions.}
    \item{Added function estimateV.asreml to form the estimated variance matrix for the observations. However, this function is not available in ASReml-R version 3.}
 }
}
\section{Changes in asremlPlus Version 4.0-0 (2017-1-1)}{
  The asremlPlus version 4.x package represents a major revamp of asremlPlus to make it compatible with ASReml-R version 4 and, at the same time, to make the functions that operate on \code{asreml}, \code{asrtests} and \code{alldiffs} objects into S3 methods for these S3 objects. This means that
  \itemize{
    \item{The function names have been modified so that the only period (`.') in the name is that before the object type for which the function is a method;}
    \item{The object part of the function name is no longer needed when calling the function;}
    \item{The arguments had to be reordered for some functions so that the argument for the object for which it is a method is the first argument; this may cause problems if matching of arguments by position has been used in calling functions;}
    \item{The names of some arguments have been changed because of changes in ASReml-R version 4; for example, \code{rcov} has been changed to \code{residual} throughout.} 
}    

The full set of changes to function names is detailed in the asremlPlus manual and can be accessed using \code{?`asremlPlus-deprecated`}.

In spite of th changes, the new version still works with ASReml-R version 3, although with changes. A particularly problematic change for version 3 is that testrcov.asrtests becomes testresidual.asrtests even though \code{rcov} is still used in calls to \code{asreml}.
}
\section{asremlPlus Version 3.0-0 (2016-11-19)}{
  \itemize{
    \item{Version 3 of asremlPlus was skipped in order to line up the version numbers of asreml and asremlPlus.}
 }
}
\section{Changes in asremlPlus Version 2.0-13 (2016-09-26)}{
  \itemize{
    \item{Add argument \code{allow.unconverged} in addrm.terms.asrtests, choose.model.asrtests, newfit.asreml, newrcov.asrtests, sig.devn.reparam.asrtests, testranfix.asreml, testrcov.asrtests, testswapran.asrtests. It allows control over whether an unconverged analysis is to be abandoned.}
 }
}
\section{Changes in asremlPlus Version 2.0-12 (2016-09-16)}{
  \itemize{
    \item{Rebuild under R 3.4.0.}
 }
}
\section{Changes in asremlPlus Version 2.0-11 (2016-08-13)}{
  \itemize{
    \item{Added check that \code{wald.tab} is a four-column data.frame to predictionplot.asreml and pred.present.asreml.}
 }
}
\section{Changes in asremlPlus Version 2.0-10 (2016-06-08)}{
  \itemize{
    \item{Added \code{ggplotFuncs} argument to predictionplot.asreml and pred.present.asreml.}
    \item{Added \code{bound.test.parameters} to reml.lrt.asreml, choose.model.asartests, testranfix.asrtests, testrcov.asrtests. testwapran.asrtests.}
 }
}
\section{Changes in asremlPlus Version 2.0-9 (2016-05-15)}{
  \itemize{
    \item{Revised some documentation.}
 }
}
\section{Changes in asremlPlus Version 2.0-8 (2016-04-18)}{
  \itemize{
    \item{Fixed a minor bug in rmboundary.asrtests.}
 }
}
\section{Changes in asremlPlus Version 2.0-7 (2016-04-06)}{
  \itemize{
    \item{Fixed some minor bugs in testswapran.asrtests.}
 }
}
\section{Changes in asremlPlus Version 2.0-6 (2016-02-16)}{
  \itemize{
    \item{Fixed some minor bugs.}
 }
}
\section{Changes in asremlPlus Version 2.0-5 (2016-01-24)}{
  \itemize{
    \item{Added WaterRunoff.dat data set and modified examples based on this data set so that the examples would run.}
    \item{Fixed some minor bugs.}
 }
}
\section{Changes in asremlPlus Version 2.0-4 (2015-11-12)}{
  \itemize{
    \item{Added means argument to variofaces.asreml.}
 }
}
\section{Changes in asremlPlus Version 2.0-3 (2015-11-11)}{
  \itemize{
    \item{Modifications to predictionplot.asreml and pred.present.asreml as a result of changes made to ggplot2 version 1.0.1.}
 }
}
\section{Changes in asremlPlus Version 2.0-2 (2015-09-22)}{
  \itemize{
    \item{Added label argument to addrm.terms.asrtests.}
    \item{Added new function newrcov.asrtests.}
 }
}
\section{Changes in asremlPlus Version 2.0-1 (2015-09-14)}{
  \itemize{
    \item{Added graphics.device argument to variofaces.asreml and changed default for this argument 
          to \code{NULL} in pred.present.asreml and predictionplot.asreml.}
 }
}
\section{Changes in asremlPlus Version 2.0-0 (2015-06-24)}{
  \itemize{
    \item The following functions were deprecated and replaced with functions that have the suffix 
          asrtests instead of asreml, to signify that they operate on asrtests objects: 
          addrm.terms.asreml, choose.model.asreml, recalc.wald.tab.asreml, rmboundary.asreml, 
          sig.devn.reparam.asreml, testranfix.asreml, testrcov.asreml and testswapran.asreml.
    \item The following functions were deprecated info.crit and reml.lrt; they are replaced 
          by info.crit.asreml and reml.lrt.asreml, repsectively.
    \item NEWS and documentation of deprecated functions added.
  }
}
\section{Changes in asremlPlus Version 1.10-16 (2015-06-06)}{
  \itemize{
    \item First version placed on CRAN.
  }
}

