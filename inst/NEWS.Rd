\name{news}
\title{News for Package asremlPlus}
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
    \item{Provided as.asrtests function, while retaining asrtests function as a pseudonym for as.asrtests. Added as.prediction.frame function.}
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

