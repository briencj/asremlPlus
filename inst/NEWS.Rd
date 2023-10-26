\name{news}
\title{News for Package asremlPlus}
\section{asremlPlus Version 4.4.17 (2023-10-26}{
  \itemize{
  \item Fix bug in investigating nugget variance for correlation models when there are \code{sections}.
  \item Further improve the handling of situations in which the residual variance or the spatial variance associated with a correlation term becomes singular or bound in fitting a spatial correlation model using addSpatialModel.asrtests, addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests.
  \item Add \code{allow.corrsJointFit} argument as a device for dealing with fits that hang when attempting to fit simultaneously both correlations for the two grid dimesnions. 
  \item Internal recoding.
}}
\section{asremlPlus Version 4.4.16 (2023-10-18}{
  \itemize{
  \item Added the \code{maxit} option for \code{asreml} to addSpatialModel.asrtests, addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests, setting its default to 30 for these \code{functions} to increase the chances that the models will converge.
  \item Extensively revise the code to improve the handling of situations in which the residual variance or the spatial variance associated with a correlation term becomes singular or bound in fitting a spatial correlation model using addSpatialModel.asrtests, addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests.
  \item Revise checking for nugget variance when there are sections and the residual has a separate terms for each section.
  \item Added a basic check for separate random and residual terms for sections and, if there are none, issue a warning that better fitting can result if they are included.
}}
\section{asremlPlus Version 4.4.15 (2023-10-02}{
  \itemize{
  \item Make the default number of cores (\code{ncores}) for parallel processing in bootREMLRT.asreml, simulate.asreml and variofaces.asreml.
  \item Add the function convAsremlobj.asreml that recreates an existing \code{asreml} object using the currently loaded \code{asreml} version so that the object is compatible with this version.
  \item Made the following modifications to addSpatialModel.asrtests, addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests:
  \itemize{
  \item Tweaked the code underlying these functions for fitting local spatial variation models to make it more robust.
  \item Added reporting the optimal rotation angles to the test.summary when \code{rotateX} is \code{TRUE}. Also, added theta(s) as an attribute to the \code{asreml.object}.
  \item Added the \code{which.rotacriterion} so that one could choose between the \code{deviance}, \code{likelihood}, \code{AIC} and \code{BIC} as the criterion to use in determining the optimal roation angle for rotating the eigenvectors of the penalty matrix.
  \item Add the capacity to parallel process, for each column angle, the search of the row angles. This option is only available for the \code{asreml.option} set to \code{grp}.
  \item Removed the need to precede a call to these functions with a call makeTPPSplineMats.data.frame when \code{asreml.option} is set to \code{mbf}.
  \item Fix a bug in fitting natural cubic smoothing spines (\code{TPNCSS}) using these functions where only the first of the \code{sections} had splines fitted to it.
  \item Change from ordinal numbers to levels in labelling the fitting of correlation \code{corr} spatial models when using these functions.
}
  \item Modify makeTPPSplineMats.data.frame so that the default enviroment of the data.frames for the \code{mbf} option of \code{asreml.opt} is that of the function that calls makeTPPSplineMats.data.frame. Also, for \code{mbf.env} set to \code{NULL}, do not assign the \code{mbf} data.frames; just include them in the \code{list} that is returned.
  \item Make the \code{data.plus} component of the \code{list} produced by makeTPPSplineMats.data.frame when \code{asreml.option} is set to \code{mbf} equivalent to that produced when \code{asreml.option} is set to \code{grp}.
  \item Corrected bug in detecting \code{Fixed (F)} variance components in rmboundary.asrtests.
  \item Added to newfit.asreml, a test for bound and singular terms in a new model when \code{update} is \code{TRUE} and, if there are any, added code to refit the model as if \code{update} was \code{FALSE}, to see if any bound terms become unbound terms.
  \item Clarified the dcoumentation of the \code{update} argument in all functions that include it.
}}
\section{asremlPlus Version 4.4.14 (2023-08-28)}{
  \itemize{
  \item Minor internal code changes.
}}
\section{asremlPlus Version 4.4.13 (2023-08-24)}{
  \itemize{
  \item Add arguments \code{rotateX} and \code{theta} to makeTPPSplineMats.data.frame so that the eigenvectors of the penalty matrix used in generating the design matrices for the random terms could be rotated by the angles in \code{theta}.
  \item Add code to carry out a grid search to find the optimal angle by which to rotate the eigenvectors of the penalty matric when fitting P-splines. Add arguments \code{rotateX} and \code{ngridlines} to addSpatialModel.asrtests, addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests to allow the use to request that grid seach be made for the optimal ritation angle and to control the size of the grid searched.
  \item Fix a bug in fitting a correlation model that occurs when addSpatialModel is used to fit the model.
  \item There remains an svd bug with the oneMKL libraries version 2023.1.0 that is causing problems for \code{IClikelhood} set to \code{"full"}.
}}
\section{asremlPlus Version 4.4.12 (2023-08-03)}{
  \itemize{
  \item Modifications to deal with classhes between multithreading and external linear algebra packages  under Windows. The environmental variable OMP_NUM_THREADS needs to be set to 1 if parallel processing is to be undertaken on a Windows machine. This can be done for individual tasks by using \code{Sys.setenv("OMP_NUM_THREADS" = 1)} prior to running the task; otherwise, the environmental variable can be set for all tasks in the Windows system via \code{Advanced system settings} > \code{Envirnomental Variables}.
  \item There remains an svd bug with the oneMKL libraries version 2023.1.0 that is causing problems for \code{IClikelhood} set to \code{"full"}.
}}
\section{asremlPlus Version 4.4.11 (2023-07-19)}{
  \itemize{
  \item This version has known issues when the Intel MKL libraries are installed in the R installation directories.
  \item Modified makeTPPSplineMats.data.frame to implement the use of mbf in fitting P-spline models for local spatial variation.
}}
\section{asremlPlus Version 4.4.10 (2022-12-28)}{
  \itemize{
  \item Modifications for ASReml-R version 4.2. Those that may impact upon users include the following:
     \itemize{
        \item The names of the columns of the design matrix has changed, as has the order of the columns.
        \item character codes are now stored in both \code{vparameters.type} and \code{vparameters.con},
              making vpt.char and vpc.char obsolete and so no longer available.
        \item For predict, the intercept-only classify and the name of the column in the 
             \code{predictions.freame}, has changed from \code{(Intercept)} to \code{Intercept}.
        \item The designation of estmates on the boder of the parameter space has changed from 
              fixed (\code{F}) to (\code{B}) e.g. correlations.
        \item When the levels of a factor are being specified in  terms with an \code{at}, as 
              opposed to an ordinal number for the position of the level in the levels list, then 
              the levels must be enclosed in single or double quotes. Further, the names of the terms 
              will nclude single or double quotes.}
  }}
\section{asremlPlus Version 4.3.56 (2023-06-20)}{
  \itemize{
  \item Correct issue with alldiffs attributes not being retained when subsetting using subset.alldiffs.
}}
\section{asremlPlus Version 4.3.55 (2023-06-13)}{
  \itemize{
  \item Remove bloated system-generated testthat-problems.rds file so that the tarball to be submitted is less than 5 Mb.
}}
\section{asremlPlus Version 4.3.54 (2023-06-12)}{
  \itemize{
  \item Fix problem with deprecated function setvarianceterms.asreml when submitted to CRAN.
}}
\section{asremlPlus Version 4.3.53 (2023-06-12)}{
  \itemize{
  \item Add new vignette \code{WheatSpatialModels} that illustrates the fitting of models for local spatial variation.
  \item Generalize correlation spatial models to allow the independent specification of any asreml variance model for rows and columns. Improve the fitting of spatial correlation models, including changing the default correlation model to \code{ar1} and postpone removing bundary terms until after all components of the correlarion model have been fitted.
  \item Revised the spatial modelling to resolve several bugs associated with them when the residual specifies heterogeneous variance using eith \code{idh} or \code{dsum}.
  \item Repurpose \code{row.factor} and \code{col.factor} arguments of \code{addSpatialModel.asrtests}, \code{addSpatialModelOnIC.asrtests} and \code{chooseSpatialModelOnIC.asrtests} to specify the row and column factors for correlation model where required and add the arguments \code{dropRowterm} and \code{dropColterm} to specify the main effect terms to be dropped in fitting spline models.
  \item Modified changeTerms.asrtests and changeModelOnIC.asrtests so that the \code{bounds} and \code{initial values} of terms specified using \code{set.terms} without changing any terms in the model.
  \item Fix problem detected during CRAN submission that found setvarianceterms.asreml was apparently an unregistered exported generic.
}}
\section{asremlPlus Version 4.3.52 (2023-05-18)}{
  \itemize{
  \item Fix bug in identifying residual terms for ASReml-R v. 4.x.
  \item Fix bug in rmboundary.asrtests when removing a bound term results in no random terms.
}}
\section{asremlPlus Version 4.3.51 (2023-04-02)}{
  \itemize{
  \item Fix bug in fitting P-splines when difforder is one.
  \item Change infoCriteria to use the deviance in calculating information criteria for \code{family} set to asr_binomial or asr_poisson.
  \item Add the argument \code{transform.function} to addBacktransforms.alldiffs, allDifferences.data.frame, predictPlus.asreml and predictPresent.asreml so that \code{predictions} and \code{error.intervals} on the scale of the \code{transform.function} can be backtransformed.
}}
\section{asremlPlus Version 4.3.50 (2023-03-9)}{
  \itemize{
  \item Add \code{classify} argument to as.predictions.frame, checking that variables in the  \code{classify} are in the \code{predictions.frame} and conversion of any variables in the \code{classify} that are not factor to factors.
  \item Update ggplot2 functions to use tidy evaluation.
  \item Fix internal anomaly causing a Note on CRAN.
}}
\section{asremlPlus Version 4.3.49 (2023-03-3)}{
  \itemize{
  \item Add the argument \code{useIntervals} to plotLSDerror.alldiffs.
  \item Fix a bug in plotLSDerrors that resulted in not being able to determine which LSDs in the LSD component apply to the pairwise comparisons. It resulted in all comparisons being shown as \code{n.a.}.
}}
\section{asremlPlus Version 4.3.48 (2023-02-18)}{
  \itemize{
  \item Rename makeTPSPlineXZMats.data.frame to makeTPPSplineMats.data.frame.
}}
\section{asremlPlus Version 4.3.47 (2023-02-17)}{
  \itemize{
  \item Tidy up code for trySpatial = "none".
  \item Revise makeTPSPlineXZMats.data.frame to cover \code{asreml.opt} set to \code{grp}.
  \item Fix bug in plotLSDerrors when LSDs have beeen computed for  \code{LSDtype} set to \code{factor.combinations}.
}}
\section{asremlPlus Version 4.3.46 (2023-02-05)}{
  \itemize{
  \item Fix bug in testing whether an ar1 fit has been achieved in addSpatialModel.asrtests and chooseSpatialModelOnIC.asrtests.
}}
\section{asremlPlus Version 4.3.45 (2023-01-24)}{
  \itemize{
  \item Modify system tests to cater for for M1Mac numerical accuracy.
  \item Fix bug in which the setting of IClikelihood was not being passed to infoCriteria calls in chooseSpatialModelOnIC.asrtests, resulting in always reporting the \code{REML} likelihood.
  \item Change that \code{spatial.IC} component of the \code{list} returned by chooseSpatialModelOnIC.asrtests to contains the information criteria for each spatial model that is fitted. Also, rename the components \code{best} and \code{best.AIC} to \code{best.spatial.mod} and \code{best.spatial.IC}.
  \item Fix bug in reporting the likelihood when a TPNCSS model does not fit.
  \item Add function addSpatialModel.asrtests that adds a spatial model to a supplied model without asessing the fit of the model.
}}
\section{asremlPlus Version 4.3.44 (2022-12-30)}{
  \itemize{
  \item Add arguments \code{nestorder}, \code{degree} and \code{difforder} to functions addSpatialModelOnIC.asrtests and  makeTPSPlineXZMats.data.frame, as well as argument \code{nestorder} to chooseSpatialModelOnIC.asrtests. Add the option \code{TPLV} to the argument \code{trySpatial} from chooseSpatialModelOnIC.asrtests to allow a linear spatial model with first-order penalty to be fitted.
  \item Fix bug in testing for fixed correlations that meant they were not always detected.
  \item Fix error in changeModelsonIC.asrtests in the action when the model is \code{Unswapped} - it was being reported as  \code{Changed random}. 
  \item Stop crashing of changeTerms.asrtests and ChangeModelOnIC.asrtests when rmboundary.asrtests unable to remove singular terms.
}}
\section{asremlPlus Version 4.3.43 (2022-11-30)}{
  \itemize{
  \item Fix bugs in choosing the best spatial model in chooseSpatialModelOnIC.asrtests. Add components \code{best.spatial} and \code{best.spatial.AIC} to the value returned by chooseSpatialModelOnIC.asrtests.
  \item Add linear-row-by-linear-column-covariates interaction to the model fitted for TPNCSS option of addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests.
}}
\section{asremlPlus Version 4.3.42 (2022-11-23)}{
  \itemize{
  \item Add functions makeTPSPlineZMats.data.frame, addSpatialModelOnIC.asrtests and chooseSpatialModelOnIC.asrtests.
}}
\section{asremlPlus Version 4.3.41 (2022-11-16)}{
  \itemize{
  \item Revise test for fixed correlation so that allow.fixedcorrelation recognises any correlation term, rather than being restricted to ar1.
}}
\section{asremlPlus Version 4.3.40 (2022-11-16)}{
  \itemize{
  \item Compress pdf files and resubmit to CRAN.
}}
\section{asremlPlus Version 4.3.39 (2022-11-13)}{
  \itemize{
  \item Rebuild for R 4.2.2 and submission to CRAN.
}}
\section{asremlPlus Version 4.3.38 (2022-10-24)}{
  \itemize{
  \item Added function pickLSDstatistics.alldiffs.Rd.
  \item Revise documentation for the arguments \code{LSDtype}, \code{LSDsupplied}, \code{LSDby} and \code{LSDstatistic }.
}}
\section{asremlPlus Version 4.3.37 (2022-09-19)}{
  \itemize{
  \item Fix bugs in function newfit.asreml that involves making changes to models that involve  the \code{at} function. This bug affected the functions changeTerms.asrtests, changeModelOnIC.asrtests, rmboundary.asrtests, testswapran.asrtests and testranfix.asrtests. In removing a term that involve a single lewvel of a factor in an at function, it is now possible to specify either the level or the level index.
  \item Modify rmboundary.asrtests to allow the removal of random term for a single level of an \code{at} function when the levels are specified singly in the random model.
  \item Add function addto.test.summary.
}}
\section{asremlPlus Version 4.3.36 (2022-08-17)}{
  \itemize{
  \item Minor correction for resubmitting to CRAN.
}}
\section{asremlPlus Version 4.3.35 (2022-08-16)}{
  \itemize{
  \item Rebuild under R 4.2.1 for new submission.
}}
\section{asremlPlus Version 4.3-34 (2022-04-20)}{
  \itemize{
  \item Add functions plotLSDerrors.alldiffs and plotLSDerrors.data.frame.
  \item Replace the use of the package \code{reshape} with the use of \code{reshape2}.
  \item Replace \code{X1} and \code{X2} with \code{Rows} and \code{Columns} as the default names for the the first two columns of the \code{data.frame} produced by plotPvalues.alldiffs, plotLSDs.alldiffs and plotLSDerrors.alldiffs and provided as input to plotPvalues.data.frame, plotLSDs.data.frame and plotLSDerrors.data.frame.
  \item Change plotPvalues.alldiffs and plotLSDs.alldiffs to return a two component list, one component of which is a \code{data.frame} that contains the statistics to be plotted and the other of which is a \code{list} containing the plots that have be produced as \code{ggplot} objects.
}}
\section{asremlPlus Version 4.3-33 (2022-04-09)}{
  \itemize{
  \item Add an error if all residual variances are bound.
}}
\section{asremlPlus Version 4.3-32 (2021-12-21)}{
  \itemize{
  \item Fix bug in chooseModelOnIC result from a change in as.vector for data.frames and lists in R 4.2.0.
}}
\section{asremlPlus Version 4.3-31 (2021-12-06)}{
  \itemize{
  \item Rebuild under R 4.1.2 for new submission.
}}
\section{asremlPlus Version 4.3-30 (2021-11-22)}{
  \itemize{
  \item Put in a trap for fixed-effects variances being inestimable when calculating the full likelihood.
}}
\section{asremlPlus Version 4.3-29 (2021-11-06)}{
  \itemize{
  \item Modify linTransform.alldiffs so that an intercept-only submodel can be fitted.
  \item Revise the manual entry for exploreLSDs.alldiffs and for pairdiffsTransform.alldiffs.
}}
\section{asremlPlus Version 4.3-28 (2021-10-26)}{
  \itemize{
  \item Fix bug in newfit.asreml that re-orders the fixed model because of R limitations.
}}
\section{asremlPlus Version 4.3-27 (2021-10-24)}{
  \itemize{
  \item Fix problem in linTransform.alldiffs when the product of the numbers of levels in the classiify factors is large and not all combinations of the levels are observed.
  \item Relax the requirement for \code{LSDby} variable to be in the \code{classify} to beng that they must be in the \code{predictions} component of the \code{alldiffs.obj}.
}}
\section{asremlPlus Version 4.3-26 (2021-10-21)}{
  \itemize{
  \item Fix bug in false.pos and false.neg for \code{LSDtype} set to \code{supplied} when \code{LSDby} is not \code{NULL}.
  \item Add quantiles 25 (\code{q25}) and 75 (\code{q75}) to the values that \code{LSDstatistic} can take and the exploreLSDs.alldiffs investigates.
  \item Allow for multiple \code{LSDstatistic} values when \code{LSDtype} is \code{factor.combination} and \code{LSDby} is set.
}}
\section{asremlPlus Version 4.3-25 (2021-10-15)}{
  \itemize{
  \item Add \code{sig.size}, \code{sig.face} and \code{sig.colour} arguments to plotPvalues.data.frame and plotPvalues.alldiffs.
}}
\section{asremlPlus Version 4.3-24 (2021-10-14)}{
  \itemize{
  \item Make sure that when \code{predicted.values} and the associated \code{vcov} and \code{sed} in an \code{\link{alldiffs.object}} have values that are close to zero as a result of the limitations of real number calculations, they are set to zero so as to avoid further calculation problems.
  \item Add \code{zero.tolerance} argument to linTransform.alldiffs.
}}
\section{asremlPlus Version 4.3-23 (2021-10-03)}{
  \itemize{
  \item Modify the test for fixed correlation from a test for \code{Fixed} (\code{F}) to cover also correlations that are \code{bound} (\code{B}) or \code{singular} (\code{S}).
  \item Add tests on incoming fits for fixed correlation.
}}
\section{asremlPlus Version 4.3-22 (2021-10-01)}{
  \itemize{
  \item Support is no longer provided in asremlPlus for versions of ASReml-R prior to version 4.0.
  \item Add a new argument \code{allow.fixedcorrelation} to changeModelOnIC.alldiffs, newfit.asreml, chooseModel.asrtests, changeTerms.asrtests, testranfix.asrtests, testresidual.asrtests, testswapran.asrtests and reparamSigDevn.asrtests so that a model with a fixed correlation parameter may be rejected.
}}
\section{asremlPlus Version 4.3-21 (2021-09-23)}{
  \itemize{
  \item Fix a bug in pairdiffsTransform.alldiffs where the factors indexing the predictions were not generated correctly.
  \item In changeModelOnIC.alldiffs, add detection of absent terms to be dropped and do not change the model if one or more are detected.
}}
\section{asremlPlus Version 4.3-20 (2021-09-03)}{
  \itemize{
  \item Fix bug in the calculation of false positives and negatives in exploreLSDs.alldiffs for \code{LSDtype} set to \code{overall}.
  \item Change default for \code{plotHistogram} in exploreLSDs.alldiffs to \code{FALSE}.
}}
\section{asremlPlus Version 4.3-19 (2021-09-01)}{
  \itemize{
  \item Rebuild under R 4.1.1 for new submission.
}}
\section{asremlPlus Version 4.3-18 (2021-08-31)}{
  \itemize{
  \item Change number of comparisons symbol in an \code{LSD.frame} from \code{n} to \code{c}.
}}
\section{asremlPlus Version 4.3-17 (2021-08-30)}{
  \itemize{
  \item Add False positives and negatives counts to LSD statistics.
}}
\section{asremlPlus Version 4.3-16 (2021-08-24)}{
  \itemize{
  \item Fix bug in redoErrorIntervals when only one factor in \code{LSDby}.
  \item Add function print.LSDdata.
}}
\section{asremlPlus Version 4.3-15 (2021-08-19)}{
  \itemize{
  \item Fixed bug in \code{halfLeastSignificant.limit} when \code{LSDtype} is \code{factor.combination} and the \code{assignedLSD} is zero to now ensure that the limit is \code{NA}.
  \item Reinstated passing of parameters from linTransform.alldiffs to redoErrorIntervals.alldiffs.
}}
\section{asremlPlus Version 4.3-14 (2021-08-18)}{
  \itemize{
  \item Fix newly introduced  bug for \code{LSDtype = "halfLeastSignificant"} when there is not a numeric covariate involved in the \code{LSDby} argument.
  \item Stop extra frequency plot produced by exploreLSDs.
  \item Minor internal modification.
}}
\section{asremlPlus Version 4.3-13 (2021-08-16)}{
  \itemize{
  \item Fix bug for \code{LSDtype = "halfLeastSignificant"} when there is a numeric covariate involved in the \code{LSDby} argument.
  \item Minor internal modification.
}}
\section{asremlPlus Version 4.3-12 (2021-08-11)}{
  \itemize{
  \item Fix bug in pairdiffsTransform when only two factors in \code{classify}.
  \item Add function subset.list.
  \item Add functions exploreLSDs, plotLSDs.alldiffs, plotLSDs.matrix and plotLSDs.data.frame.
  \item Incorporate arguments \code{retain.zerLSDs} and \code{zero.tolerance} into AllDifferences.data.frame, redoErrotIntervals.alldiffs amd exploreLSDs.alldiffs to allow control of the retention of the zero LSD values when calculating LSD summaries.
}}
\section{asremlPlus Version 4.3-11 (2021-07-28)}{
  \itemize{
  \item Deprecate the function asrtests, which has been replaced by the function as.asrtests.
  \item There has been a substantial revision of the calculation and saving of the LSD component of an \code{alldiffs} object. 
    \itemize{
  \item The argument \code{meanLSD.type} has been changed to \code{LSDtype}. All functions that are S3 methods for \code{alldiffs} objects will check for the \code{meanLSD.type} attribute and, if found, change it to \code{LSDtype}.
  \item The option \code{supplied} has been added to the argument \code{LSDtype} so that the user can supply LSD values to replace those in the \code{assignedLSD} column of the \code{LSD} component.
  \item The two new arguments \code{LSDstatistic} and \code{LSDaccuracy} have been added. The argument \code{LSDstatistic} allows a choice of \code{minimum}, \code{mean}, \code{median} or \code{maximum}. The values of the nominated LSD statistics are stored in a column of the LSD component named \code{assignedLSD}. The argument \code{LSDaccuracy} allows the specification of a statistic summarizing the differences between the the LSDs and the assigned LSD value for each row of the LSD component. The values of the chosen measure are stored in a column of the \code{LSD} component named \code{accuracyLSD}. These two arguments have been added to the functions allDifferences.data.frame, as.alldiffs.object, as.alldiffs, linTransform.alldiffs, pairdiffsTransform, predictPlus.asreml, predictPresent.asreml, ratioTransform.alldiffs, recalcLSD.alldiffs and redoErrorIntervals.alldiffs.
  \item The \code{LSD} component is always calculated and stored as a component when \code{pairwise} is \code{TRUE}.
  }
  \item A new argument \code{accuracy.threshold} has been introduced to control the addition of a column names \code{LSDwarning} to the \code{predictions} of an \code{alldiffs} object. THis column indicates whether or not the an \code{halfLeastSignificant} \code{error.interval} is reliable as determined using the supplied value of  \code{accuracy.threshold}.
  \item The attributes of an alldiffs object  and of the \code{predictions} component have been modified. The arguments \code{LSDtype}, \code{LSDstatistic} and \code{LSDby}, if set, are attributes of an \code{alldiffs} object. They are also attributes of a \code{prediction.frame} and the \code{backtransforms} component of an \code{alldiffs} object, as is \code{LSDvalues}, when \code{halfLeastSignificant} \code{error.intervals} are specified. The \code{LSDvalues} attribute has the values used in calculating half-LSD intervals. The \code{meanLSD} attribute has been removed.
  \item A separate entry for \code{LSD.frame} has been introduced into the manual to dedscribe the \code{data.frame} that is stored in the \code{LSD} component of an \code{alldiffs} object.
}}
\section{asremlPlus Version 4.3-10 (2021-07-5)}{
  \itemize{
  \item Fix bug in facRecast that meant that the argument \code{levels.order} was being ignored.
  \item Fix a bug in facRecast where only the predictions component was appropriately re-ordered.
  \item Add new functions ratioTransform.alldiffs and pairdiffsTransform.
  \item Add new method to sort for a \code{predictions.frame}: function sort.predictions.frame.
  \item Change the argument \code{sortWithinVals} to \code{sortParallelToCombo} in sort.alldiffs.
  \item Add the new argument \code{sortNestinFactor} to sort.alldiffs.
  \item Add new function getTestEntry.asrtests.
  \item Fix bug in infoCriteria.asreml that resulted in values of \code{Inf} for the \code{full} likelihood.
}}
\section{asremlPlus Version 4.2-36 (2021-05-16)}{
  \itemize{
  \item Deprecate facRecode.alldiffs because it only make sense to redo the current levels or labels.
  \item Tidy up function facRecast.alldiffs.
}}
\section{asremlPlus Version 4.2-35 (2021-05-16)}{
  \itemize{
  \item Add new function facRecast.alldiffs that supercedes facRecode.alldiffs.
  \item Revise sort.alldiffs to use facRecast.alldiffs to reorder the alldiffs components. This fixed a bug when not all levels of the \code{saortFactor} are present in the subset on which the sort is based.
  \item Add argument \code{alpha} to plotPvalues.data.frame and plotPvalues.alldiffs.
}}
\section{asremlPlus Version 4.2-34 (2021-05-12)}{
  \itemize{
    \item Add argument \code{axis.text.size} to plotPvalues.data.frame and plotPvalues.alldiffs.
}}
\section{asremlPlus Version 4.2-33 (2021-05-06)}{
  \itemize{
    \item Modify testresidual to not change the variance model when there is no difference in the variance DF between the models being compared.
    \item In as.asrtests, add variance DF to \code{test.summary} when forming an \code{asrtests} object.
}}
\section{asremlPlus Version 4.2-32 (2021-03-19)}{
  \itemize{
    \item Rebuild vignerttes.
    \item Modify more links to asreml.
}}
\section{asremlPlus Version 4.2-31 (2021-03-19)}{
  \itemize{
    \item Modify link to asreml.
}}
\section{asremlPlus Version 4.2-30 (2021-03-18)}{
  \itemize{
    \item Revise doi specification in a link to a reference in the manual.
    \item Increment version number for resubmission.
}}
\section{asremlPlus Version 4.2-29 (2021-03-18)}{
  \itemize{
    \item Increment version number for resubmission.
}}
\section{asremlPlus Version 4.2-28 (2021-03-17)}{
  \itemize{
    \item Recompile for CRAN submission.
}}
\section{asremlPlus Version 4.2-27 (2021-03-13)}{
  \itemize{
    \item Fix bugs in redoErrorIntervals.alldiff concerning whether or not meanLSD.type and LSDby settings result in the LSD component being recalculated.
    \item Add rowname \code{overall} to the LSD component when \code{meanLSD.type} is set to \code{overall}.
}}
\section{asremlPlus Version 4.2-26 (2020-11-11)}{
  \itemize{
    \item Modify URLs for CRAN submission.
}}
\section{asremlPlus Version 4.2-25 (2020-11-11)}{
  \itemize{
    \item Try to neutralize problem that occurs with \code{colourise} in print.predictions.frame.
    \item Rebuild for R 4.0.3 and recompile vignettes for CRAN submission.
}}
\section{asremlPlus Version 4.2-24 (2020-09-25)}{
  \itemize{
    \item Include a \code{key} option in the \code{which} that produces \code{vparametersummary} rather than the full \code{asremlsummary}.
    \item Remove attributes from formuale when executing printFormulae.
    \item Add the ChickpeaEnd data set.
}}
\section{asremlPlus Version 4.2-23 (2020-09-07)}{
  \itemize{
    \item Modify code so that error intervals are set to NA when the LSDs are zero.
}}
\section{asremlPlus Version 4.2-22 (2020-09-04)}{
  \itemize{
    \item Modify code for calculating LSDs so that all zero LSDs are allowed; this is needed when a linear transformation has resulted in predictions that are equal.
}}
\section{asremlPlus Version 4.2-21 (2020-07-11)}{
  \itemize{
    \item Rebuild the vignettes prior to CRAN submission.
    \item Add \code{envir} argument to getFormulae.asreml and printFormulae.asreml to make them more robust to where the objects in the \code{asreml.obj} have been defined.
}}
\section{asremlPlus Version 4.2-20 (2020-06-29)}{
  \itemize{
    \item Fix bug in getFormulae.asreml and printFormulae.asreml that occurs when formula(e) have been supplied to asreml in an object.
}}
\section{asremlPlus Version 4.2-19 (2020-04-02)}{
  \itemize{
    \item Change the Maintainer email address.
    \item Add ellipsis to as.asrtests call within an rmboundary.asrtests call so that arguments in the calling routine are passed through.
    \item Add a trap in recalcWaldTab.asrtests for a \code{NULL} \code{wald.tab} before dealing with \code{denDF} so that recalcWaldTab.asrtests throws a warning rather than an error.
    \item Modify recalcWaldTab.asrtests to deal with a Wald table that has Chi-square instead of F statistics (this happens when \code{denDF = "none"}.)
}}
\section{asremlPlus Version 4.2-18 (2020-03-16)}{
  \itemize{
    \item Increment version number for resubmission to CRAN.
}}
\section{asremlPlus Version 4.2-17 (2020-03-16)}{
  \itemize{
    \item Fix bug in changeTerms when both \code{addFixed} and \code{dropFixed} are \code{NULL}.
}}
\section{asremlPlus Version 4.2-16 (2020-02-25)}{
  \itemize{
    \item Revisions to the manual.
    \item Fix bug in forming test.summary.
}}
\section{asremlPlus Version 4.2-15 (2020-02-17)}{
  \itemize{
    \item Added the \code{IClikelihood} argument to chooseModel.asrtests, reparamSigDevn.asrtests, rmboundary.asrtests, testranfix.asrtests, testresidual.asrtests, and testswapran.asrtests.
    \item Fix bug in printFormulae when the formula is long.
    \item Various bug fixes associated with asreml-R version 3.
}}
\section{asremlPlus Version 4.2-14 (2020-02-14)}{
  \itemize{
    \item Fix bug in processing a formula that includes an \code{at} function.
}}
\section{asremlPlus Version 4.2-13 (2020-02-11)}{
  \itemize{
    \item Removed the \code{material.diff} argument and the \code{both} option from the \code{which.IC} argument of changeModelOnIC.
}}
\section{asremlPlus Version 4.2-12 (2020-02-09)}{
  \itemize{
    \item Change \code{likelihood} in infoCriteria to \code{IClikelihood} to make it consistent with other functions.
    \item Replace \code{tol.diff} with \code{material.diff} and set the default to 0.5 in changeModelOnIC.
    \item Revise all vignettes.
}}
\section{asremlPlus Version 4.2-11 (2020-02-08)}{
  \itemize{
    \item Add a new vignette for the use of information criteria with the wheat experiment.
    \item Fixed bugs in changeModelOnIC associated with \code{addFixed} and \code{dropFixed} and with dealing with unconverged models.
}}
\section{asremlPlus Version 4.2-10 (2020-02-07)}{
  \itemize{
    \item This version represents a major revision in that the facilities for employing information criteria (AIC and BIC) have been considerabley expanded. The infoCriteria method has been modified to add the possibility of calculating not only the the information criteria based on the REML likelihood, but those based on the full likelihood instead. There are now two infoCriteria methods, one for \code{asreml} objects and the other for \code{list}s of \code{asreml} objects. An alternative is to use changeTerms.asrtests wiith \code{IClikelihood} set to \code{REML} or \code{full} and the information criteria will be included in the The \code{test.summary} of the \code{asrtest} object. A new function changeModelOnIC has been included to allow model selection to be based on the information criteria, rather than on hypothesis tests. However, the two methods of selecting models can be interweaved.
    \item Add  argument \code{omit.columns} to print.test.summary that allows the specification of columns to be omitted from the printout.
    \item Various modifications to the internal code for the handling of the \code{test.summary} component of an \code{asrtest.object} to allow more flexibility in the make up of \code{test.summary}.
    \item Add function iterate.asrtests to execute additional iterations of the fit stored in an  \code{asreml.obj} within an \code{asrtests} object; the model itself cannot be changed.
    \item Add functions getFormulae.asreml and printFormulae.asreml to extract and print the \code{formula}e stored in the \code{call} of an \code{asreml} object.
    \item Add \code{trace} argument to newfit.asreml and activated it for asreml-R version 4.
    \item Made R >= 3.5.0 a requirement for asremlPlus.
}}
