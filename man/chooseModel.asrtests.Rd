\name{chooseModel.asrtests}
\alias{chooseModel.asrtests}
\title{Determines and records the set of significant terms using an \code{\link{asrtests.object}}, 
       taking into account the hierarchy or marginality relations of the terms.}
\description{Performs a series of hypothesis tests on a set of fixed and/or random terms taking 
  into account the marginality of terms. In particular, a term will not be tested if it is 
  marginal to (or nested in) one that is significant. For example, if A:B is significant, then 
  neither A nor B will be tested. For a random term, the term is removed from 
  the model fit, any boundary terms are removed using \code{\link{rmboundary.asrtests}} 
  and a REML likelihood ratio test is performed using \code{\link{REMLRT.asreml}}.  
  If it is not significant and \code{drop.ran.ns} is TRUE, the term is permanently removed 
  from the model. Note that if boundary terms are removed, the reduced model may not 
  be nested in the full model in which case the test is not valid. For fixed terms, 
  the Wald tests are performed and the p-value for the term obtained. If it is not 
  significant and \code{drop.fix.ns} is TRUE, the term is permanently removed 
  from the model. A row that records the outcome of a test is added to 
  \code{test.summary} for each term that is tested. }
\usage{
\method{chooseModel}{asrtests}(object, terms.marginality=NULL, 
            alpha = 0.05, allow.unconverged = TRUE, 
            allow.fixedcorrelation = TRUE,  
            checkboundaryonly = FALSE, drop.ran.ns=TRUE, 
            positive.zero = FALSE, bound.test.parameters = "none", 
            drop.fix.ns=FALSE, denDF = "numeric",  dDF.na = "none", 
            dDF.values = NULL, trace = FALSE, update = TRUE, 
            set.terms = NULL, ignore.suffices = TRUE, 
            bounds = "P", initial.values = NA, 
            IClikelihood = "none", ...)
}
\arguments{
 \item{object}{an \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{terms.marginality}{A square matrix of ones and zeros with row and column names 
                being the names of the terms to be tested. 
               The names of fixed terms must match those in the \code{wald.tab} component of the 
               \code{object}, while the names of random terms must match those 
               in the \code{vparameters} component of the 
               \code{asreml.obj} component in the \code{object}.
               The diagonal elements of the matrix should be one, indicating 
               that a term is marginal to itself. Elements should be one if the row 
               term is marginal to the column term. All other elements should be zero.}
 \item{alpha}{The significance level for the test.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
                      even when it does not converge. If \code{FALSE} and a fit when 
                      a term is removed does not converge, the term will not be removed.}
 \item{allow.fixedcorrelation}{A \code{logical} indicating whether to accept a new model 
              even when it contains correlations in the model whose values have been 
              designated as fixed, bound or singular. If \code{FALSE} and the new model 
              contains correlations whose values have not been able to be estimated, 
              the supplied \code{asrtests.obj} is returned. The fit in the 
              \code{asreml.obj} component of the supplied \code{asrtests.obj} will 
              also be tested and a warning issued if both fixed correlations are found 
              in it and \code{allow.fixedcorrelation} is \code{FALSE}.}
 \item{checkboundaryonly}{If \code{TRUE} then boundary and singular terms are not removed by 
              \code{\link{rmboundary.asrtests}}; a warning is issued instead.}
 \item{drop.ran.ns}{A logical indicating whether to drop nonsignificant random 
                    terms from the model.}
 \item{positive.zero}{Indicates whether the hypothesized values for the 
                      variance components being tested are on the boundary 
                      of the parameter space. For example, this is true 
                      for positively-constrained variance components that, 
                      under the reduced model, are zero. This argument does 
                      not need to be set if \code{bound.test.parameters} is set.}
 \item{bound.test.parameters}{Indicates whether for the variance components 
                      being tested, at least some of the hypothesized values
                      are on the boundary of the parameter space. 
                      The possibilities are \code{"none"}, \code{"onlybound"} 
                      and \code{"one-and-one"}. The default is \code{"none"}, 
                      although if it is set to \code{"none"} and 
                      \code{positive.zero} is \code{TRUE} then 
                      \code{bound.test.parameters} is taken to be \code{"onlybound"}. 
                      When \code{bound.test.parameters} is set to 
                      \code{"one-and-one"}, it signifies that there are two 
                      parameters being tested, one of which is bound and the 
                      other is not. For example, the latter is true for testing 
                      a covariance and a positively-constrained variance component 
                      that, under the reduced model, are zero.}
 \item{drop.fix.ns}{A logical indicating whether to drop a fixed  
                   term from the model when it is nonsignificant }
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to automatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{dDF.na}{The method to use to obtain substitute denominator degrees of freedom. 
              when the numeric or algebraic methods produce an \code{NA}. If 
              \code{dDF.na = "none"}, no substitute denominator degrees of freedom 
              are employed; if \code{dDF.na = "residual"}, 
              the residual degrees of freedom from \code{asreml.obj$nedf} are used; 
              if \code{dDF.na = "maximum"}, the maximum of those denDF that are available, 
              excluding that for the Intercept, is used; if all denDF are \code{NA}, 
              \code{asreml.obj$nedf} is used. If \code{dDF.na = "supplied"}, 
              a \code{vector} of values for the denominator degrees of freedom is to 
              be supplied in \code{dDF.values}. Any other setting is ignored and 
              a warning message produced. Generally, substituting these degrees of freedom is 
              anticonservative in that it is likely that the degrees of freedom used will be
              too large.}
 \item{dDF.values}{A \code{vector} of values to be used when \code{dDF.na = "supplied"}. 
              Its values will be used when \code{denDF} in a test for a fixed effect
              is \code{NA}. This vector must be the same length as the number of 
              fixed terms, including (Intercept) whose value could be \code{NA}.}
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{update}{If \code{TRUE}, and \code{set.terms} is \code{NULL}, then 
              \code{\link{newfit.asreml}} is called to fit the model to be tested, 
              using the values of the variance parameters stored in 
              the \code{asreml.object}, that is stored in \code{asrtests.obj}, as starting values. 
              If \code{FALSE} or \code{set.terms} is not \code{NULL}, then 
              \code{\link{newfit.asreml}} will not use the stored variance parameter 
              values as starting values when fitting the new model, the only 
              modifications being (i) to the terms in the fixed and random models 
              corresponding to terms in \code{terms.marginality} and 
              (ii) those specified via \code{\dots}.}
 \item{set.terms}{A character vector specifying the terms that are to have 
                  bounds and/or initial values set prior to fitting. 
                  The names must match those in the \code{vparameters} component of the 
                  \code{asreml.obj} component in the new \code{\link{asrtests.object}}.}
 \item{ignore.suffices}{A logical vector specifying whether the suffices of the 
                    \code{asreml}-assigned names of the variance terms (i.e. the 
                    information to the right of an "!", other than "R!") is to 
                    be ignored in matching elements of \code{terms}. 
                    If \code{TRUE} for an element of \code{terms}, the suffices 
                    are stripped from the \code{asreml}-assigned names. 
                    If \code{FALSE} for an element of \code{terms}, the element 
                    must exactly match an \code{asreml}-assigned name for a 
                    variance term. This vector must be of length one or the 
                    same length as \code{terms}. If it is of length one then 
                    the same action is applied to the \code{asreml}-assigned 
                    suffices for all the terms in \code{terms}.}
 \item{bounds}{A \code{\link{character}} vector specifying the bounds to be applied 
                    to the terms specified in \code{set.terms}. This vector 
                    must be of length one or the same length as \code{set.terms}. 
                    If it is of length one then the same constraint is 
                    applied to all the terms in \code{set.terms}. 
                    If any of the bounds are equal to NA then they are 
                    left unchanged for those terms.}
 \item{initial.values}{A character vector specifying the initial values for  
                    the terms specified in \code{terms}. This vector 
                    must be of length one or the same length as \code{terms}. 
                    If it is of length one then the same initial value is 
                    applied to all the terms in \code{terms}. 
                    If any of the initial.values are equal to NA then they are 
                    left unchanged for those terms.}
 \item{IClikelihood}{A \code{character} that controls both the occurrence and the type 
                    of likelihood for information criterion in the \code{test.summary} 
                    of the new \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML}, then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full}, 
                    then the AIC and BIC based on the full likelihood, evaluated using
                    REML estimates, are included. 
                    (See also \code{\link{infoCriteria.asreml}}.)}
 \item{\dots}{further arguments passed to \code{asreml}, \code{wald.asreml} and
                   \code{\link{as.asrtests}} via \code{\link{testranfix.asrtests}}.}
}
\value{A list containing:
       \enumerate{
         \item \code{asrtests.obj}: an \code{\link{asrtests.object}} containing the 
               components (i) \code{asreml.obj}, (ii) \code{wald.tab}, and 
               (iii) \code{test.summary}.;
         \item \code{sig.tests}: a \code{character vector} whose elements are the 
                significant terms amongst those tested.}
}
\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}
\seealso{\code{\link{chooseModel}}, \code{\link{chooseModel.data.frame}}, 
         \code{\link{as.asrtests}}, \code{\link{testranfix.asrtests}}, \cr
         \code{\link{testresidual.asrtests}}, \code{\link{REMLRT.asreml}}, 
         \code{\link{rmboundary.asrtests}}, \code{\link{newfit.asreml}}, \cr
         \code{\link{changeModelOnIC.asrtests}}, \code{\link{changeTerms.asrtests}}, 
         \code{\link{reparamSigDevn.asrtests}}}
\examples{\dontrun{
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(log.Turbidity ~ Benches + (Sources * (Type + Species)) * Date, 
                      random = ~Benches:MainPlots:SubPlots:spl(xDay), 
                      data = WaterRunoff.dat, keep.order = TRUE)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
terms.treat <- c("Sources", "Type", "Species", 
                 "Sources:Type", "Sources:Species")
terms <- sapply(terms.treat, 
                FUN=function(term){paste("Date:",term,sep="")}, 
                simplify=TRUE)
terms <- c("Date", terms)
terms <- unname(terms)
marginality <-  matrix(c(1,0,0,0,0,0, 1,1,0,0,0,0,  1,0,1,0,0,0, 
                         1,0,1,1,0,0, 1,1,1,0,1,0, 1,1,1,1,1,1), nrow=6)
rownames(marginality) <- terms
colnames(marginality) <- terms
choose <- chooseModel(current.asrt, marginality) 
current.asrt <- choose$asrtests.obj
sig.terms <- choose$sig.terms
}}
\keyword{asreml}
\keyword{htest}