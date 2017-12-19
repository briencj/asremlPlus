\name{rmboundary.asrtests}
\alias{rmboundary.asrtests}
\alias{rmboundary}
\title{Removes any boundary or singular variance components from the fit stored in 
       \code{asrtests.obj} and records their removal in a data.frame.}
\description{Any terms specified in the random model that are 
  estimated on the boundary or are singular and can be removed are removed from the fit stored in the 
  \code{asreml} object stored in the \code{asrtests} object. Terms that specify multiple parameters in the random model 
  cannot be removed (e.g. terms specified using the \code{at} function with more 
  than one level of the factor) and terms in \code{residual} model are not removed. 
  Terms that can be removed are selected for removal in the following order based on whether they involve:  
  (i) a \code{dev} function, (ii) only factors, (iii) an \code{spl} function, 
  (iv) a \code{pol} function and (v) a \code{lin} function or a variable that
  is an \code{integer} or a \code{numeric}. It should be noted that this order of removal 
  presumes that random deviation terms are specified via the \code{dev} function 
  rather than via a random factor. Once the earliest of the above classes 
  with a boundary term is identified, a term within this class is selected for 
  removal. For all classes, except for factor-only terms, the smallest term 
  with the largest number of  variables/factors is removed. 
  Amongst factor-only terms, the smallest term with the smallest number of 
  variables/factors is removed. After each 
  variance component is removed, a row for it is added to the 
  \code{test.summary} \code{data.frame} and the model refitted. If there are further 
  boundary or singular terms, one is removed using the above strategy. 
  This process continues until there are no further boundary or singular 
  variance components that are removable. Other types of boundary or singular terms, 
  which cannot be removed, are reported in warning messages.}
\usage{
\method{rmboundary}{asrtests}(asrtests.obj, checkboundaryonly = FALSE, 
           trace = FALSE, update = TRUE,
           set.terms = NULL, ignore.suffices = TRUE, 
           bounds = "P", initial.values = NA, ...)}
\arguments{
 \item{asrtests.obj}{an \code{\link{asrtests}} object for a fitted model that is a list 
                containing an \code{asreml} object, a \code{wald.tab} 
                \code{data.frame} with 4 columns, and a \code{data.frame} with 
                5 columns that records any previous changes and tests in the 
                fitted model.}
 \item{checkboundaryonly}{if \code{TRUE} then boundary and singular terms are not removed by 
              \code{\link{rmboundary.asrtests}}; a warning is issued instead.}
 \item{trace}{if \code{TRUE} then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{update}{if \code{TRUE} then \code{update.asreml} is called to fit the model with 
              any boundary terms removed. In doing this the arguments 
              \code{R.param} and \code{G.param} are set to those in the \code{asreml} 
              object stored in \code{asrtests.obj} so that the values from the previous 
              model are used as starting values. If \code{FALSE} then a call is made to 
              \code{asreml} in which the only changes from the previous call are that 
              (i) the terms  for boundary variance components are removed 
              from the models and (ii) modifications specified via \code{\dots} are made.}
 \item{set.terms}{a character vector specifying the terms that are to have 
                  bounds and/or initial values set prior to fitting.}
 \item{ignore.suffices}{a logical vector specifying whether the suffices of the 
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
 \item{bounds}{a \code{\link{character}} vector specifying the bounds to be applied 
                    to the terms specified in \code{set.terms}. This vector 
                    must be of length one or the same length as \code{set.terms}. 
                    If it is of length one then the same constraint is 
                    applied to all the terms in \code{set.terms}. 
                    If any of the bounds are equal to NA then they are 
                    left unchanged for those terms.}
 \item{initial.values}{a character vector specifying the initial values for  
                    the terms specified in \code{terms}. This vector 
                    must be of length one or the same length as \code{terms}. 
                    If it is of length one then the same initial value is 
                    applied to all the terms in \code{terms}. 
                    If any of the initial.values are equal to NA then they are 
                    left unchanged for those terms.}
 \item{\dots}{further arguments passed to \code{asreml}.}
}
\value{An \code{\link{asrtests}} object, which is a list containing:
       \enumerate{
         \item \code{asreml.obj}: an \code{asreml} object containing the fit 
               of the model after all boundary and singular terms have been removed;
         \item \code{wald.tab}: a 4-column \code{data.frame} containing a 
               pseudo-anova table for the fixed terms produced by \code{wald.asreml};
         \item \code{test.summary}: a \code{data.frame} with columns \code{term}, 
                \code{DF}, \code{denDF}, \code{p} and \code{action}. A row is added to 
                it for each of the boundary terms removed, the row containing the 
                name of the term, one for the DF, NA for the p-value and \code{Boundary} for 
                the action.}
}

\author{Chris Brien}
\seealso{\code{\link{asrtests}}, \code{\link{changeTerms.asrtests}}, 
         \code{\link{testranfix.asrtests}}, \code{\link{testresidual.asrtests}},  \cr
         \code{\link{newfit.asreml}}, \code{\link{reparamSigDevn.asrtests}}, 
         \code{\link{chooseModel.asrtests}}}
\examples{\dontrun{
current.asrt <- rmboundary(current.asrt)
}}
\keyword{asreml}
\keyword{htest}