\name{newfit.asreml}
\alias{newfit.asreml}
\alias{newfit}
\title{Refits an \code{asreml} model with changed arguments by extracting, modifying 
             and evaluating its \code{call} - an alternate to update.asreml.}
\description{Extracts the call from the \code{asreml.obj} and evaluates that call, 
             replacing supplied \code{asreml} arguments with changed values. 
             If \code{update} is \code{TRUE} and \code{set.terms} is \code{NULL}, 
             the call is evaluated using the variance parameter estimates 
             from the \code{call} stored in \code{asreml.obj}; if some 
             variance terms in the newly fitted model are singular (S) or 
             bound (B), a refit of the model will be tried in which the 
             variance parameter estimates in \code{asreml.obj} are not used 
             and will become the fitted model if its boundary terms are a 
             subset of those in the fitted model stored in \code{asreml.obj}. 
             If \code{update} is \code{FALSE} or \code{set.terms} is not 
             \code{NULL}, the previous values of the variance parameters 
              are not used as initial values for fitting the new model; 
             \code{G.param} and \code{R.param} are set to \code{NULL} or 
             to values as specified for \code{set.terms}. The \code{\dots} 
             argument can be used to pass \code{G.param} and/or 
             \code{R.param}, provided \code{update} is \code{FALSE} and 
             \code{set.terms} is \code{NULL}. 

             Irrespective of whether \code{set.terms} is \code{NULL}, 
             variance parameter names, bounds and initial.values stored in 
             the \code{setvparameters} \code{\link{data.frame}} located in the 
             \code{call} component of the \code{asreml.obj} are added to any 
             \code{set.terms} supplied in the current call. except for those 
             in \code{setvparameters} that are specified in the 
             \code{set.terms} argument. In the process of fitting the model, 
             the \code{setvparameters} \code{\link{data.frame}} stored in the 
             supplied \code{asreml.obj} is updated to include the prior and 
             current \code{set.terms}. Then, just before exiting 
             \code{newfit.asreml}, a search for changes in the \code{bound}s 
             of the stored terms is made. If any have changed, an attempt is 
             made to force the values back to their values on entry.}
\usage{
\method{newfit}{asreml}(asreml.obj, fixed., random., sparse., 
       residual., rcov., update = TRUE, trace = FALSE, 
       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
       keep.order = TRUE, 
       set.terms = NULL, ignore.suffices = TRUE, 
       bounds = "P", initial.values = NA, ...)}
\arguments{
 \item{asreml.obj}{A valid \code{asreml} object with 
                   a component named \code{call} (from a previous 
                   call to either \code{asreml} or \code{update.asreml}).}
 \item{fixed.}{A character or formula specifying changes to the fixed 
               formula. This is a two-sided formula 
               where "." is substituted for existing components in the 
               fixed component of \code{asreml.obj$call}. If changes are 
               specified, the fixed terms will be reordered so that 
               single-variable terms come first, followed by two-variable 
               terms and so on.}
 \item{random.}{A character or formula specifying changes to the random 
                formula. This is a one-sided formula 
                where "." is substituted for existing components in the 
                random component of \code{asreml.obj$call}.}
 \item{sparse.}{A character or formula specifying changes to the sparse 
                formula. This is a one-sided formula 
                where "." is substituted for existing components in the 
                sparse component of \code{asreml.obj$call}.}
 \item{residual.}{A character or formula specifying changes to the error 
                formula, used when version 4 or later of ASReml-R is loaded.
                This is a one-sided formula 
                where "." is substituted for existing components in the 
                residual component of \code{asreml.obj$call}.}
 \item{rcov.}{A character or formula specifying changes to the error 
                formula, used when version 3 of ASReml-R is loaded. 
                This is a one-sided formula 
                where "." is substituted for existing components in the 
                residual component of \code{asreml.obj$call}.}
 \item{update}{A logical indicating whether to use the variance parameter 
               estimates in fitting a new model using \code{\link{newfit.asreml}}. 
               If \code{update} is \code{TRUE} and \code{set.terms} is \code{NULL}, 
               the call stored in the \code{asreml.obj} is evaluated using the 
               variance parameter estimates stored in \code{R.param} and 
               \code{G.param}, as well as the changes specified in the 
               arguments to \code{newfit.asreml}. If \code{FALSE} or   
               \code{set.terms} is not \code{NULL}, then the use of previous 
               values of the variance parameters as initial values is not 
               enforced; \code{G.param} and \code{R.param} are set to 
               \code{NULL} or to values as specified for \code{set.terms}.}
 \item{trace}{A \code{\link{logical}} that control output from ASReml-R. 
              If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
              even when it does not converge. If \code{FALSE} and the fit does not 
              converge, the supplied \code{asreml.obj} is returned.}
 \item{allow.fixedcorrelation}{A \code{logical} indicating whether to accept a new model 
              even when it contains correlations in the model whose values have been 
              designated as fixed, bound or singular. If \code{FALSE} and the new model 
              contains correlations whose values have not been able to be estimated, 
              the supplied \code{asreml.obj} is returned. The fit in the supplied 
              the \code{asreml.obj} will also be tested and a warning issued if 
              both fixed correlations are found in it and 
              \code{allow.fixedcorrelation} is \code{FALSE}.}
 \item{keep.order}{A logical value indicating whether the terms should keep 
               their positions. If \code{FALSE} the terms are reordered so 
               that main effects come first, followed by the interactions, 
               all second-order, all third-order and so on. Effects of a 
               given order are kept in the order specified.}
 \item{set.terms}{A character vector specifying the terms that are to have 
               bounds and/or initial values set prior to fitting.
               The names must match those in the \code{vparameters} component of the 
               new \code{asreml.obj}.}
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
 \item{\dots}{additional arguments to the call, or arguments with changed 
            values.}
}
\value{An \code{asreml} object.}


\references{Butler, D. G., Cullis, B. R., Gilmour, A. R., Gogel, B. J. and 
            Thompson, R. (2023). \emph{{ASReml-R} Reference Manual Version 4.2.}  
            VSN International Ltd, \url{https://asreml.kb.vsni.co.uk/}.}

\author{Chris Brien}
\seealso{\code{\link{convAsremlobj.asreml}}, \code{update.asreml}, \code{\link{setvarianceterms.call}}}
\examples{\dontrun{
    m2.asreml <- newfit(m1.asreml, random. = "~ . - Blocks:Plots", maxiter=75)
}}
\keyword{asreml}