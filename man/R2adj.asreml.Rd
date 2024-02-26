\name{R2adj.asreml}
\alias{R2adj.asreml}
\alias{R2adj}
\title{Calculates the adjusted coefficient of determination for a specified combination of fixed and random terms.}
\description{Calculates the adjusted coefficient of determination (R2) that measures the contributions to 
             the total variance exhibited by the observations of a specified combination of fixed and random 
             terms in a fitted linear mixed model.  
             
             Note that the adjusted R2 can be negative, which indicates that the contribution of the terms 
             to the total variance is very small relative to the sum of the degrees of freedom of the terms.
          
             Piepho's (2023) method for GLMMs has not been implemented. This function is not available for 
             ASReml-R version 3.
}
\usage{
\method{R2adj}{asreml}(asreml.obj, 
      include.which.fixed = ~ ., orthogonalize = "hybrid", 
      include.which.random = NULL, 
      bound.exclusions = c("F","B","S","C"), ...)}
\arguments{
 \item{asreml.obj}{An \code{asreml} object returned from a call to \code{asreml}.}
 \item{include.which.fixed}{A \code{\link{formula}} specifying the fixed terms whose joint contribution 
               to the total variance is to be measured.  If it is \code{NULL}, no fixed term 
               is to be included in the terms whose joint contribution is to be assessed. The 
               \code{\link{formula}} \code{~ .} indicates that the joint contribution of all fixed 
               terms are to be measured. Otherwise, the joint contribution of the set of terms 
               specified by the \code{\link{formula}} will be assessed. The \code{\link{formula}} can 
               include a "\code{.}", which means all fixed terms currently fitted, and is most likely 
               followed by a "\code{-}" with a bracketed set of terms to be removed that can be 
               specified using \code{\link{formula}} operators.  The names of the resulting 
               terms must be the same as those in either the \code{terms} attribute of the 
               \code{fixed} component of the \code{coefficient} component of the supplied 
               \code{asreml.obj}, or the Wald table produced by \code{wald.asreml}. 
               
               Note that the contribution of a subset of the fixed terms is only unique if the effects 
               for the fixed terms are orthogonal; if the effects are not orthogonal then the 
               contributions will depend on the order of the terms in the formula. Also, determining 
               the joint contribution of a subset of the fixed terms in the model may be 
               computationally demanding because the projection matrices have to be formed for all 
               fixed terms and these projections matrices have to be orthogonalized. A heavy 
               computational burden is most likely when the effects for the fixed terms are not 
               orthogonal, for example, when numeric covariates are included amongst the terms.}
 \item{orthogonalize}{A \code{\link{character}} vector indicating the method for orthogonalizing a 
                      projector to those for terms that occurred previously in the 
                      \code{\link{formula}} for \code{include.which.fixed}. Orthogonalizing the 
                      projectors of fixed terms is not performed for the default setting of 
                      \code{. ~}.  WHen required, two options are available for orthogonalizing: 
                      \code{hybrid} and \code{eigenmethods}. The \code{hybrid} option is the most 
                      general and uses the relationships between the projection operators for the 
                      terms in the \code{formula} to decide which \code{\link{projector}}s to 
                      subtract and which to orthogonalize using eigenmethods. The \code{eigenmethods} 
                      option recursively orthogonalizes the \code{\link{projector}}s 
                      using an eigenanalysis of each \code{\link{projector}} 
                      with previously orthogonalized \code{\link{projector}}s. 
                      See the documentation for \code{porthogonalize.list} from the R package 
                      \code{dae} for more information.}
 \item{include.which.random}{A \code{\link{formula}} specifying the random terms whose joint 
               contribution  to the total variance is to be measured.  If it is \code{NULL}, no 
               random term is to be included in the terms whose joint contribution is to be assessed. 
               The formula \code{~ .} indicates that the joint contribution of all random terms is 
               to be measured. Otherwise, the joint contribution of the set of terms specified by 
               the \code{\link{formula}}  will be assessed. The \code{\link{formula}} can 
               include a "\code{.}", which means all random terms currently fitted, and is most likely 
               followed by a "\code{-}" with a bracketed set of terms to be removed that can be 
               specified using \code{\link{formula}} operators. The resulting terms must be one of 
               those occurring in either the \code{vparameters} component of the supplied 
               \code{asreml.obj}, or in the \code{terms} attribute of the \code{random} component of 
               the \code{coefficient} component of the supplied \code{asreml.obj}.}
 \item{bound.exclusions}{A \code{character} specifying one or more bound codes that 
                    will result in a variance parameter in the \code{random} model being excluded 
                    from contributing to the variance. If set to \code{NULL} then none will
                    be excluded.}
 \item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
}
 \value{A \code{numeric} that is the adjusted R2, expressed as a percentage. It has attributes 
        \code{include.which.fixed}, \code{include.which.random} and \code{missing.termmatrix} 
        (use \code{attr(x, which = "name")} to access the attribute \code{name}). 
        The \code{missing.termmatrix} attribute will be \code{NULL}, unless the design matrix 
        could not be obtained for one or more model terms. If is is not \code{NULL}, it will be 
        a list of terms whose design matices could not be produced  and so are not included in 
        the variance matrix estimate. An \code{NA} will be returned for the adjusted R2 if 
        \code{missing.termmatrix} is not \code{NULL} or a generalized inverse could 
        not be computed for the variance matrix estimate.}

 \details{The method used to compute the adjusted R2 under a linear mixes model (LMM) is that 
          described by Piepho (2023). Here, the method has been extended to allow computation  
          of the adjusted R2 for a subset of the fixed terms. A set of orthogonalized 
          projectors for all of the fixed terms in the model (a set of 
          \eqn{\mathbf{Q}_i\mathrm{s}}{Q_i}) is obtained and the combined contribution of the 
          fixed terms nominated in \code{include.which.fixed} is obtained by computing the 
          average semisquared bias, ASSB, for the nominated fixed terms as:
\deqn{\Sigma_i \{(\mathbf{Q}_i \mathbf{X}\boldsymbol{\beta})^\mathrm{T}\mathbf{Q}_i \mathbf{X}\boldsymbol{\beta} + \textnormal{trace}(\mathbf{X}^\mathrm{T} \mathbf{Q}_i \mathbf{X} \mathrm{var}(\boldsymbol{\beta})) \} / (n - 1)}{sum((Q_i X\beta)^T Q_1 X\beta + trace(X^T Q_i X var(\beta))) / (n-1).}

        Of the two methods, \code{eigenmethods} is least likely to fail, but it does 
        not establish the marginality between the terms. It is often needed when 
        there is nonorthogonality between terms, such as when there are several linear 
        covariates. It can also be more efficient in these circumstances.
        
        The process can be computationally expensive, particularly for a large data set 
        (500 or more observations) and/or when many terms are to be orthogonalized, 
        particularly if they are not orthogonal. 
             
        If the error "\code{Matrix is not idempotent}" should occur then, especially if 
        there are  many terms, one might try using \code{set.daeTolerance} from the 
        \code{dae} package to reduce the tolerance used in determining if values are 
        either the same or are zero; it may be necessary to lower the tolerance to as 
        low as 0.001. Also, setting \code{orthogonalize} to \code{eigenmethods} is 
        worth a try.

        In doing the computations, no changes are made to the fitted model, nor is the 
        \code{\link{formula}} stored in \code{asreml.obj} referred to. Instead, the 
        names of the terms referred to are those stored in the \code{coefficients} 
        component of the \code{asreml.obj}.  Use 
        \code{attr(asreml.obj$coefficients$fixed, which = "terms")} to access the attribute 
        for fixed terms; substitute \code{random} for \code{fixed} to see the names of the 
        random terms. For fixed terms. the term names are the same as those in the 
        Wald table produced by \code{wald.asreml}, and, for random terms, the same as those 
        in the \code{vparameters} component of the \code{asreml.obj}. Two \code{asreml} 
        \code{\link{formula}} functions whose terms can differ from their formulation in a 
        model \code{\link{formula}} are \code{at} and \code{str}.)
             
        The function \code{estimateV.asreml} is used to calculate the variance matrices 
        required in calculating the adjusted R2.
}

\references{Piepho, H.-P. (2023). An adjusted coefficient of determination (R2) for generalized 
            linear mixed models in one go. \emph{Biometrical Journal}, \bold{65(7)}, 2200290. 
            \doi{10.1002/bimj.202200290}.}

\author{Chris Brien}
\seealso{\code{asreml}, \code{\link{estimateV.asreml}}.}
\examples{\dontrun{
  data(Oats.dat)
  
  current.asr <- asreml(Yield ~ Nitrogen*Variety, 
                        random=~Blocks/Wplots,
                        data=Oats.dat)
  R2.adj.fix <- R2adj.asreml(current.asr)
  R2.adj.ran <- R2adj.asreml(current.asr, 
                             include.which.fixed = NULL, include.which.random = ~ .)
  R2.adj.tot <- R2adj.asreml(current.asr, include.which.random = ~ .)
  R2.adj.tot <- R2adj.asreml(current.asr, include.which.random = ~ Blocks)
  R2.adj.add <- R2adj.asreml(current.asr, include.which.fixed = ~ Nitrogen + Variety)
  R2.adj.int <- R2adj.asreml(current.asr, 
                             include.which.fixed = ~ . - (Nitrogen + Variety))
  R2.adj.int <- R2adj.asreml(current.asr, include.which.fixed = ~ Nitrogen:Variety)
}}
\keyword{asreml}
