\name{R2adj.asreml}
\alias{R2adj.asreml}
\alias{R2adj}
\title{Calculates the adjusted coefficient of determination for a specified combination of fixed and random terms.}
\description{Calculates the adjusted coefficient of determination (R2) that measures the contributions to the total 
             variance exhibited by the observations of a specified combination of fixed and random terms in a 
             fitted linear mixed model. This function is not available in ASReml-R version 3.}
\usage{
\method{R2adj}{asreml}(asreml.obj, 
      include.which.fixed = ~ ., orthogonalize = "hybrid", 
      include.which.random = NULL, 
      bound.exclusions = c("F","B","S","C"), ...)}
\arguments{
 \item{asreml.obj}{An \code{asreml} object returned from a call to \code{asreml}.}
 \item{include.which.fixed}{A \code{\link{formula}} specifying the fixed terms whose joint contribution to 
               the total variance is to be measured.  If it is \code{NULL}, no fixed term 
               is to be included in the terms whose contributions is to be assessed. 
               The formula \code{~ .} indicates that the joint contribution of all fixed terms are to 
               be measured. Otherwise, the joint contribution of the set of terms in the formula will 
               be assessed. The names of the terms in the \code{\link{formula}} must be the same as 
               those in either the Wald table produced by \code{wald.asreml} or the \code{terms} attribute 
               of the \code{fixed} component of the \code{coefficient} component of the supplied 
               \code{asreml.obj}. 
               
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
               to be measured. Otherwise, the joint contribution of the set of terms in the formula 
               will be assessed. The terms in the \code{\link{formula}} must be on of those occurring 
               in the \code{vparameters} component of the supplied \code{asreml.obj}.}
 \item{bound.exclusions}{A \code{character} specifying one or more bound codes that 
                    will result in a variance parameter in the \code{random} model being excluded 
                    from contributing to the variance. If set to \code{NULL} then none will
                    be excluded.}
 \item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
}
 \value{A \code{numeric} that is the adjusted R2, expressed as a percentage.}

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
             
        If the error \code{Matrix is not idempotent} should occur then, especially if 
        there are  many terms, one might try using \code{set.daeTolerance} from the 
        \code{dae} package to reduce the tolerance used in determining if values are 
        either the same or are zero; it may be necessary to lower the tolerance to as 
        low as 0.001. Also, setting \code{orthogonalize} to \code{eigenmethods} is 
        worth a try.

        Note that the adjusted R2 can be negative, which indicates that the contribution 
        of the terms to the total variance is very small relative to the sum of the 
        degrees of freedom of the terms.
          
        The function \code{estimateV.asreml} is used to calculate the variance matrices 
        required in calculating the adjusted R2.
          
        Piepho's (2023) method for GLMMs has not been implemented. 
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
  R2.adj.add <- R2adj.asreml(current.asr, include.which.fixed = ~ Nitrogen + Variety)
}}
\keyword{asreml}
