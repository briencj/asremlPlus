\name{convEffectNames2DataFrame.asreml}
\alias{convEffectNames2DataFrame.asreml}
\alias{convEffectNames2DataFrame}
\title{Converts the effects names for a term stored in the component of an \code{asreml} object into a \code{\link{data.frame}}.}
\description{Converts the effects names for a term stored in the component of an \code{asreml} object into 
             a \code{\link{data.frame}} that has a column for each factor and variable in the term. 
             It facilitates adding the effects to the \code{\link{data.frame}} supplied to \code{asreml} 
             for an analysis. This function can only be used with \code{asreml} v4.2 or later.}
\usage{
\method{convEffectNames2DataFrame}{asreml}(asreml.obj, term, use = "design.matrix", sep = ":", ...)}
\arguments{
 \item{asreml.obj}{An object resulting from the fitting of 
               a model using \code{asreml} v4.2.}
 \item{term}{A \code{\link{character}} nominating a single model term. It should be the 
             name of the \code{term} in the component, of  the \code{asreml.obj}, nominated 
             in the \code{use} argument. The variables/factors in \code{term} must match those 
             in this component. Thus, it will not include \code{asreml} variance or 
             correlation functions, but may include functions such at \code{at}. Generally, 
             each element of the \code{term}, separated from other elements by colons, should 
             involve  a single factor or variable.}
 \item{use}{A \code{\link{character}} specifying a component in the \code{asreml.obj}. 
            The default is the \code{design.matrix}. Other possibilities are the 
            \code{random.coeffs}  or \code{G.aom}. For the option \code{design.matrix}, 
            the \code{design} argument of \code{asreml.options} must be set to \code{TRUE}, 
            before executing either the \code{asreml} or \code{update.asreml} call 
            produces the \code{asreml.obj}. The option \code{G.aom} specifies the \code{G} 
            component of the \code{aom} component that is included in the \code{asreml.obj} 
            if the \code{aom} argument is set to \code{TRUE}  in either \code{asreml.options} 
            before calling \code{asreml} or \code{update.asreml}, or in a call to 
            \code{asreml} or \code{update.asreml} that is involved in producing the 
            \code{asreml.obj}.}
 \item{sep}{A \code{\link{character}} specifying the separator of factors/variables 
            in the \code{term}.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A \code{\link{data.frame}} with columns for the factors and variables in \code{term}. 
            It includes the attribute \code{effect.names} that contains the extracted effects names
            for the \code{term}}


\references{Butler, D. G., Cullis, B. R., Gilmour, A. R., Gogel, B. J. and 
            Thompson, R. (2023). \emph{{ASReml-R} Reference Manual Version 4.2.}  
            VSN International Ltd, \url{https://asreml.kb.vsni.co.uk/}.} 

\author{Chris Brien}
\examples{\dontrun{
    G.dat <- convEffectNames2DataFrame(m1.asr, term = "Row:Column", use = "G.aom")
    
    G.dat <- lapply(c("at(Smarthouse, 'SW'):Lane:Position", 
                      "at(Smarthouse, 'SE'):Lane:Position"), 
                   function(term, asreml.obj)
                     tmp <- convEffectNames2DataFrame.asreml(asreml.obj, term = term),
                   asreml.obj = m1.asr)
    G.dat <- do.call(rbind, G.dat)
}}
\keyword{asreml}
