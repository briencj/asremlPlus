\name{powerTransform}
\alias{powerTransform}
\title{Performs a combination of a linear and a power transformation on a variable. The 
       transformed variable is stored in the \code{data.frame data}.}
\description{Perform a combination of a linear and a power transformation on a variable 
       whose name is given as a character string in \code{var.name}. The transformed 
       variable is stored in the \code{data.frame data}. The name of the transformed 
       variable is made by prepending to the original \code{var.name} a combination of 
       (i) \code{.offset}, if \code{offset} is nonzero, 
       (ii) \code{neg.}, if \code{scale} is -1, or \code{scaled.}, if \code{abs(scale)} 
       is other than one, and (iii) either \code{log.}, \code{sqrt.}, \code{recip.} or 
       \code{power.}, if \code{power} is other than one. 
       No action is taken if there is no transformation (i.e. \code{offset = 0}, 
       \code{scale = 1} and \code{power = 1}). Also, the \code{titles} list is 
       extended to include a component with a generated title for the 
       transformed variable with text indicating the transformation prepended 
       to the \code{title} for the \code{var.name} obtained from the 
       \code{titles list}.  For nonzero \code{offset}, `Offset ' is prepended, 
       For \code{scaled} not equal to one, the possible prepends are `Negative of '  
       and `Scaled '. The possible prepended texts for \code{power} not equal to one 
       are `Logarithm of', `Square root of ', `Reciprocal of ' and `Power nnnn of ', 
       where \code{nnn} is the power used.}
\usage{powerTransform(var.name, power = 1, offset = 0, scale = 1, titles = NULL, data)}
\arguments{
 \item{var.name}{A \code{character string} specifying the name of the variable 
                 in the \code{data.frame data} that is to be transformed.}
 \item{power}{A number specifying the power to be used in the transformation. 
              If equal to 1, the default, no power transformation is applied. 
              Otherwise, the variable is raised to the specified \code{power}, 
              after scaling and applying any nonzero \code{offset}.
              If \code{power = 0}, the natural logarithm is used to transform 
              the response; however, if the smallest value to be log-transformed 
              is less than 1e-04, an error is generated. A log-transformation in 
              this situation may be possible if a nonzero \code{offset} and/or a 
              \code{scale} not equal to one is used.}
 \item{offset}{A number to be added to each value of the variable, after any scaling 
               and before applying any power transformation.}
 \item{scale}{A number to multiply each value of the variable, before adding 
               any offset and applying any power transformation.}
 \item{titles}{A \code{character vector}, each element of which is named for a 
               variable in \code{data} and is a \code{character string} giving a title 
               to use in output (e.g. tables and graphs) involving the variable. If
               titles are not supplied, the column name of the variable in \code{data} 
               is used.}
 \item{data}{A \code{data.frame} containing the variable to be transformed and 
             to which the transformed variable is to be appended.}
}
\value{A \code{list} with a component named \code{data} that is the \code{data.frame} containing the 
       transformed variable, a component named tvar.name that is a \code{character string} 
       that is the name of the transformed variable in \code{data}, and a component named 
       \code{titles} that extends the \code{list} supplied in the \code{titles} argument to 
       include a generated title for the transformed title, the name of the new component being 
       \code{tvar.name}.}
\author{Chris Brien}
\seealso{\code{\link{angular}}, \code{\link{angular.mod}}.}
\examples{
## set up a factor with labels
x.dat <- data.frame(y = c(14, 42, 120, 150))
 
## transform y to logarithms
trans <- powerTransform("y", power = 0, titles=list(y = "Length (cm)"), data = x.dat)
x.dat <- trans$data
tvar.name <- trans$tvar.name

## transform y to logarithms after multiplying by -1 and adding 1.
z.dat <- data.frame( y = c(-5.25, -4.29, -1.22, 0.05))
trans <- powerTransform("y", power = 0, scale = -1, offset = 1 ,
                        titles=list(y = "Potential"), data = z.dat)
z.dat <- trans$data
tvar.name <- trans$tvar.name
}
\keyword{manip}