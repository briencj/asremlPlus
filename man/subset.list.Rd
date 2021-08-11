\name{subset.list}
\alias{subset.list}
\title{Forms a \code{\link{list}} that contains a subset of the components of the supplied 
       \code{\link{list}}.}
\description{Select components of a \code{\link{list}} specified by a list of numbers or names,
             or by a logical indicating for each component of the list whether or not it is to 
             be retained.}
\usage{\method{subset}{list}(x, select = 1:length(x), ...)}
\arguments{
 \item{x}{An \code{\link{list}} object.}
 \item{select}{A \code{\link{numeric}} or \code{\link{character}} that lists or names the 
               components of the \code{\link{list}} \code{x} that are to be retained in the 
               subset. It can also be a \code{\link{logical}} that is the same length as 
               \code{x} and indicates whether or not a component is to be retained.}
 \item{\dots}{further arguments passed to or from other methods. Not used at present.}
}
\value{A \code{\link{list}} with just the subset of the components from \code{x}. If the 
       components of \code{x} are named, then these names are retained in the subset 
       \code{\link{list}}.}

\author{Chris Brien}
\seealso{\code{\link{subset.alldiffs}}}
\examples{
x <- list(1:3,letters[1:3],LETTERS[1:3])
y <- subset.list(x, select = c(1,3))
y <- subset.list(x, select = c(TRUE,FALSE,TRUE))

names(x) <- LETTERS[1:3]
y <- subset.list(x, select = c(1,3))
z <- subset.list(x, select = LETTERS[c(1,3)])
x <- list(1:3,letters[1:3],LETTERS[1:3])
names(x)[c(1,3)] <- LETTERS[c(1,3)]
z <- subset.list(x, select = c(1,2))
v <- subset.list(x)
}
\keyword{asreml}