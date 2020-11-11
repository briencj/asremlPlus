\name{Wheat.dat}
\docType{data}
\alias{Wheat.dat}
\title{Data for a 1976 experiment to investigate 25 varieties of wheat}
\description{
The data appears in Gilmour et al. (1995) and is from a field experiment designed to compare the
performance of 25 varieties of spring wheat. An analysis of it using asreml is presented by 
Butler et al. (2018, Section 7.6), although they suggest that it is a barley experiment. 
It is used in the Wheat vignettes  [Enter \code{vignette(package = "asremlPlus")}] as an 
executable example of the use of the \code{asremlPlus} to analyse a data set.

The experiment was conducted at Slate Hall Farm, UK, in 1976 and was designed as 
a balanced lattice square with 6 replicates laid out in a 10 \eqn{\times}{x} 15 rectangular
grid. The columns in the data frame are: Rep, Row, Column, WithinColPairs, Variety, yield.
The response variable is the grain yield.

}
\usage{data(Wheat.dat)}
\format{A data.frame containing 150 observations of 6 variables.}
\source{
Butler, D. G., Cullis, B. R., Gilmour, A. R., Gogel, B. J. and 
            Thompson, R. (2018). \emph{{ASReml-R} Reference Manual Version 4.}  
            VSN International Ltd, \url{https://asreml.org/}.
            
Gilmour, A. R., et al. (1995) Average Information REML: An efficient algorithm for 
            variance parameter estimation in linear mixed models. \emph{Biometrics}, 
            \bold{51}, 1440-1450.
}
\author{Chris Brien}
\keyword{datasets}