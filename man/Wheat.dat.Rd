\name{Wheat.dat}
\docType{data}
\alias{Wheat.dat}
\title{Data for an experiment to investigate 25 varieties of wheat}
\description{
The data appears in Gilmour et al. [1995] and is from a field experiment designed to compare the
performance of 25 varieties of wheat. An analysis of it using asreml is presented by 
Butler et al. (2010, Sectoion 8.6), although they suggest that it is a barley experiment. It is used in \code{\link{asremlPlus-package}} as an executable example of the use of the \code{asremlPlus} to analyse a data set.

The experiment was conducted at Slate Hall Farm, UK, in
1976 and was designed as a balanced lattice square with 6 replicates laid out in a 10 \eqn{\times}{x} 15 rectangular
grid. The columns in the data frame are: Rep, Row, Column, WithinColPairs, Variety, yield.
The response variable is the grain yield.

}
\usage{data(Wheat.dat)}
\format{A data.frame containing 150 observations of 6 variables.}
\source{
Butler, D. G., et al. (2010). \emph{Analysis of Mixed Models 
            for {S} language environments: {ASReml-R} reference manual.}  
            Brisbane, DPI Publications.
            
            Gilmour, A. R., et al. (1995) Average Information REML: An efficient algorithm for 
            variance parameter estimation in linear mixed models. \emph{Biometrics}, 
            \bold{51}, 1440-1450.
}
\author{Chris Brien}
\keyword{datasets}