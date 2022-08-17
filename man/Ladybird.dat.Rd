\name{Ladybird.dat}
\docType{data}
\alias{Ladybird.dat}
\title{Data for an experiment to investigate whether ladybirds transfer aphids}
\description{
Welham et al. (2015, Example 8.2) describe a three-factor factorial experiment to investigate whether ladybirds transfer fungus to live aphids on plants. The three factors are Host plant (beans, trefoil), infected Cadavers (5, 10, 20), and Ladybird (-, +). A generalized randomized complete-block design is used to assign the three factors  to 2 Runs, each of which involves 36 containers with a plant and live aphids. The response to be analyzed is the logit of the proportion of live aphids that were infected.

The columns in the data frame are: ID, Run, Plant, Host, Ladybird, Cadavers, Live, Infected, logitP, Prop. The column ID numbers the observations. 
Live, Infected, logitP, Prop are response variables.

}
\usage{data(Ladybird.dat)}
\format{A data.frame containing 72 observations of 10 variables.}
\source{
Welham, S. J., Gezan, S. A., Clark, S. J., & Mead, A. (2015). \emph{Statistical Methods in Biology: Design and Analysis of Experiments and Regression}. Boca Raton: Chapman and Hall/CRC.. 
}
\author{Chris Brien}
\keyword{datasets}