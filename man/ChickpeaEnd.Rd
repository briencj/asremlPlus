\name{ChickpeaEnd.dat}
\docType{data}
\alias{ChickpeaEnd.dat}
\title{A large data set comprising the end of imaging data from a chick pea experiment conducted in high-throughput greenhouses}
\description{
The data collected after imaging  had been completed on the 1056 plants in the experiment reported by Atieno et al. (2017). The design employed for the experiment was a split-plot design in which two consecutive carts formed a main plot. The split-plot design assigned 245 genotypes to main plots, the genotypes being unequally replicated 2 or 3 times. Treatments (non-saline, saline) were randomized to the two subplots (carts) within each main plot. 

The columns in the \code{data.frame} are: Smarthouse, Lane, Position, Zone, Mainplot, Subplot, Replicate, xLane, xPosition, Genotypes, Treatments, Biomass, PlantHeight, SenescenceRank, TotalPods, FilledPods, EmptyPods, SeedNo, TotalSeedWt, SeedWt100.   

The columns Smarthouse, Lane and Position uniquely identify the rows of observations. Zones are groups of 4 Lanes, Mainplots are the 44 pairs of consecutive Subplots within each Zone, and a Subplot is a cart containing a single plant. The columns xLane and xPosition are numeric covariates for location within a Smarthouse. Genotypes and Treatments indicate the genotype and treatment that each plant was allocated.  The response variables are Biomass, PlantHeight, SenescenceRank, TotalPods, FilledPods, EmptyPods, SeedNo, TotalSeedWt and SeedWt100.
}

\usage{data(ChickpeaEnd.dat)}
\references{Atieno, J., Li, Y., Langridge, P., Dowling, K., Brien, C., Berger, B., Varshney, R. K., and Sutton, T. (2017). Exploring genetic variation for salinity tolerance in chickpea using image-based phenotyping. \emph{Scientific Reports}, \bold{7}, 1300. \doi{10.1038/s41598-017-01211-7}
}
\format{A \code{data.frames} with 1056 rows by 20 columns.}
\keyword{datasets}