\name{makeTPPSplineMats.data.frame}
\alias{makeTPPSplineMats.data.frame}
\alias{makeTPPSplineMats}
\title{Make the spline basis matrices and data needed to fit Tensor Product P-splines.}
\description{Prepares the fixed and random P-spline basis matrices, and associated 
       information, that are needed for fitting of Tensor Product P-splines (TPPS) 
       as described by Rodriguez-Alvarez et al. (2018). This function is called 
       internally by \code{\link{addSpatialModelOnIC.asrtests}}, 
       \code{\link{addSpatialModelOnIC.asrtests}} and 
       \code{\link{chooseSpatialModelOnIC.asrtests}} when fitting \code{TPPS} models 
       for local spatial variation. There are two methods available, controlled by 
       \code{asrem.option} for creating and storing the basis functions. 
       This function is most likely to be called directly when \code{mbf} has been used 
       in creating an \code{asreml.object} and it is desired to use the object in a 
       session subsequent to the session in which the object was created.
}
\usage{
\method{makeTPPSplineMats}{data.frame}(data, sections = NULL, 
                  row.covar, col.covar, 
                  nsegs = NULL, nestorder = c(1,1), 
                  degree = c(3,3), difforder = c(2,2),
                  rotateX = FALSE, theta = c(0,0), 
                  asreml.option = "mbf", mbf.env = sys.frame(), 
                  ...)
}
\arguments{
 \item{data}{An \code{\link{data.frame}} that holds the spline bases for a section. 
             It is indexed by columns named \code{col} and \code{row}.}
 \item{sections}{A single \code{character} string that species the name of the column 
               in the \code{\link{data.frame}} that contains the \code{\link{factor}} 
               that identifies different sections of the data 
               to which separate spatial models are to be fitted.}
 \item{row.covar}{A single \code{\link{character}} string nominating a \code{\link{numeric}} 
               column in the \code{\link{data.frame}} that contains the values of a 
               covariate indexing the rows of the grid.}
 \item{col.covar}{A single \code{\link{character}} string nominating a \code{\link{numeric}} 
               column in the \code{\link{data.frame}} that contains the values of a 
               covariate indexing the columns of the grid.}
 \item{nsegs}{A pair of \code{\link{numeric}} values giving the number of segments into
              which the column and row ranges are to be split, respectively, for fitting 
              a P-spline model (\code{TPPS}) (each value specifies the number of internal 
              knots + 1). If not specified, then (number of unique values - 1) is used 
              in each dimension; for a grid layout with equal spacing, this gives a knot 
              at each data value. If \code{sections} is not \code{NULL} and the grid 
              differs between the \code{sections}, then \code{nsegs} will differ between 
              the \code{sections}.}
 \item{nestorder}{A \code{\link{numeric}} of length 2. The order of nesting for column 
              and row dimensions, respectively, in fitting a P-spline model (\code{TPPS}). 
              A value of 1 specifies no nesting, a value of 2 generates a spline with 
              half the number of segments in that 
              dimension, etc. The number of segments in each direction must be a 
              multiple of the order of nesting.}
 \item{degree}{A \code{\link{numeric}}  of length 2. The degree of polynomial spline to 
              be used for column and row dimensions respectively, in fitting a P-spline 
              (\code{TPSS}).}
 \item{difforder}{A \code{\link{numeric}} of length 2. The order of differencing for 
              column and row dimensions, respectively, in fitting a P-spline 
              (\code{TPSS}).}
 \item{rotateX}{A \code{\link{logical}} indicating whether to rotate the eigenvectors of the 
              penalty matrix, as described by Piepho, Boer and Williams (2022), when 
              fitting a P-spline (\code{TPSS}). Setting \code{rotateX} to \code{TRUE} 
              results in an optimized rotation produced by exploring 
              a two-dimensional grid of rotation angle pairs, and for each pair 
              analyzing the data under a model that omits the random interaction terms. 
              The angle pair with the minimum deviance is used to apply an optimized 
              rotation. Rotation of the eigenvectors is only relevant for 
              \code{difforder} values greater than 1 and has only been implemented for 
              \code{difforder} equal to 2.}
 \item{theta}{A \code{\link{numeric}} of length 2. The angle (in degrees) to be used in 
              rotating the eignevectors of the penalty matrix of a P-spline (\code{TPSS}).}
 \item{asreml.option}{A single \code{character} string specifying whether the \code{grp} or 
              \code{mbf} methods are to be used to supply externally formed covariate 
              matrices to \code{asreml} when fitting a P-spline (\code{TPSS}). Compared 
              to the \code{mbf} method, the \code{grp} method is somewhat faster, but 
              creates large \code{\link{asrtests.object}}s for which the time it takes 
              to save them can exceed any gains in execution speed. The \code{grp} 
              method adds columns to the \code{\link{data.frame}} containing the \code{data}. 
              On the other hand, the \code{mbf} method adds only the fixed covariates 
              to \code{data} and stores the random covariates in the environment of the 
              internal function that calls the spline-fitting function; there are three smaller 
              \code{\link{data.frame}s} for each \code{section} that are not stored in the 
              \code{asreml.object} resulting from the fitted model.}
 \item{mbf.env}{A \code{\link{environment}} specifying the environment to which the 
              \code{\link{data.frame}s} containing the spline bases are to be assigned. 
              If \code{mbf.env} is \code{NULL}, the data.frames will not be assigned.}
\item{\dots}{Further arguments passed to \code{tpsmmb} from package \code{TPSbits}.}
}
\details{
The objects are formed using the function \code{tpsmmb} from the \code{R} package \code{TPSbits} 
authored by Sue Welham (2022). This function has been extended to allow for \code{sections} (see below) 
and to allow rotation of the penalty matrix for the linear component of the interaction terms in 
\code{TPPCS} models (for more information about rotation see Piepho, Boer and Williams, 2022).

Each combination of a \code{row.covar} and a \code{col.covar} does not have to specify a 
single observation; for example, to fit a local spatial variation model to the main units 
of a split-unit design, each combination would correspond to a main unit and all subunits 
of the main unit would would have the same combination.

The data for experiment can be divided \code{sections} and the spline bases and associated 
data will be produced for each section. If there is more than one \code{sections}, then a 
\code{\link{list}} is returned that has a component for each section. The component for 
a section contains:
}

\value{A \code{\link{list}} of length equal to the number of sections is produced. Each of 
       these components is a \code{\link{list}} with 8 or 9 compnents. The component 
       named \code{data.plus}, being the input \code{\link{data.frame}} to which has been added 
       the columns required to fit the \code{TPPS} model (the \code{\link{data.frame}} stored 
       in the \code{data} component holds only the covariates from \code{data}).
       
       List of length 8 or 9 (according to the \code{asreml.option}).
 \enumerate{
    \item \code{data} = the input data frame augmented with structures required 
          to fit tensor product splines in \code{asreml-R}. This data frame can be used
          to fit the TPS model.

         Added columns:
       \itemize{
       \item \code{TP.col}, \code{TP.row} = column and row coordinates
       \item \code{TP.CxR} = combined index for use with smooth x smooth term
       \item \code{TP.C.n} for n=1:{diff.c} = X parts of column spline for use
             in random model (where diff.c is the order of column differencing)
       \item \code{TP.R.n} for n=1:{diff.r} = X parts of row spline for use in
             random model (where diff.r is the order of row differencing)
       \item \code{TP.CR.n} for n=1:{(diff.c*diff.r)} = interaction between the
             two X parts for use in fixed model. The first variate is 
             a constant term which should be omitted from the model when the constant 
             (1) is present. If all elements areincluded in the model then the 
             constant term should be omitted,eg. 
             \code{y ~ -1 + TP.CR.1 + TP.CR.2 + TP.CR.3 + TP.CR.4 + other terms...}
       \item when \code{asreml="grp"} or \code{"sepgrp"}, the spline basis 
             functions are also added into the data frame. Column numbers for each
             term are given in the \code{grp} list structure.
       }
   \item \code{mbflist} = list that can be used in call to asreml (so long as Z 
          matrix data frames extracted with right names, eg BcZ<stub>.df)
   \item \code{BcZ.df} = mbf data frame mapping onto smooth part of column
          spline, last column (labelled \code{TP.col}) gives column index
   \item \code{BrZ.df} = mbf data frame mapping onto smooth part of row spline,
         last column (labelled \code{TP.row}) gives row index
   \item \code{BcrZ.df} = mbf data frame mapping onto smooth x smooth term, last
         column (labelled \code{TP.CxR}) maps onto col x row combined index
   \item \code{dim} = list structure, holding dimension values relating to the
         model:
     \enumerate{
       \item \code{"diff.c"} = order of differencing used in column dimension
       \item \code{"nbc"} = number of random basis functions in column dimension
       \item \code{"nbcn"} = number of nested random basis functions in column dimension
              used in smooth x smooth term
       \item \code{"diff.r"} = order of differencing used in column dimension
       \item \code{"nbr"} = number of random basis functions in column dimension
       \item \code{"nbrn"} = number of nested random basis functions in column dimension
              used in smooth x smooth term
       }
  \item \code{trace} = list of trace values for ZGZ' for the random TPspline
        terms, where Z is the design matrix and G is the known diagonal variance
        matrix derived from eigenvalues. This can  be used to rescale the spline 
        design matrix (or equivalently variance components).
  \item \code{grp} = list structure, only added for setting \code{asreml="grp"}. 
        For \code{asreml="grp"}, provides column indexes for each of the 5 
        random components of the 2D splines in \code{data.plus}. Dimensions of 
        the components can be derived from the values in the \code{dim} item.  
  \item \code{data.plus} = the input \code{\link{data.frame}} to which has been  
       added the columns required to fit tensor product splines in \code{asreml-R}. 
       This \code{\link{data.frame}} can be used to fit the TPS model.  
       FOr multiple \code{sections}, this \code{\link{data.frame}} will occur 
       in the component for each section. If \code{asreml.option} is set to 
       \code{mbf}, then this component will have the \code{attribute} \code{mbf.env} 
       that specifies the environment to which the \code{\link{data.frame}s} 
       containing the spline bases have been assigned. 
 }
}

\references{Piepho, H-P, Boer, M. P. & Williams, E. R. (2022) Two-dimensional P-spline smoothing for spatial analysis of plant breeding trials. \emph{Biometrical Journal}, \bold{64}, 835-857.)

Rodriguez-Alvarez, M. X., Boer, M. P., van Eeuwijk, F. A., & Eilers, P. H. C. (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. \emph{Spatial Statistics}, \bold{23}, 52-71. 

Welham, S. J. (2022) \code{TPSbits}: \emph{Creates Structures to Enable Fitting and Examination of 2D Tensor-Product Splines using ASReml-R.} Version 1.0.0 \url{https://mmade.org/tpsbits/}} 

\author{Chris Brien}

\seealso{\code{\link{addSpatialModel.asrtests}}, \code{\link{addSpatialModelOnIC.asrtests}}, \cr 
         \code{\link{chooseSpatialModelOnIC.asrtests}}, \code{tpsmmb} from \code{TPSbits}}
\examples{\dontrun{

data(Wheat.dat)

#Add row and column covariates
Wheat.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })

#Set up the matrices
tps.XZmat <- makeTPPSplineMats(wheat.dat, 
                                row.covar = "cRow", col.covar = "cColumn")
}}
\keyword{asreml}
