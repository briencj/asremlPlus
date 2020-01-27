\name{chooseModel.data.frame}
\alias{chooseModel.data.frame}
\title{Determines the set of significant terms from results stored in a \code{data.frame}, 
       taking into account the marginality relations of terms and recording the tests used 
       in a \code{\link{data.frame}}.}
\description{Uses the \code{p.values} from a set of hypothesis tests that are stored in 
  the supplied \code{\link{data.frame}} to choose a model to describe the effects of the 
  terms corresponding to the p-values, taking into account the hierarchy or marginality 
  of terms. In particular, a term will not be tested if it is marginal to (or nested in) 
  one that is significant. For example, if A:B is significant, then neither A nor B will 
  be tested. The tests used in choosing the selected model are listed in the 
  \code{\link{data.frame}} \code{choose.summary}. 

 No change is made to the \code{p.values}, the \code{DF} and \code{denDF} being for 
 information only.
}
\usage{
\method{chooseModel}{data.frame}(object, terms=NULL, p.values = "Pr", 
            DF = "Df", denDF = "denDF", omit.DF = FALSE, 
            terms.marginality=NULL, alpha = 0.05, ...)
}
\arguments{
 \item{object}{a \code{\link{data.frame}} object containing the results of hypothesis tests for a 
              set of terms. It components shoud include \code{terms}, \code{p.values}, 
              and, if not set to \code{NA}, \code{DF} and \code{denDF}.}
 \item{terms}{A \code{character} giving the name of the column in \code{object} containing the 
              terms corresponding to the \code{p.values}. If \code{NULL}, it is assumed that 
              the row names of \code{object} give the terms.}
 \item{p.values}{A \code{character} giving the name of the column in \code{object} containing 
              the p-values to use in deciding whether or not terms are significant.}
 \item{DF}{Can be a \code{character} or a \code{numeric} that specifies the numerator 
              degrees of freedom for the tests. 
              If it is a \code{character}, it must be the name of a column in \code{object} 
              containing the numerator degrees of freedom that are to be included in the 
              \code{choose.summary} \code{data.frame}. if it is a \code{numeric}, its length 
              must equal 1 or the number of rows in \code{object}. In either case, a column labelled 
              \code{DF} will be included in the \code{choose.summary} \code{data.frame}. It will 
              contain either the replicated single value (which can be \code{NA}) or the values supplied.}
 \item{denDF}{Can be a \code{character} or a \code{numeric} that specifies the denominator 
              degrees of freedom for the tests. 
              If it is a \code{character}, it must be the name of a column in \code{object} 
              containing the denominator degrees of freedom that are to be included in the 
              \code{choose.summary} \code{data.frame}. if it is a \code{numeric}, its length 
              must equal 1 or the number of rows in \code{object}. In either case, a column labelled 
              \code{denDF} will be included in the \code{choose.summary} \code{data.frame}. It will 
              contain either the replicated single value (which can be \code{NA}) or the values supplied.}
 \item{omit.DF}{A \code{logical} indicating whether or not both the numerator and denominator 
              degrees of freedom are to be omitted from \code{choose.summary}. Doing so will 
              mean that the \code{choose.summary} no longer has the same columns as a 
              \code{test.summary} from an \code{\link{asrtests.object}}.}
 \item{terms.marginality}{A square matrix of ones and zeros with row and column names 
              being the names of the those terms in the \code{terms} column of \code{object} 
              that are to be tested. The diagonal elements should be one, indicating 
              that a term is marginal to itself. Elements should be one if the row 
              term is marginal to the column term. All other elements should be zero. 
              The names of the rows and columns should match the those elements of 
              \code{terms} that are to be tested.}
 \item{alpha}{The significance level for the hypothesis testing.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A list containing:
       \enumerate{
         \item \code{choose.summary}: a \code{\link{data.frame}} summarizing 
               the tests carried out in choosing the significant terms; 
               provided \code{omit.DF =} \code{FALSE}, it has the same columns as a 
              \code{test.summary} from an \code{\link{asrtests.object}}
         \item \code{sig.tests}: a \code{character vector} whose elements are the 
                the significant terms amongst those tested.}
}

\author{Chris Brien}
\seealso{\code{\link{chooseModel}}, \code{\link{chooseModel.asrtests}}}
\examples{  
  data("Ladybird.dat")
  
  ## Use asreml to get the table of p-values

  \dontrun{
  m1.asr <- asreml(logitP ~ Host*Cadavers*Ladybird, 
                   random = ~ Run,
                   data = Ladybird.dat)
  current.asrt <- as.asrtests(m1.asr)
  fixed.tab <-  current.asrt$wald.tab
  col.p <- "Pr"
  df = "Df"
  den.df = "denDF"
  }
  
  ## Use lmeTest to get the table of p-values
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(logitP ~ Host*Cadavers*Ladybird + (1|Run),
                              data=Ladybird.dat)
    fixed.tab <- anova(m1.lmer, type = "II")
    col.p <- "Pr(>F)"
    df = "NumDF"
    den.df = "DenDF"
  }
  
  ## Select a model using the table of p-values obtained with either asreml or lmerTest
  if (exists("fixed.tab"))
  {
    term.marg <- dae::marginality(dae::pstructure(~ Host*Cadavers*Ladybird, 
                                                  data = Ladybird.dat))
    chosen <- chooseModel(fixed.tab, p.values = col.p, DF = df, denDF = den.df, 
                          terms.marginality = term.marg)
  }
}
\keyword{htest}