\name{plotPredictions.data.frame}
\alias{plotPredictions.data.frame}
\alias{plotPredictions}
\title{Plots the predictions for a term, possibly with error bars.}
\description{This function plots the predictions \code{y} that are based on 
             \code{classify} and stored in the \code{data.frame data}. The package 
             \code{ggplot2} is used to produce the plots. Line plots are 
             produced when variables involving \code{x.num} or \code{x.fac} are 
             involved in \code{classify} for the predictions; otherwise, bar 
             charts are produced. Further, for line charts, the argument \code{panels} 
             determines whether a single plot or multiple plots in a single window 
             are produced; for bar charts, the argument \code{panels} is 
             ignored.}
\usage{
\method{plotPredictions}{data.frame}(data, classify, y, 
                x.num = NULL, x.fac = NULL, nonx.fac.order = NULL,  
                colour.scheme = "colour", panels = "multiple", 
                graphics.device = NULL,
                error.intervals = "Confidence",  titles = NULL, 
                y.title = NULL, filestem = NULL, ggplotFuncs = NULL, ...)
}
\arguments{
 \item{data}{a \code{data.frame} containing the values of the variables to be plotted. It 
                 should be consistent with an object of class \code{asremlPredict} or 
                 \code{asreml.predict} such as is stored in the \code{pvals} component 
                 of the value produced by \code{predict.asreml}; that is, in addition to 
                 variables classifying the predictions, it will include a column 
                 with the name specified in the \code{y} argument, usually 
                 \code{predicted.value} or \code{backtransformed.predictions}; 
                 each row contains a single predicted value. If \code{error.intervals} 
                 is not \code{"none"}, then the \code{predictions} component and, 
                 if present, the \code{backtransforms} component should contain 
                 columns for the lower and upper values of the limits for the interval 
                 with names that begin with \code{lower} and \code{upper}, 
                 respectively. The second part of the name must be one of 
                 \code{Confidence}, \code{StandardError} or 
                 \code{halfLeastSignificant}. The last part needs to be consistent 
                 between the lower and upper limits.}
 \item{classify}{a character string giving the combinations of the independent 
                 variables on which the predictions are based. It is 
                 an interaction type term formed from the independent variables, 
                 that is, separating the variable names with the \code{:} operator.}
 \item{y}{a character string giving the name of the variable that is to be plotted 
                 on the Y axis.}
 \item{x.num}{A character string giving the name of the numeric covariate that 
              corresponds to \code{x.fac}, is potentially included in terms in the 
              fitted model and which corresponds to the x-axis variable. It should 
              have the same number of unique values as the number of levels in 
              \code{x.fac}.}
 \item{x.fac}{A character string giving the name of the factor that corresponds to 
              \code{x.num}, is potentially included in terms in the fitted model and 
              which corresponds to the x-axis variable. It should have the same 
              number of levels as the number of unique values in \code{x.num}. 
              The levels of \code{x.fac} must be in the order in which they are to 
              be plotted - if they are dates, then they should be in the form 
              yyyymmdd, which can be achieved using \code{as.Date}. However, the levels
              can be non-numeric in nature, provided that \code{x.num} is also set.}
 \item{nonx.fac.order}{A character vector giving the order in which factors other 
              than \code{x.fac} are to be plotted in facetted plots (i.e. where the 
              number of non x factors is greater than 1). The first factor in the vector 
              will be plotted on the X axis (if there is no \code{x.num} or 
              \code{x.fac}. Otherwise, the order of plotting the factors is in 
              columns (X facets) and then rows (Y facets). By default the order is 
              in decreasing order for the numbers of levels of the non x factors.}
 \item{colour.scheme}{A character string specifying the colour scheme for the plots. 
              The default is \code{"colour"} which produces coloured lines and bars, 
              a grey background and white gridlines. A value of \code{"black"} 
              results in black lines, grey bars and gridlines and a white background.}
 \item{panels}{Possible values are \code{"single"} and \code{"multiple"}. 
              When line plots are to be produced, because variables involving 
              \code{x.num} or \code{x.fac} are involved in \code{classify} for 
              the predictions, \code{panels} determines whether or not a single 
              panel or multiple panels in a single window are produced. The 
              \code{panels} argument is ignored for bar charts.}
 \item{graphics.device}{A \code{character} specifying a graphics device for plotting. 
              The default is \cr\code{graphics.device = NULL}, which will result 
              in plots being produced on the current graphics device. Setting it to 
              \code{"windows"}, for example,  will result in a windows graphics 
              device being  opened.}
 \item{error.intervals}{A character string indicating the type of error interval, if any, 
             to plot in order to indicate uncertainty in the results. 
             Possible values are \code{"none"}, \code{"StandardError"}, \code{"Confidence"} 
             and \code{"halfLeastSignificant"}. Here, any option other than \code{"none"} 
             will result in the interval limits contained in \code{data} being plotted. }
 \item{titles}{A \code{list}, each component of which is named for a column in 
               the \code{data.frame} for the \code{asreml.obj} used in making the predictions 
               and contains a \code{character string} giving a title to use 
               in output (e.g. tables and graphs). Here they will be used for 
              axis labels for nonresponse variables. For response variable labels 
              see \code{y.title}.}
 \item{filestem}{A character sting giving the beginning of the name of the file in which to 
              save the plot. If \code{filestem = NULL}, the plot is not saved. The remainder 
              of the file name will be generated automatically and consists of the 
              following elements separated by full stops: the classify term, 
              \code{Bar} or \code{Line} and, if \code{error.intervals} 
              is not \code{"none"}, one of \code{SE}, \code{CI} or \code{LSI}. The
              file will be saved as a `png' file in the current work directory.}
 \item{y.title}{The title to be displayed on the y axis of any plot.} 
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.}
 \item{\dots}{further arguments passed to \code{ggplot}.}
}
\value{no values are returned.}

\author{Chris Brien}
\seealso{\code{\link{allDifferences.data.frame}}, \code{\link{predictPresent.asreml}}, 
         \code{\link{redoErrorIntervals.alldiffs}},  \cr 
         \code{\link{recalcLSD.alldiffs}}, \code{ggplot}, \code{Devices}}
\examples{\dontrun{
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE) #required for asreml4 only
current.asr <- asreml(fixed = log.Turbidity ~ Benches + Sources + Type + Species +
                                Sources:Type + Sources:Species + 
                                Sources:xDay + Species:xDay + Species:Date,
                      data = WaterRunoff.dat, keep.order = TRUE)
current.asrt <- asrtests(current.asr, NULL, NULL)

#### Get the observed combinations of the factors and variables in classify
class.facs <- c("Species","Date","xDay")
levs <- as.data.frame(table(WaterRunoff.dat[class.facs]))
levs <- as.list(levs[levs$Freq != 0, class.facs])
levs$xDay <- as.numfac(levs$xDay)

predictions <- predict(current.asr, classify="Species:Date:xDay", 
                       parallel = TRUE, levels = levs, 
                       present = c("Type","Species","Sources"))
#### for asreml3
predictions <- predictions$predictions$pvals
predictions <- predictions[predictions$est.status == "Estimable",]
#### for asreml4
predictions <- predictions$pvals
predictions <- predictions[predictions$status == "Estimable",]

plotPredictions(classify="Species:Date:xDay", y = "predicted.value", 
                data = predictions, wald.tab = current.asrt$wald.tab, 
                x.num = "xDay", x.fac = "Date", 
                x.title = "Days since first observation",
                y.title = "Predicted log(Turbidity)",
                present = c("Type","Species","Sources"),
                error.intervals = "none",
                ggplotFuncs = list(ggtitle("Transformed turbidity over time")))
  
diffs <- predictPlus(classify="Species:Date:xDay", 
                     present=c("Type","Species","Sources"), 
                     asreml.obj = current.asr, tables = "none", 
                     x.num = "xDay", x.fac = "Date", 
                     parallel = TRUE, levels = levs, 
                     x.plot.values=c(0,28,56,84),
                     wald.tab = current.asrt$wald.tab)
x.title <- "Days since first observation"
names(x.title) <- "xDay"
plotPredictions(classify="Species:Date:xDay", y = "predicted.value", 
               data = diffs$predictions, wald.tab = current.asrt$wald.tab, 
               x.num = "xDay", x.fac = "Date", 
               titles = x.title,
               y.title = "Predicted log(Turbidity)")
}}
\keyword{asreml}