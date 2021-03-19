# asremlPlus

[![Project Status: Active:  The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-2.10.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/asremlPlus)](https://cran.r-project.org/package=asremlPlus)
[![packageversion](https://img.shields.io/badge/Package%20version-4.2--32-orange.svg?style=flat-square)](/commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--03--19-yellowgreen.svg)](/commits/master)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)
[![Downloads](https://cranlogs.r-pkg.org/badges/last-week/asremlPlus)](commits/master)


`asremlPlus` is an R package that augments the use of `ASReml-R` in fitting mixed models and packages generally in exploring prediction differences. This version is compatible with both `ASReml-R` versions 3 and 4.1, but not 4.0.

Versions 4.x-xx of `asremlPlus` are a major revamp of the package and include substantial syntax changes. In particular, most functions are S3 methods and so the type of the object can be omitted from the function name when calling the function.  

## More information

For more information install the package and run the R command `news(package = “asremlPlus”)` or consult the [manual](./vignettes/asremlPlus-manual.pdf). 

An overview can be obtained using `?asremlPlus`. In particular, an example of its use is given towards the bottom of the help information and this is avalable as the [Wheat.analysis vignette](./vignettes/Wheat.analysis.pdf). It that shows how to select the terms to be included in a mixed model for an experiment that involves spatial variation; it also illustrates diagnostic checking and prediction production and presentation for this example. A second vignette is the [Wheat.infoCriteria vignette](./vignettes/Wheat.infoCriteria.pdf) that illustrates the facilities in `asremlPlus` for producing and using information criteria. Two further vignettes show how to use `asremlPlus` for exploring and presenting predictions from a linear mixed model analysis in the context of a three-factor factorial experiment on ladybirds: one vignette, [Ladybird.asreml vignette](./vignettes/Ladybird.asreml.pdf), uses `asreml` and `asremlPlus` to produce and present  predictions; the other vignette, [Ladybird.lm vignette](./vignettes/Ladybird.asreml.pdf), uses `lm` to produce the predictions and `asremlPlus` to present the predictions. The vignettes can be accessed via `vignette(name, package = "asremlPlus")`, where `name` is one of `"Wheat.analysis"`, `"Wheat.infoCriteria"`, `"Ladybird.asreml"` or `"Ladybird.lm"`.

## Installing the package

### From a repository using `drat`

Windows binaries and source tarballs of the latest version of `asremlPlus` are available for installation from my [repository](http://chris.brien.name/rpackages). Installation instructions are available there.

### Directly from  GitHub

`asremlPlus` is an R package available on GitHub, so it can be installed from the RStudio console or an R command line session using the `devtools` command `install_github`. First, make sure `devtools` is installed, which, if you do not have it, can be done as follows:

`install.packages("devtools")`

Next, install `asremlPlus` from GitHub by entering:

`devtools::install_github("briencj/asremlPlus")`.

Version 2.0-12 of the package is available from CRAN so that you could first install it and its dependencies using:

`install.packages("asremlPlus")`

If you have not previously installed `asremlPlus` then you could first install it and its dependencies from CRAN using:

`install.packages("asremlPlus")`

Otherwise, you will need to install its dependencies manually:

`install.packages(c("dae", "ggplot2", "reshape", "plyr", "dplyr", "stringr", "RColorBrewer", `
`                   "foreach", "parallel", "doParallel"))`

## What is does

It assists in automating the testing of terms in mixed models when `asreml-R` is used 
to fit the models. A history of the fitting of a sequence of models is kept in a data frame. 
Procedures are available for choosing models that conform to the hierarchy or marginality principle. It can also be used to display, in tables and graphs, predictions obtained from a mixed model using your favourite model fitting functions and to explore differences between predictions. As a general rule functions that are methods for `asreml` and `asrtests` objects require `asreml-R`; on the other hand, functions that are methods for `alldiffs` and `data.frame` objects do not require `asreml-R`.

The use of the package is exemplified in four vignettes: [Wheat.analysis vignette](./vignettes/Wheat.analysis.pdf), [Wheat.infoCriteria vignette](./vignettes/Wheat.infoCriteria.pdf), [Ladybird asreml vignette](./vignettes/Ladybird.asreml.pdf) and [Ladybird lm vignette](./vignettes/Ladybird.asreml.pdf). They can be accessed via `vignette(name, package = "asremlPlus")`, where `name` is one of `"Wheat.analysis"`, `"Wheat.infoCriteria"`, `"Ladybird.asreml"` or `"Ladybird.lm"`.

The content falls into the following natural groupings: 

(i) Data, 

(ii) Object manipulation functions, 

(iii) Model modification functions, 

(iv) Model testing functions, 

(v) Model diagnostics functions, 

(vi) Prediction production and presentation functions, 

(vii) Response transformation functions, and 

(viii) Miscellaneous functions. 

For a list of the functions for each group, see the help for `asremlPlus-package` or the entry in the manual for `asremlPlus-package`.  
  
## What it needs  
  
To use those functon in `asremlPlus` that are methods for `asreml` or `asrtests` objects, you must have a licensed version of the package `asreml`. It provides a computationally efficient algorithm for fitting mixed models using Residual Maximum Likelihood. It can be purchased from 'VSNi' <http://www.vsni.co.uk/> as `asreml-R`, who will supply a zip file for local installation/updating.
  
  It also imports [dae](<https://CRAN.R-project.org/package=dae>), [ggplot2](<https://CRAN.R-project.org/package=ggplot2>), `stats`, `methods`, `utils`, [reshape](<https://CRAN.R-project.org/package=reshape>), [plyr](<https://CRAN.R-project.org/package=plyr>), [dplyr](<https://CRAN.R-project.org/package=dplyr>), [stringr](<https://CRAN.R-project.org/package=stringr>), [RColorBrewer](<https://CRAN.R-project.org/package=RColorBrewer>), `grDevices`, 
[foreach](<https://CRAN.R-project.org/package=foreach>), [parallel](<https://CRAN.R-project.org/package=parallel>), [doParallel](<https://CRAN.R-project.org/package=doParallel>).

## License

The `asremlPlus` package is distributed under the [MIT licence](<https://opensource.org/licenses/MIT>) -- for details see [LICENSE.md](https://github.com/briencj/asremlPlus/blob/master/LICENSE.md).
