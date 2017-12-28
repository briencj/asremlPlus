# asremlPlus
asremlPlus is an R package that augments the use of `ASReml-R` and `ASReml4-R` in fitting mixed models.

This version is compatible with both `ASReml-R` versions 3 and 4. `ASReml-R` version 4 is currently undergoing beta-testing and 
has some changes in syntax that necessitate changes in asremlPlus. 

Versions 4.x-xx of `asremlPlus` are a major revamp of the package and include substantial syntax changes. 

## More information

For more information install the package and run the R command `news(package = “asremlPlus”)` or consult the [manual](./inst/doc/asremlPlus-manual.pdf). 

An overview can be obtained using `?asremlPlus`. In particular, an example of its use is given towards the bottom of the help information.

## Installing the package

`asremlPlus` is an R package available on GitHub, so it can be installed from the RStudio console or an R command line session using the `devtools` command `install_github`. First, make sure `devtools` is installed, which, if you do not have it, can be done as follows:

`install.packages("devtools")`

Next, install `asremlPlus` from GitHub by entering:

`devtools::install_github("briencj\asremlPlus")`.

Version 2.0-12 of the package is available from CRAN so that you could first install it and its dependencies using:

`install.packages("asremlPlus")`


If you have not previously installed `asremlPlus` then you will need to install it dependencies:

`install.packages(c("dae", "ggplot2", "reshape", "plyr", "stringr", "RColorBrewer", `
`                   "foreach", "parallel", "doParallel"))`

<!--However, the version here may be somewhat newer than that on CRAN.-->

## What is does

It assists in automating the testing of terms in mixed models when 'asreml' is used 
to fit the models. A history of the fitting of a sequence of models is kept in a data frame. 
Procedures are available for choosing models that conform to the hierarchy or marginality principle 
and for displaying predictions for significant terms in tables and graphs. 

The content falls into the following natural groupings: 

(i) Data, 

(ii) Object manipulation functions, 

(iii) Model modification functions, 

(iv) Model testing functions, 

(v) Model diagnostics functions, 

(vi) Prediction production and presentation functions, 

(vii) Response transformation functions, and 

(viii) Miscellaneous functions. 
  
## What it needs  
  
To use `asremlPlus`, you must have a licensed version of either of the packages `asreml` and `asreml4`. 
They provide a computationally efficient algorithm for fitting mixed models using Residual Maximum 
  Likelihood. They can be purchased from 'VSNi' <http://www.vsni.co.uk/> as `asreml-R`, 
  who will supply a zip file for local installation/updating.
  
  It also imports [dae](<https://CRAN.R-project.org/package=dae>), [ggplot2](<https://CRAN.R-project.org/package=ggplot2>), `stats`, `methods`, `utils`, [reshape](<https://CRAN.R-project.org/package=reshape>), [plyr](<https://CRAN.R-project.org/package=plyr>), [stringr](<https://CRAN.R-project.org/package=stringr>), [RColorBrewer](<https://CRAN.R-project.org/package=RColorBrewer>), `grDevices`, 
[foreach](<https://CRAN.R-project.org/package=foreach>), [parallel](<https://CRAN.R-project.org/package=parallel>), [doParallel](<https://CRAN.R-project.org/package=doParallel>).
