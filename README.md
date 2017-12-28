# asremlPlus
asremlPlus is an R package that augments the use of 'ASReml-R' and 'ASReml4-R' in fitting mixed models

This version is compatible with both ASReml-R versions 3 and 4. ASReml-R version 4 is currently undergoing beta-testing and 
has some changes in syntax that necessitate changes in asremlPlus. 

Versions of asremlPlus 4.x-xx are a major revamp of the package and also include substantial syntax changes. 

## More information

For more information install the package and run the R command news(package = “asremlPlus”) or consult the [manual](./inst/doc/asremlPlus-manual.pdf). 

An overview can be obtained using `?asremlPlus`. In particular, an example of its use is given towards the bottom of the help information.

## Installing the package

1. If you do not already have it, `install(devtools)` (<https://CRAN.R-project.org/package=devtools>)
2. Execute the following in R: `devtools::install_github("briencj\asremlPlus")`

The package is also be available from CRAN 
(<https://CRAN.R-project.org/package=asremlPlus>). However, the version here may be somewhat newer than that on CRAN. 

## What is does

It assists in automating the testing of terms in mixed models when 'asreml' is used 
to fit the models. A history of the fitting of a sequence of models is kept in a data frame. 
Procedures are available for choosing models that conform to the hierarchy or marginality principle 
and for displaying predictions for significant terms in tables and graphs. 

The content falls into the following natural groupings: (i) Data, (ii) Object 
  manipulation functions, (iii) Model modification functions, (iv) Model testing functions, 
  (v) Model diagnostics functions, (vi) Prediction production and presentation functions, 
  (vii) Response transformation functions, and (viii) Miscellaneous functions. 
  
## What it needs  
  
You must have a licensed version of either of the packages 'asreml' and 'asreml4'. 
They provide a computationally efficient algorithm for fitting mixed models using Residual Maximum 
  Likelihood. They can be purchased from 'VSNi' <http://www.vsni.co.uk/> as 'asreml-R', 
  who will supply a zip file for local installation/updating.
  
  It also imports [dae](<https://CRAN.R-project.org/package=dae>), [ggplot2](<https://CRAN.R-project.org/package=ggplot2>), `stats`, `methods`, `utils`, [reshape](<https://CRAN.R-project.org/package=reshape>), [plyr](<https://CRAN.R-project.org/package=plyr>), [stringr](<https://CRAN.R-project.org/package=stringr>), [RColorBrewer](<https://CRAN.R-project.org/package=RColorBrewer>), `grDevices`, 
[foreach](<https://CRAN.R-project.org/package=foreach>), [parallel](<https://CRAN.R-project.org/package=parallel>), [doParallel](<https://CRAN.R-project.org/package=doParallel>).
