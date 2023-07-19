.onAttach <- function(...)
{ 
  if (!any(c("asreml","asreml41") %in% loadedNamespaces()))
  { 
#    if ("asreml4" %in% .packages(all.available = TRUE))
#      asreml.loaded <- requireNamespace("asreml41", quietly=TRUE)
#    else
#      asreml.loaded <- requireNamespace("asreml", quietly=TRUE)
#    if (!asreml.loaded)  
    { 
packageStartupMessage("ASReml-R needs to be loaded if the mixed-model functions are to be used.

ASReml-R is available from VSNi. Please visit http://www.vsni.co.uk/ for more information.\n")
    }
  }
  
  if (!interactive() || sample.int(2, 1) == 1) 
    return()
  tips <- c("Need help? The manual is a vignette and is in the vignettes subdirectory of the package's install directory.", 
            "Find out what has changed in asremlPlus: enter news(package = 'asremlPlus').",
            "Need help getting started? Enter vignette(package = 'asremlPlus').", 
            "To avoid start-up message that ASReml-R is needed, load asreml before asremlPlus.",
            "The methods for alldiffs and data.frame do not require asreml",
            "Use suppressPackageStartupMessages() to eliminate all package startup messages.", 
            "To see all the intermittent, randomly-presented, startup tips enter ?asremlPlusTips.",
            "To install the latest version: go to http://chris.brien.name/rpackages.",
            "For versions between CRAN releases (and more) go to http://chris.brien.name/rpackages.")
  tip <- sample(tips, 1)
  packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
}
