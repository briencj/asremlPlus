.onAttach <- function(...)
{ 
  if (!any(c("asreml","asreml4") %in% loadedNamespaces()))
  { 
#    if ("asreml4" %in% .packages(all.available = TRUE))
#      asreml.loaded <- requireNamespace("asreml4", quietly=TRUE)
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
  tips <- c("Need help? The manual is in the doc subdirectory of the package's install directory.", 
            "Find out what has changed in asremlPlus: enter news(package = 'asremlPlus').",
            "Need help getting started? Look at the example in ?`asremlPlus-package`.", 
            "To avoid start-up message that ASReml-R is needed, load asreml before asremlPlus.",
            "Use suppressPackageStartupMessages() to eliminate all package startup messages.", 
            "To see all the intermittent, randomly-presented, startup tips enter ?asremlPlusTips.",
            "For versions between CRAN releases (and more) go to http://chris.brien.name/rpackages.")
  tip <- sample(tips, 1)
  packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
}
