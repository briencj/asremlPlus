
cat("#### Test for parallel processing with asreml42\n")
test_that("Parallel_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(tictoc)
  Sys.setenv("OMP_NUM_THREADS" = 1)
  
  #'## load a data set consisting of 1056 observations
  data("ChickpeaEnd.dat")
  n <- 1000
  
  (cl <- makeCluster(detectCores()))
  registerDoParallel(cl)
  if (requireNamespace("tictoc", quietly = TRUE))
    tictoc::tic(paste("parallel", n, "analyses"))
  setTimeLimit(elapsed = 900)
  fits <- foreach (i = 1:n, .packages = c("asreml","asremlPlus"))  %dopar%
    { 
      current.asr <- asreml(Biomass ~ Smarthouse + Genotypes*Treatments , 
                            random = ~ Smarthouse:(Lane + Position),
                            residual = ~ dsum(~ar1(Lane):ar1(Position) | Smarthouse), 
                            data=ChickpeaEnd.dat)
    }
  stopCluster(cl)
  if (requireNamespace("tictoc", quietly = TRUE))
  { 
    timing <- tictoc::toc(log = TRUE)
    testthat:::expect_true((timing$toc - timing$tic) < 420)
  }
})  
