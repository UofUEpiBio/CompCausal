###################
## load packages ##
###################
library(tidyverse)
library(splines)
library(gridExtra)
library(parallel)

source(paste0(getwd(), "/HelperFunction.R"))
source(paste0(getwd(), "/singleindexmodelfunctions.R"))
source(paste0(getwd(), "/est_s_t_y.R"))
source(paste0(getwd(), "/est_exchange.R"))
source(paste0(getwd(), "/SIDR_Ravinew.R"))
source(paste0(getwd(), "/SIDRnew.R"))
source(paste0(getwd(), "/SensIAT_sim_outcome_modeler_mave.R"))


# Preparing the data
iteration_list <- 1:2000
set.seed(312)
seeds <- sample.int(1e7, length(iteration_list))

# Creating the cluster
cl <- makeCluster(10)
clusterExport(cl, c("seeds", "est_s_t_y_create_containers", "IF_trunc_func", "SIM", "est_psi", 
                    "gam.variables", "single.index.variables", "eq", "fit_one_analysis", "pboot", 
                    "fit_SensIAT_single_index_norm1coef_model", "SIDR_Ravinew", "SIDRnew", 
                    "SIDRnew_fixed_bandwidth", "cumuSIR", "SIDRnew_h2", "NW_new"))
clusterEvalQ(cl, {
  library(tidyverse)
  library(splines)
  library(gridExtra)
  library(parallel)
})


# Running the simulation
ans <- parLapply(cl, iteration_list, \(i) {
  
  # Check if simulation result exist
  fn_i <- file.path("results", paste0("result_", i, ".rds") )
  
  if (file.exists(fn_i))
    return(readRDS(fn_i))
  
  ans_i <- tryCatch({
    
    # Grab the seed for the individual simulation
    set.seed(seeds[i])
    
    # Run the simulation
    rnorm(1)
  }, error = \(e) e)
  
  # Save if it didn't failed
  if (!inherits(ans_i, "error"))
    saveRDS(ans_i, fn_i)

return(ans_i)

})

# Just in case some result failed
saveRDS(ans, "ans.rds")

# Combining results (assuming all went well)
ans <- do.call(rbind, ans)

stopCluster(cl)

