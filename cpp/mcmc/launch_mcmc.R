#!/usr/bin/env Rscript

# Script header
header <- "#!/bin/sh
#SBATCH --mem=500
#SBATCH -q fast"

# Script footer
footer <- "exit 0"

# Parameters
library(dplyr, quietly = T)

for (a in c("yes", "no")) {
  for (m in c("model1", "model2")) {
    for (curr in 1:100) {
      # Redirect output and error messages
      reOut <- paste0("#SBATCH -o mcmc_sim_", m, "_", curr, "_", a, ".log")
      reErr <- paste0("#SBATCH -e mcmc_sim_", m, "_", curr, "_", a, ".err")
      
      # Command to run the script
      command <- paste("srun program.out 1 simulation", curr, m, a)
      
      # Create script to submit
      scriptToLaunch <- paste0("mcmc_sim_", m, "_", curr, "_", a, ".sh")
      
      script <- paste0(header, "\n",
                       reOut, "\n", 
                       reErr, "\n\n",
                       command, "\n\n",
                       # paste0("rm ", scriptToLaunch, "\n"),
                       footer)
      
      write(script, scriptToLaunch)
      
      # Submit to cluster queue
      system(paste("sbatch", scriptToLaunch)) 
    }  
  }  
}


