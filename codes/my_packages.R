suppressMessages({
  
  load_lib <- c("dplyr", "plyr", "xlsx", "readxl", "reshape", "fs", "rmarkdown",
                "reticulate", "tidyr","tidyverse", "ggplot2", "pheatmap", 
                "RColorBrewer", "viridis", "FactoMineR", "BiocManager", 
                "magrittr", "ggalluvial", "ggfortify", "factoextra", "ggExtra",
                "paletteer", "svglite")
  
  install_lib <- load_lib[!(load_lib %in% installed.packages())] # check package
  if (length(install_lib)) for (i in install_lib) utils::install.packages(i) # install
  
  cat("Loaded Packages:\n")
  print(sapply(load_lib, require, character = TRUE)) # load
  cat("\n\n")
  
})
