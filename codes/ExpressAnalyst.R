# Rascunho


## ExpressAnalystR attempt

library(pacman)
library(devtools)

pacman::p_load(igraph, RColorBrewer, qs, rjson, RSQLite)

# Step 2: Install ExpressAnalystR WITHOUT documentation
devtools::install_github("xia-lab/ExpressAnalystR", build = TRUE, 
                         build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))

# Step 2: Install ExpressAnalystR WITH documentation
devtools::install_github("xia-lab/ExpressAnalystR", build = TRUE,
                         build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)


# Skipping 6 packages not available: multtest, metagenomeSeq, sva, limma, edgeR, fgsea

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("sva", "multtest", "metagenomeSeq", 'limma', 'edgeR', 'fgsea'))


library(ExpressAnalystR)

setwd("C:/Users/Windows 10/Documents/Faculdade/IC_dengue/interferome-how-to/ExpressAnalyst/Teste1/")
getwd()

#boolean FALSE indicates local mode (vs web mode);
Init.Data(FALSE);

# Set analysis type to meta-analysis
SetAnalType("metadata");

#Read dataset text file
dataSets <- ReadOmicsData("E-GEOD-25713.txt");
dataSets <- SanityCheckData("E-GEOD-25713.txt");

#Map gene id to entrez id
dataSets <- AnnotateGeneData("E-GEOD-25713.txt", "mmu", "entrez");