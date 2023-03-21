# Import plink data and convert to genind
# Initialized 2022-12-05, B. Sutherland

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("adegenet")

library("rstudioapi")
library("adegenet")


# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/ms_scallop_popgen\\/01_scripts", replacement = "", x = current.path)
current.path <- paste0(current.path, "/stacks_workflow")
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
input.FN <- "05-stacks/popn_out_single_snp/populations_single_snp.raw"

#### 1. Import data ####
print(paste0("Loading data from ", input.FN))
my.data <- read.PLINK(file = input.FN)
my.data

# Update population names
pop(my.data) <- gsub(pattern = 1, replacement = "BC", x = pop(my.data))
pop(my.data) <- gsub(pattern = 2, replacement = "VIU", x = pop(my.data))
pop(my.data) <- gsub(pattern = 3, replacement = "JPN", x = pop(my.data))
# unique(pop(my.data)) # What pops?

#### 2. Convert genlight to genind ####
# Convert genlight to matrix
my.data.mat <- as.matrix(my.data)
my.data.mat[1:5,1:5]

# Translate the number of minor allele to genind format
my.data.mat[my.data.mat == 0] <- "1/1" #homozygote reference
my.data.mat[my.data.mat == 1] <- "1/2" #heterozygote
my.data.mat[my.data.mat == 2] <- "2/2" #homozygote alternate
my.data.mat[1:5,1:5]

# Convert matrix to genind
my.data.gid <- df2genind(my.data.mat, sep = "/", ploidy = 2) # convert df to genind

# Transfer pop attributes
pop(my.data.gid) <- pop(my.data) 
unique(pop(my.data.gid))

# Data is now a genind, and therefore can be used with simple_pop_stats
save(my.data.gid, file="../ms_scallop_popgen/03_results/prepared_genind.RData")

# Next go to "ms_scallop_popgen/01_scripts/02_sps_char_and_filt.R"