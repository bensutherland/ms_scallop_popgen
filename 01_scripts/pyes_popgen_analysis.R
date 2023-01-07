# Main popgen analysis script for Yesso scallop 
# Initialized 2022-12-05, B. Sutherland

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("adegenet")

library("rstudioapi")
library("adegenet")


# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/ms_scallop_popgen\\/01_scripts", replacement = "", x = current.path)
current.path <- paste0(current.path, "/stacks_workflow")
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
output.dir <- "12-results"
input.FN <- "05-stacks/populations_single_snp.raw"


#### 1. Import data ####
print(paste0("Loading data from ", input.FN))
my.data <- read.PLINK(file = input.FN)
my.data

# ## Explore data attributes
# nInd(my.data) # How many individuals?
# nLoc(my.data) # How many loci?
# indNames(my.data)[1:20] # What individuals?
# nPop(my.data) # How many pops?
# unique(pop(my.data)) # What pops?
# locNames(my.data, withAlleles = T)[1:20] # What allele names?


# Update population names
pop(my.data) <- gsub(pattern = 1, replacement = "BC", x = pop(my.data))
pop(my.data) <- gsub(pattern = 2, replacement = "VIU", x = pop(my.data))
pop(my.data) <- gsub(pattern = 3, replacement = "JPN", x = pop(my.data))
# unique(pop(my.data)) # What pops?


#### 2. Data exploration (allele frequency) ####
# Plot instances of minor allele across individuals and loci
png(file = paste0(output.dir, "/", "glPlot.png"), width = 924, height = 600)
glPlot(x = my.data, posi="topleft") 
dev.off()

# Create density plot of minor allele frequencies
pdf(file = paste0(output.dir, "/", "maf_hist.pdf"), width = 6, height = 4)
myFreq <- glMean(my.data)
hist(myFreq
     #, proba=T # note: does not sum to 1, not worth using
     , col="gold", xlab = "Minor allele frequency (MAF)"
     , main = ""
     , ylim = c(0, 2500)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
     )
text(x = 0.4, y = 1500, labels = paste("n = ", nLoc(my.data), " loci", sep = "" ))
dev.off()


#### 3. Convert genlight to genind ####
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

# Next go to "ms_scallop_popgen/01_scripts/pyes_popgen_simple_pop_stats.R"