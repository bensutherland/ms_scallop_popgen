# Main popgen analysis script for Yesso scallop 
# Initialized 2022-12-05, B. Sutherland

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
install.packages("rstudioapi")
install.packages("adegenet")

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


#### 2. Minor allele freq ####
# Plot instances of minor allele across individuals and loci
png(file = paste0(output.dir, "/glPlot_all.png"), width = 924, height = 600)
glPlot(x = my.data, posi="topleft") 
dev.off()

# Create density plot of minor allele frequencies
pdf(file = paste0(output.dir, "/maf_hist_all.pdf"), width = 6, height = 4)
myFreq <- glMean(my.data)
hist(myFreq, proba=T, col="gold", xlab = "Allele frequencies"
     , main = ""
     , ylim = c(0,50)
     , ylab = "Density of second allele frequencies"
     )
text(x = 0.4, y = 7, labels = paste(nLoc(my.data), " loci", sep = "" ))
temp <- density(myFreq)
lines(temp$x, temp$y, lwd=3)
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

# Data is now a genind, and therefore can be used with simple_pop_stats


