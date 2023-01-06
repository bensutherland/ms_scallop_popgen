# simple_pop_stats analysis component of the Yesso scallop RADseq analysis
# B. Sutherland
# Initialized 2022-12-05
# Requires running "ms_scallop_popgen/01_scripts/pyes_popgen_simple_pop_stats.R" first

#### 01. Load Data ####
load(file = "03_results/post_all_filters.RData") # loaded from prerequisite script above

# Data is present in
obj

#### 04. Analysis ####
## Multivariate
# PCA from genind
pca_from_genind(data = obj
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/formatted_cols.csv"
                )

# Determine variance per axis
str(pca.obj) # $eig holds the absolute variance, this can be standardized
pca.obj$eig[1:10]
sum(pca.obj$eig)
(100 * pca.obj$eig / sum(pca.obj$eig))[1:10]
#16 / 251 * 100 # testing

# DAPC from genind
dapc_from_genind(data = obj
                 , plot_allele_loadings = TRUE
                 , colour_file = "00_archive/formatted_cols.csv"
                 , n.pca = 10, n.da = 2
                 , scree.da = TRUE
                 , posi.pca = "topright"
                 , scree.pca = TRUE
                 , dapc.width = 7
                 , dapc.height = 5
                 )


## Genetic differentiation
calculate_FST(format = "genind", dat = obj, separated = FALSE, bootstrap = TRUE)

## Private alleles (but recall that pa's were also evaluated by Stacks populations module)
regional_obj <- obj

# Combine related pops to query private alleles at regional level
unique(pop(regional_obj))
# pop(regional_obj) <- gsub(pattern = "VIU_offspring|VIU_parent", replacement = "VIU", x = pop(regional_obj)) # combine VIU
# pop(regional_obj) <- gsub(pattern = "PEN|FRA|JPN", replacement = "JPN", x = pop(regional_obj))              # combine JPN lineage
# unique(pop(regional_obj))

pa <- private_alleles(gid = regional_obj)
str(pa)
pa.t <- t(pa)
head(pa.t)
table(pa.t[,"BC"] > 0)
table(pa.t[,"JPN"] > 0)
table(pa.t[,"VIU"] > 0)

table(pop(regional_obj))

write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)

# Downsample and reassess
regional_obj_downsampled <- downsample_pops(data = regional_obj, subset_method = "min")
pa_downsampled <- private_alleles(gid = obj_subset)
pa_downsampled_t <- t(pa_downsampled)
head(pa_downsampled_t)
table(pa_downsampled_t[,"BC"] > 0)
table(pa_downsampled_t[,"JPN"] > 0)
table(pa_downsampled_t[,"VIU"] > 0)

table(pop(obj_subset))

write.csv(x = pa_downsampled, file = "03_results/private_alleles_downsampled.csv", quote = F)


## Inbreeding
# Estimating inbreeding (from adegenet tutorial)
obj_BC <- seppop(x = obj)$BC
obj_JPN <- seppop(x = obj)$JPN
obj_VIU <- seppop(x = obj)$VIU

# Use likelihood-based estimate of inbreeding to compute inbreeding coefficient of an individual (F)
# estimate inbreeding and return a sample of F values (# Note: warnings occur)
F_coeff_BC  <- inbreeding(x = obj_BC, N = 200)   # Calculates 100 values for each sample and outputs as a list
F_coeff_JPN <- inbreeding(x = obj_JPN, N = 200) 
F_coeff_VIU <- inbreeding(x = obj_VIU, N = 200) 

# ## plot the first 10 results (for first ten individuals) (example using BC)
# invisible(sapply(F_coeff_BC[1:10], function(e) plot(density(e)
#                                                     , xlab="F"
#                                                     , xlim=c(0,1)
#                                                     , main="Density of the sampled F values")
#                  )
#           )


pdf(file = "03_results/per_popn_mean_per_indiv_F_val.pdf", width = 6, height = 7)
par(mfrow=c(3,1))

## Compute means for all individuals
Fmean_BC=sapply(F_coeff_BC, mean)
hist(Fmean_BC, col="grey", xlab="mean value of F",
     main="Per-indiv average F (BC)"
     , xlim = c(0,1)
     , las = 1
)

text(x = 0.8, y = 10, label = paste0("mean = ", round(mean(sapply(F_coeff_BC, mean)), digits = 3)))

Fmean_JPN=sapply(F_coeff_JPN, mean)
hist(Fmean_JPN, col="grey", xlab="mean value of F",
     main="Per-indiv average F (JPN)"
     , xlim = c(0,1)
     , las = 1
)

text(x = 0.8, y = 15, label = paste0("mean = ", round(mean(sapply(F_coeff_JPN, mean)), digits = 3)))


Fmean_VIU=sapply(F_coeff_VIU, mean)
hist(Fmean_VIU, col="grey", xlab="mean value of F",
     main="Per-indiv average F (VIU)"
     , xlim = c(0,1)
     , las = 1
)

text(x = 0.8, y = 10, label = paste0("mean = ", round(mean(sapply(F_coeff_VIU, mean)), digits = 3)))
dev.off()

# Could potentially use related would be good to run after here
# uses function relatedness_calc.r

#### Relatedness ####
require("dartR")
require("related")
# run relatedness_calc.r interactively


# then
require(tidyr)

obj

relatedness_plot(file = "03_results/kinship_analysis_2023-01-05.Rdata", same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)

#relatedness_plot(file = "03_results/kinship_analysis_<date>.Rdata"
#, same_pops = TRUE
#, plot_by = "names")
#...where you can use either "names" or "codes" if using only same-on-same.
#...and if you set same_pops to FALSE, you will get all pops pairwise comparisons. (but can't use names)
file <- "03_results/kinship_analysis_2023-01-05.Rdata"
