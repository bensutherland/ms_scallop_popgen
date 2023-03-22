# simple_pop_stats analysis component of the Yesso scallop RADseq analysis
# B. Sutherland
# Initialized 2022-12-05
# Requires first running "ms_scallop_popgen/01_scripts/02_sps_char_and_filt.R"

#### 01. Load Data ####
load(file = "03_results/post_all_filters.RData") # loaded from prerequisite script above

# Data is present in
obj


#### 02. Compare MAF ####
# Create dataset per population
obj.VIU <- obj.sep$VIU
obj.JPN <- obj.sep$JPN

maf_filt(data = obj.VIU, maf = 0.01)
obj.VIU <- obj_maf_filt
myFreq.VIU <- myFreq

maf_filt(data = obj.JPN, maf = 0.01)
obj.JPN <- obj_maf_filt
myFreq.JPN <- myFreq

# How many variants and what percentage are under 0.1? 
table(myFreq.VIU < 0.1)
table(myFreq.JPN < 0.1)

table(myFreq.VIU < 0.1)[2] / length(myFreq.VIU) # 49.7%
table(myFreq.JPN < 0.1)[2] / length(myFreq.JPN) # 61.8%

# Plot
pdf(file = paste0("03_results/MAF_hist_VIU_JPN.pdf"), width = 7, height = 4)
par(mfrow=c(1,2))
hist(myFreq.VIU
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF (VIU)"
     , main = ""
     , ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 800, labels = paste("n = ", length(myFreq.VIU), " loci", sep = "" ))

hist(myFreq.JPN
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF (JPN)"
     , main = ""
     , ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 800, labels = paste("n = ", length(myFreq.JPN), " loci", sep = "" ))
dev.off()


#### 03. Multivariate analysis ####
obj

## Multivariate
# PCA from genind
pca_from_genind(data = obj
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/formatted_cols.csv"
                )

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

#### 04. Genetic differentiation and private alleles  ####
calculate_FST(format = "genind", dat = obj, separated = FALSE, bootstrap = TRUE)

## Private alleles
pa <- private_alleles(gid = obj)
pa.t <- t(pa)
head(pa.t)
table(pa.t[,"BC"] > 0)
table(pa.t[,"JPN"] > 0)
table(pa.t[,"VIU"] > 0)

table(pop(regional_obj))

write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)

#### 05. Relatedness and Inbreeding ####
obj

## Relatedness
# Calculate inter-individual relatedness
relatedness_calc(data = obj, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"

# Plot
relatedness_plot(file = "03_results/kinship_analysis_2023-01-06.Rdata", same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)

## Inbreeding
# Estimating inbreeding (from adegenet tutorial)
obj_BC  <- seppop(x = obj)$BC
obj_JPN <- seppop(x = obj)$JPN
obj_VIU <- seppop(x = obj)$VIU

# Use likelihood-based estimate of inbreeding to compute inbreeding coefficient of an individual (F)
# estimate inbreeding and return a sample of F values (# Note: warnings occur)
F_coeff_BC  <- inbreeding(x = obj_BC, N = 200)   # Calculates 100 values for each sample and outputs as a list
F_coeff_JPN <- inbreeding(x = obj_JPN, N = 200) 
F_coeff_VIU <- inbreeding(x = obj_VIU, N = 200) 

# Note: receiving warnings in all the above

# ## plot the first 10 results (for first ten individuals) (example using BC)
pdf(file = "03_results/per_popn_mean_per_indiv_F_val.pdf", width = 6, height = 7)
par(mfrow=c(3,1))

## Compute means for all individuals
Fmean_BC <- sapply(F_coeff_BC, mean) # Provides the average value per individual
hist(Fmean_BC, col="grey", xlab="mean value of F",
     main="Per-indiv average F (BC)"
     , xlim = c(0.4, 0.6)
     , las = 1
)

#text(x = 0.8, y = 5, label = paste0("mean = ", round(mean(sapply(F_coeff_BC, mean)), digits = 3)))

Fmean_JPN <- sapply(F_coeff_JPN, mean)
hist(Fmean_JPN, col="grey", xlab="mean value of F",
     main="Per-indiv average F (JPN)"
     , xlim = c(0.4, 0.6)
     , las = 1
)

#text(x = 0.8, y = 15, label = paste0("mean = ", round(mean(sapply(F_coeff_JPN, mean)), digits = 3)))

Fmean_VIU <- sapply(F_coeff_VIU, mean)
hist(Fmean_VIU, col="grey", xlab="mean value of F",
     main="Per-indiv average F (VIU)"
     , xlim = c(0.4, 0.6)
     , las = 1
)

#text(x = 0.8, y = 10, label = paste0("mean = ", round(mean(sapply(F_coeff_VIU, mean)), digits = 3)))
dev.off()


# single SNP per locus analysis is complete

# Comparison between VIU and JPN HOBS is available here: 
#"~/Documents/00_sutherland_bioinformatics/VIU_scallop/ms_scallop_popgen/01_scripts/per-locus_Hobs_compare_JPN_VIU.R"
# and will use "03_results/post_all_filters.RData"