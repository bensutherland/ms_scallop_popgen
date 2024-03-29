# simple_pop_stats analysis component of the Yesso scallop RADseq analysis
# B. Sutherland
# Initialized 2022-12-05
# Requires first running "ms_scallop_popgen/01_scripts/02_sps_char_and_filt.R"

#install.packages("ggpubr")
require("ggpubr")

#### 01. Load Data ####
load(file = "03_results/post_all_filters.RData") # loaded from prerequisite script above

# Data is present in
obj


#### 02. Compare MAF ####
# Create dataset per population
obj.VIU <- obj.sep$VIU
obj.JPN <- obj.sep$JPN
obj.BC  <- obj.sep$BC

# Filter to remove markers with population MAF < 0.01
maf_filt(data = obj.VIU, maf = 0.01)
obj.VIU <- obj_maf_filt
myFreq.VIU <- myFreq

maf_filt(data = obj.JPN, maf = 0.01)
obj.JPN <- obj_maf_filt
myFreq.JPN <- myFreq

maf_filt(data = obj.BC, maf = 0.01)
obj.BC <- obj_maf_filt
myFreq.BC  <- myFreq

# How many variants and what percentage are between 0.01 and 0.1? 
table(myFreq.VIU < 0.1)
table(myFreq.JPN < 0.1)
table(myFreq.BC < 0.1)

table(myFreq.VIU < 0.1)[2] / length(myFreq.VIU) # 49.7%
table(myFreq.JPN < 0.1)[2] / length(myFreq.JPN) # 61.8%
table(myFreq.BC < 0.1)[2] / length(myFreq.BC)   # 39.8%

# Plot
pdf(file = paste0("03_results/MAF_hist_popn_sp.pdf"), width = 7.4, height = 3.5)
par(mfrow=c(1,3))
hist(myFreq.JPN
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF (JPN)"
     , main = ""
     , ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
abline(v = 0.1, lty = 3)
text(x = 0.4, y = 800, labels = paste("n = ", length(myFreq.JPN), " loci", sep = "" ))
text(x = 0.375, y = 700
     , labels = paste0("MAF<0.1: "
                       , round(x = as.numeric(table(myFreq.JPN < 0.1)[2] / length(myFreq.JPN)) * 100
                               , digits = 1)
                       , "%"))

hist(myFreq.VIU
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF (VIU)"
     , main = ""
     , ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
abline(v = 0.1, lty = 3)
text(x = 0.4, y = 800, labels = paste("n = ", length(myFreq.VIU), " loci", sep = "" ))
text(x = 0.375, y = 700
     , labels = paste0("MAF<0.1: "
                       , round(x = as.numeric(table(myFreq.VIU < 0.1)[2] / length(myFreq.VIU)) * 100
                              , digits = 1)
                       , "%"))

hist(myFreq.BC
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF (BC farm)"
     , main = ""
     , ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
abline(v = 0.1, lty = 3)
text(x = 0.4, y = 800, labels = paste("n = ", length(myFreq.BC), " loci", sep = "" ))
text(x = 0.375, y = 700
     , labels = paste0("MAF<0.1: "
                       , round(x = as.numeric(table(myFreq.BC < 0.1)[2] / length(myFreq.BC)) * 100
                               , digits = 1)
                       , "%"))

dev.off()


#### 03. Per locus statistics #####
### Save out previous runs of per_locus_stats, where pops were not separated
dir.create(path = "03_results/all_popn_per_loc_stats")
filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
file.copy(from = paste0("03_results/", filename) , to = paste0("03_results/all_popn_per_loc_stats/", filename))

### Make new working directory for population-specific per locus stats
dir.create(path = "03_results/popn_sp_hobs_calc")

### Population-specific per locus stats
# JPN
obj.JPN
per_locus_stats(data = obj.JPN) # note: since there is only one popn, expect no Fit or Fst calc
per_loc_stats_JPN.df <- per_loc_stats.df
head(per_loc_stats_JPN.df)
filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
filename <- tail(filename, n = 1) # take only the most recent
file.copy(from = paste0("03_results/", filename) , to = "03_results/popn_sp_hobs_calc/per_locus_stats_JPN.txt")

# VIU
obj.VIU
per_locus_stats(data = obj.VIU)
per_loc_stats_VIU.df <- per_loc_stats.df
filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
filename <- tail(filename, n = 1) # take only the most recent
file.copy(from = paste0("03_results/", filename) , to = "03_results/popn_sp_hobs_calc/per_locus_stats_VIU.txt")

## Direct comparison compare Hobs
# Rename cols for clarity
colnames(per_loc_stats_JPN.df)[which(colnames(per_loc_stats_JPN.df)=="Hobs")] <- "Hobs.JPN"
colnames(per_loc_stats_VIU.df)[which(colnames(per_loc_stats_VIU.df)=="Hobs")] <- "Hobs.VIU"

all_per_loc_data.df <- merge(x = per_loc_stats_JPN.df, y = per_loc_stats_VIU.df, by = "mname")
head(all_per_loc_data.df)
dim(all_per_loc_data.df) 
# note: 3,539 polymorphic loci shared between both populations

# Plot
pdf(file = "03_results/popn_sp_hobs_calc/per_locus_hobs_VIU_vs_JPN.pdf", width = 5, height = 5)
plot(x = all_per_loc_data.df$Hobs.JPN, y = all_per_loc_data.df$Hobs.VIU
     , xlab = expression(per ~ locus ~ H[OBS] ~ (VIU))
     , ylab = expression(per ~ locus ~ H[OBS] ~ (JPN))
     , las = 1
)
abline(h = 0.5, lty = 2)
abline(v = 0.5, lty = 2)
dev.off()

## Plot 2d density plot
# example as basic scatterplot
ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) +
  geom_point()

# # Contour only
# ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) + 
#   geom_density_2d()

# # with colour and level
# ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) + 
#   stat_density_2d(aes(fill = ..level..), geom = "polygon")
  
# # Contour only, coloured by level
# ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) + 
#   stat_density_2d(aes(color = ..level..))

# # Contour only, coloured by level, modify BINS
# ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) + 
#   stat_density_2d(bins = 20, aes(color = ..level..))

# Contour only, coloured by level, modify BINS, add points
plot <-ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) + 
        stat_density_2d(bins = 17, aes(color = ..level..)) + 
        geom_point(alpha = 0.05) + 
        theme_bw() + 
        xlab(bquote(H[OBS]*" VIU")) +
        ylab(bquote(H[OBS]*" JPN"))
plot

pdf(file = "03_results/popn_sp_hobs_calc/per_locus_hobs_VIU_vs_JPN_contours.pdf", width = 6, height = 5)
plot
dev.off()

# ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) + 
#   geom_bin2d() +
#   theme_bw()
# 
# ggplot(data = all_per_loc_data.df, aes(x = Hobs.VIU, y = Hobs.JPN)) + 
#   geom_bin2d(bins = 20) +
#   scale_fill_continuous(type = "viridis") +
#   theme_bw()

# How many loci in the specific range of interest? 
table(all_per_loc_data.df$Hobs.JPN > 0.3 & all_per_loc_data.df$Hobs.VIU > 0.3) # 756
table(all_per_loc_data.df$Hobs.JPN < 0.2 & all_per_loc_data.df$Hobs.VIU < 0.2) # 1,209


## Write out results
write.csv(x = all_per_loc_data.df, file = "03_results/popn_sp_hobs_calc/popn_sp_JPN_VIU_HOBS.csv", quote = F)

#### 03. Multivariate analysis ####
obj

## Multivariate
# For an unknown reason, sps currently requires a
#   manual sourcing of the function to properly use the retain_pca_obj function
source("~/Documents/pyes/simple_pop_stats/01_scripts/utilities/pca_from_genind.r", echo=TRUE)

# PCA from genind
pca_from_genind(data = obj
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/formatted_cols.csv"
                )

file.copy(from = "03_results/pca_scores_per_sample.txt", to = "03_results/pca_scores_per_sample_sibs_incl.txt")

## Prepare an eigenvalue plot for inset
num_eigenvals <- 10
vals.df <- as.data.frame(pca.obj$eig[1:num_eigenvals])
colnames(vals.df)[1] <- "vals"
vals.df$pc <- seq(1:num_eigenvals)
vals.df
colnames(vals.df) <- c("PVE", "PC")

# Express eigenvalues as a percentage of total variation explained
tot.var <- sum(pca.obj$eig)
vals.df$PVE <- vals.df$PVE/tot.var *100

# Barplot
eig.plot <- ggplot(data = vals.df, aes(x=PC, y=PVE)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_blank() #remove x axis labels
        , axis.ticks.x=element_blank() #remove x axis ticks
        , panel.background = element_blank()
  )


## Plot
# Remove legend pc1 v pc2
pc1_v_pc2.plot  <- pc1_v_pc2.plot + theme(legend.position = "none")
pc1_v_pc2.plot  <- pc1_v_pc2.plot + annotation_custom(ggplotGrob(eig.plot)
                                                       , xmin = 1, xmax = 5
                                                       , ymin = -10, ymax = -3.5
                    )

# Legend inside panel second plot
pc3_v_pc4.plot <- pc3_v_pc4.plot + theme(legend.justification = c(1,0), legend.position = c(1,0)
                                         , legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")
                                )

final.figure <- ggarrange(pc1_v_pc2.plot, pc3_v_pc4.plot
                  , labels = c("A", "B")
                  , ncol = 2, nrow = 1
                  )

pdf(file = "03_results/pca_composite_figure_w_close_kin.pdf", width = 12, height = 6.5)
print(final.figure)
dev.off()


#### 04. Genetic differentiation and private alleles  ####
# This is moved down now, because it should be calculated with putative close-kin removed
# calculate_FST(format = "genind", dat = obj, separated = FALSE, bootstrap = TRUE)

## Private alleles
pa <- private_alleles(gid = obj)
pa.t <- t(pa)
head(pa.t)
table(pa.t[,"BC"] > 0)  #   58
table(pa.t[,"JPN"] > 0) # 1401
table(pa.t[,"VIU"] > 0) #  270

write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)

## Regional private alleles
regional_obj <- obj
unique(pop(regional_obj))
pop(regional_obj) <- gsub(pattern = "VIU|BC", replacement = "CAN", x = pop(regional_obj))
#pop(obj) # to check

pa <- private_alleles(gid = regional_obj)
pa.t <- t(pa)
head(pa.t)
table(pa.t[,"CAN"] > 0) # 690
table(pa.t[,"JPN"] > 0) # 1401

pa.df <- as.data.frame(pa.t)
head(pa.df)
CAN.pa <- rownames(pa.df[pa.df$CAN!=0,])
CAN.pa <- gsub(pattern = "\\..*", replacement = "", x = CAN.pa)
head(CAN.pa)

JPN.pa <- rownames(pa.df[pa.df$JPN!=0,])
JPN.pa <- gsub(pattern = "\\..*", replacement = "", x = JPN.pa)
head(JPN.pa)

# Compare frequency of regional pa to all loci
maf_filt(data = obj, maf = 0) # calculate MAF, no cutoff keeps 6547 variants

# Obtain frequencies of PA
myFreq.JPN.pa <- myFreq[names(myFreq) %in% JPN.pa]
length(myFreq.JPN.pa) # 1401
myFreq.CAN.pa <- myFreq[names(myFreq) %in% CAN.pa]
length(myFreq.CAN.pa) # 690

# Obtain frequencies of non-PA
myFreq.other <- myFreq[!(names(myFreq) %in% c(JPN.pa, CAN.pa))]
length(myFreq.other) # 4456

# Plot
pdf(file = paste0("03_results/MAF_hist_pa.pdf"), width = 7.4, height = 3.5)
par(mfrow=c(1,3))
hist(myFreq.other
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF, shared"
     , main = ""
     , ylim = c(0, 700)
     , xlim = c(0, 0.5)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 450, labels = paste("n = ", length(myFreq.other), " loci", sep = "" ))
text(x = 0.375, y = 400, labels = paste0("median = ", round(median(myFreq.other), digits = 3)))

hist(myFreq.JPN.pa
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF (JPN, private alleles)"
     , main = ""
     , ylim = c(0, 700)
     , xlim = c(0,0.5)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 450, labels = paste("n = ", length(myFreq.JPN.pa), " loci", sep = "" ))
text(x = 0.375, y = 400, labels = paste0("median = ", round(median(myFreq.JPN.pa), digits = 3)))

hist(myFreq.CAN.pa
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF (CAN, private alleles)"
     , main = ""
     , ylim = c(0, 700)
     , xlim = c(0, 0.5)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 450, labels = paste("n = ", length(myFreq.CAN.pa), " loci", sep = "" ))
text(x = 0.375, y = 400, labels = paste0("median = ", round(median(myFreq.CAN.pa), digits = 3)))

dev.off()


#### 05. Relatedness ####
obj

## Relatedness
# Calculate inter-individual relatedness
relatedness_calc(data = obj, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"
gc()

# Plot
relatedness_plot(file = "03_results/kinship_analysis_2023-09-22.Rdata", same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)

# Inspect relatedness
head(output$relatedness)
dim(output$relatedness) # 14,878 pairs
output_reduced.df <- output$relatedness[output$relatedness$group=="BCBC" 
                                        | output$relatedness$group=="JPJP"
                                        | output$relatedness$group=="VIVI", ]

dim(output_reduced.df)  # 5,231 pairs

# Write out the same-on-same output
write.csv(x = output_reduced.df, file = "03_results/relatedness_results_same_only.csv", quote = F
          , row.names = F)

# Use function to identify the putative close kin
source("~/Documents/pyes/simple_pop_stats/01_scripts/dev/id_close_kin.R", echo=TRUE)
id_close_kin(cutoff = 0.25, statistic = "ritland")
drop.list

drop.inds <- c(drop.list$BCBC, drop.list$VIVI)

obj
table(pop(obj))

keep.inds <- setdiff(x = indNames(obj), y = drop.inds)
obj_sibs_purged <- obj[keep.inds]
obj_sibs_purged
table(pop(obj_sibs_purged))

maf_filt(data = obj_sibs_purged, maf = 0.01)

obj_maf_filt

#### 06. Redone analyses after sibs purged ####
# PCA from genind
pca_from_genind(data = obj_maf_filt
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/formatted_cols.csv"
)

file.copy(from = "03_results/pca_scores_per_sample.txt", to = "03_results/pca_scores_per_sample_sibs_purged.txt")

## Prepare an eigenvalue plot for inset
num_eigenvals <- 10
vals.df <- as.data.frame(pca.obj$eig[1:num_eigenvals])
colnames(vals.df)[1] <- "vals"
vals.df$pc <- seq(1:num_eigenvals)
vals.df
colnames(vals.df) <- c("PVE", "PC")

# Express eigenvalues as a percentage of total variation explained
tot.var <- sum(pca.obj$eig)
vals.df$PVE <- vals.df$PVE/tot.var *100

# Barplot
eig.plot <- ggplot(data = vals.df, aes(x=PC, y=PVE)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_blank() #remove x axis labels
        , axis.ticks.x=element_blank() #remove x axis ticks
        , panel.background = element_blank()
  )


## Plot
# Remove legend pc1 v pc2
pc1_v_pc2.plot  <- pc1_v_pc2.plot + theme(legend.position = "none")
pc1_v_pc2.plot  <- pc1_v_pc2.plot + annotation_custom(ggplotGrob(eig.plot)
                                                      , xmin = -4, xmax = -0.5
                                                      , ymin = -10, ymax = -5
)

# Legend inside panel second plot
pc3_v_pc4.plot <- pc3_v_pc4.plot + theme(legend.justification = c(1,0), legend.position = c(1,0)
                                         , legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")
)

final.figure <- ggarrange(pc1_v_pc2.plot, pc3_v_pc4.plot
                          , labels = c("A", "B")
                          , ncol = 2, nrow = 1
)

pdf(file = "03_results/pca_composite_figure_close_kin_purged.pdf", width = 12, height = 6.5)
print(final.figure)
dev.off()

# FST recalculated
calculate_FST(format = "genind", dat = obj_maf_filt, separated = FALSE, bootstrap = TRUE)

# Per locus stats recalculated
dir.create(path = "03_results/purged_sibs_per_locus_stats")
per_locus_stats(data = obj_maf_filt)
file.copy(from = "03_results/per_locus_stats_2023-09-22.txt", to = "03_results/purged_sibs_per_locus_stats/per_locus_stats_purged_sibs.txt")


# ## Inbreeding
# # Estimating inbreeding (from adegenet tutorial)
# obj_BC  <- seppop(x = obj)$BC
# obj_JPN <- seppop(x = obj)$JPN
# obj_VIU <- seppop(x = obj)$VIU
# 
# # Use likelihood-based estimate of inbreeding to compute inbreeding coefficient of an individual (F)
# # estimate inbreeding and return a sample of F values (# Note: warnings occur)
# F_coeff_BC  <- inbreeding(x = obj_BC, N = 200)   # Calculates 100 values for each sample and outputs as a list
# F_coeff_JPN <- inbreeding(x = obj_JPN, N = 200) 
# F_coeff_VIU <- inbreeding(x = obj_VIU, N = 200) 
# 
# # Note: receiving warnings in all the above
# 
# # ## plot the first 10 results (for first ten individuals) (example using BC)
# pdf(file = "03_results/per_popn_mean_per_indiv_F_val.pdf", width = 6, height = 7)
# par(mfrow=c(3,1))
# 
# ## Compute means for all individuals
# Fmean_BC <- sapply(F_coeff_BC, mean) # Provides the average value per individual
# hist(Fmean_BC, col="grey", xlab="mean value of F",
#      main="Per-indiv average F (BC)"
#      , xlim = c(0.4, 0.6)
#      , las = 1
# )
# 
# #text(x = 0.8, y = 5, label = paste0("mean = ", round(mean(sapply(F_coeff_BC, mean)), digits = 3)))
# 
# Fmean_JPN <- sapply(F_coeff_JPN, mean)
# hist(Fmean_JPN, col="grey", xlab="mean value of F",
#      main="Per-indiv average F (JPN)"
#      , xlim = c(0.4, 0.6)
#      , las = 1
# )
# 
# #text(x = 0.8, y = 15, label = paste0("mean = ", round(mean(sapply(F_coeff_JPN, mean)), digits = 3)))
# 
# Fmean_VIU <- sapply(F_coeff_VIU, mean)
# hist(Fmean_VIU, col="grey", xlab="mean value of F",
#      main="Per-indiv average F (VIU)"
#      , xlim = c(0.4, 0.6)
#      , las = 1
# )
# 
# #text(x = 0.8, y = 10, label = paste0("mean = ", round(mean(sapply(F_coeff_VIU, mean)), digits = 3)))
# dev.off()


# single SNP per locus analysis is complete

# Comparison between VIU and JPN HOBS is available here: 
#"~/Documents/00_sutherland_bioinformatics/VIU_scallop/ms_scallop_popgen/01_scripts/per-locus_Hobs_compare_JPN_VIU.R"
# and will use "03_results/post_all_filters.RData"