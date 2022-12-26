# Analyze population-specific per-locus Hobs
# Sutherland, B.J.G. 2022-12-26

# Assumes that you have a filtered R object for loading after running "01_scripts/pyes_popgen_simple_pop_stats.R"
load(file = "03_results/post_all_filters.RData")

# Post hwe filter, could sep pop, calc per_locus_stats per pop, then correlate Hobs
sep.obj <- seppop(x = obj)

dir.create(path = "03_results/popn_sp_hobs_calc")

# Back up the original
dir.create(path = "03_results/all_popn_hobs_calc")
filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
file.copy(from = paste0("03_results/", filename) , to = paste0("03_results/all_popn_hobs_calc/", filename))

# Per population, run for
# JPN
per_locus_stats(data = sep.obj$JPN)  # Not clear why getting NAs for all loci, Fit and Fst
per_loc_stats_JPN.df <- per_loc_stats.df

filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
file.copy(from = paste0("03_results/", filename) , to = "03_results/popn_sp_hobs_calc/per_locus_stats_JPN.txt")

# VIU
per_locus_stats(data = sep.obj$VIU)
per_loc_stats_VIU.df <- per_loc_stats.df
filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
file.copy(from = paste0("03_results/", filename) , to = "03_results/popn_sp_hobs_calc/per_locus_stats_VIU.txt")

# Restore the original
filename <- basename(Sys.glob("03_results/all_popn_hobs_calc/per_locus_stats_*"))
file.copy(from = paste0("03_results/all_popn_hobs_calc/", filename) , to = paste0("03_results/", filename), overwrite = TRUE)

# Compare Hobs
head(per_loc_stats_JPN.df)
colnames(per_loc_stats_JPN.df)[which(colnames(per_loc_stats_JPN.df)=="Hobs")] <- "Hobs.JPN"
head(per_loc_stats_VIU.df)
colnames(per_loc_stats_VIU.df)[which(colnames(per_loc_stats_VIU.df)=="Hobs")] <- "Hobs.VIU"

all_per_loc_data.df <- merge(x = per_loc_stats_JPN.df, y = per_loc_stats_VIU.df, by = "mname")
head(all_per_loc_data.df)

# Plot
plot(x = all_per_loc_data.df$Hobs.JPN, y = all_per_loc_data.df$Hobs.VIU)

# Lots of shared polymorphism to select from for the top markers
