# Analyze population-specific per-locus Hobs
# Sutherland, B.J.G. 2022-12-26

# Assumes that you have a filtered R object for loading after running "01_scripts/pyes_popgen_simple_pop_stats.R"
load(file = "03_results/post_all_filters.RData")

# Post filter, could sep pop, calc per_locus_stats per pop, then correlate Hobs
sep.obj <- seppop(x = obj
                  #, drop = TRUE   # note: drop = T removes alleles that are not observed in the pop
                  ) 

dir.create(path = "03_results/popn_sp_hobs_calc")

# Back up the original
dir.create(path = "03_results/all_popn_hobs_calc")
filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
file.copy(from = paste0("03_results/", filename) , to = paste0("03_results/all_popn_hobs_calc/", filename))

# Per population, run for
# JPN
drop_loci(sep.obj$JPN, drop_monomorphic = T) # Remove monomorphs; # 5857 polymorphic loci
per_locus_stats(data = obj_filt) # note: since there is only one popn, expect no Fit or Fst calc
per_loc_stats_JPN.df <- per_loc_stats.df

filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
filename <- tail(filename, n = 1) # take only the most recent
file.copy(from = paste0("03_results/", filename) , to = "03_results/popn_sp_hobs_calc/per_locus_stats_JPN.txt")

# VIU
drop_loci(sep.obj$VIU, drop_monomorphic = T) # Remove monomorphs # 4856 polymorphic loci
per_locus_stats(data = obj_filt)
per_loc_stats_VIU.df <- per_loc_stats.df
filename <- basename(Sys.glob("03_results/per_locus_stats_*"))
filename <- tail(filename, n = 1) # take only the most recent
file.copy(from = paste0("03_results/", filename) , to = "03_results/popn_sp_hobs_calc/per_locus_stats_VIU.txt")

# Restore the original
filename <- basename(Sys.glob("03_results/all_popn_hobs_calc/per_locus_stats_*"))
filename <- head(filename, n = 1)
file.copy(from = paste0("03_results/all_popn_hobs_calc/", filename) , to = paste0("03_results/", filename), overwrite = TRUE)
per_loc_stats.df <- read.delim(file = Sys.glob("03_results/all_popn_hobs_calc/per_locus_stats_*"))
head(per_loc_stats.df)
# TODO: this could be done better #

# Compare Hobs
head(per_loc_stats_JPN.df)
colnames(per_loc_stats_JPN.df)[which(colnames(per_loc_stats_JPN.df)=="Hobs")] <- "Hobs.JPN"
dim(per_loc_stats_JPN.df)
head(per_loc_stats_VIU.df)
colnames(per_loc_stats_VIU.df)[which(colnames(per_loc_stats_VIU.df)=="Hobs")] <- "Hobs.VIU"
dim(per_loc_stats_VIU.df)

all_per_loc_data.df <- merge(x = per_loc_stats_JPN.df, y = per_loc_stats_VIU.df, by = "mname")
head(all_per_loc_data.df)
dim(all_per_loc_data.df) 
# note: 4,224 polymorphic loci shared between both populations
# however, we will run with all loci, as HOBS = 0 is still a value

# Plot
pdf(file = "03_results/per_locus_hobs_VIU_vs_JPN.pdf", width = 5, height = 5)
plot(x = all_per_loc_data.df$Hobs.JPN, y = all_per_loc_data.df$Hobs.VIU
     , xlab = expression(per ~ locus ~ H[OBS] ~ (VIU))
     , ylab = expression(per ~ locus ~ H[OBS] ~ (JPN))
     , las = 1
     )
abline(h = 0.3, lty = 2)
abline(h = 0.5, lty = 2)
abline(v = 0.3, lty = 2)
abline(v = 0.5, lty = 2)
dev.off()

# How many 0.3 < HOBS < 0.5 in JPN && VIU?
nrow(all_per_loc_data.df[all_per_loc_data.df$Hobs.JPN > 0.3 
                         & all_per_loc_data.df$Hobs.JPN < 0.5
                         & all_per_loc_data.df$Hobs.VIU > 0.3
                         & all_per_loc_data.df$Hobs.VIU < 0.5
                         , ]) # 538

# TODO: merge all stats with HOBS VIU and JPN including HOBS global # 


# Lots of shared polymorphism to select from for the top markers
