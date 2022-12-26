# simple_pop_stats component of the Yesso scallop RADseq analysis
# B. Sutherland
# Initialized 2022-12-05
# Requires running "ms_scallop_popgen/01_scripts/pyes_popgen_analysis.R" first

# Prior to running the following, source simple_pop_stats and choose Yesso scallop

#### 01. Load Data ####
load(file = "02_input_data/yesso_scallop_genind_2022-12-05.RData") # loaded from prerequisite script above
datatype <- "SNP" # normally assigned by load_genepop() when input is a genepop

obj <- my.data.gid
obj

#### 02. Prepare Data ####
unique(pop(obj))

characterize_genepop(obj)

## Define population colours
pops_in_genepop <- unique(pop(obj))
pops_in_genepop.df <- as.data.frame(pops_in_genepop)

## Download colours file from previous git repo to see previously used colours in Pacific oyster
# url = "https://raw.githubusercontent.com/bensutherland/ms_oyster_popgen/master/00_archive/my_cols.csv" # only need to run once
# destfile <- "00_archive/my_cols.csv"
# download.file(url, destfile)   # only need to run once
#my_colours <- read.csv(destfile)

# Manually write, based on the above where possible
pop_colours <- matrix(c("BC", "JPN", "VIU", "purple1", "turquoise4", "red"), nrow = 3, ncol = 2)
colnames(pop_colours) <- c("my.pops", "my.cols")

# Connect colours to empirical populations
colours <- merge(x = pops_in_genepop.df, y =  pop_colours, by.x = "pops_in_genepop", by.y = "my.pops"
                 #, sort = F
                 , all.x = T
)
colours


#### 03. Characterize missing data (indiv) and filter ####
##### 03.1 Individuals - missing data #####

# If file already exists, do not re-run percent_missing_data
if(file.exists("03_results/missing_data_per_indiv.csv")){
  
  print("Missing data information available, loading")
  
  missing_data.df <- read.csv(file = "03_results/missing_data_per_indiv.csv")
  
}else{
  
  print("Missing data information is not available, generating")
  
  percent_missing_by_ind(df = obj)
  
  write.csv(x = missing_data.df, file = "03_results/missing_data_per_indiv.csv", row.names = F)
  
}


# Provide population IDs to missing data, based on names
missing_data.df$pop <- rep(x = NA, times = nrow(missing_data.df))

# Provide population based on the individual name
missing_data.df$pop[grep(pattern = "BC", x = missing_data.df$ind)] <- "BC"
missing_data.df$pop[grep(pattern = "JPN", x = missing_data.df$ind)] <- "JPN"
missing_data.df$pop[grep(pattern = "VIU", x = missing_data.df$ind)] <- "VIU"
table(missing_data.df$pop)

head(missing_data.df)

# Combine colours to dataframe for plotting, don't sort, as it is still in the same order as the obj
colours
plot_cols.df <- merge(x = missing_data.df, y = colours, by.x = "pop", by.y = "pops_in_genepop", all.x = T
                      , sort = F
)

# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = 8, height = 5)
plot(1 - plot_cols.df$ind.per.missing, ylab = "Genotyping rate (%)"
     , col = plot_cols.df$my.cols
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,1)
)

abline(h = 0.7, lty = 3)

legend("bottomright", legend = unique(plot_cols.df$pop)
       , fill = unique(plot_cols.df$my.cols)
       , cex = 1.0
       , bg = "white"
)
dev.off()


## Filter individuals
# Identify which samples to retain based on genotyping rate
#  to keep inds with 70% genotyping rate (% missing < 0.3)
keep <- missing_data.df[missing_data.df$ind.per.missing < 0.3, "ind"]
print(paste0("Retaining ", length(keep), " of the total ", nInd(obj), " individuals"))

# Retain only the keep indiv
obj.filt <- obj[(keep)]
obj.filt
table(pop(obj.filt))





##### 03.4 Post-QC data filter #####
obj <- obj.filt
obj

## View the ind or loc names
inds <- indNames(obj)
loci <- locNames(obj)

# Save out which individuals have passed the filters
write.table(x = inds, file = "03_results/retained_individuals.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)

write.table(x = loci, file = "03_results/retained_loci.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)


##### 03.5 per marker stats and filters #####
## Per locus statistics
per_locus_stats(data = obj)
head(per_loc_stats.df)

pdf(file = "03_results/per_locus_Hobs.pdf", width = 6, height = 5) 
plot(x = per_loc_stats.df$Hobs
     , xlab = "Marker (index)"
     , ylab = "Observed Heterozygosity (Hobs)"
     , las = 1
)

abline(h = 0.5, lty = 3)
dev.off()


## Per locus, per population Hardy-Weinberg proportion statistics
hwe_eval(data = obj, alpha = 0.01)

hwe_outlier_mname_BC.vec    <- per_locus_hwe_BC.df[per_locus_hwe_BC.df$`Pr(chi^2 >)` < 0.01, "mname"]
hwe_outlier_mname_JPN.vec    <- per_locus_hwe_JPN.df[per_locus_hwe_JPN.df$`Pr(chi^2 >)` < 0.01, "mname"]
hwe_outlier_mname_VIU.vec    <- per_locus_hwe_VIU.df[per_locus_hwe_VIU.df$`Pr(chi^2 >)` < 0.01, "mname"]

# How many outliers (p < 0.01) per population
length(hwe_outlier_mname_BC.vec)    #  705 markers out of HWE
length(hwe_outlier_mname_JPN.vec)   # 1780 markers out of HWE
length(hwe_outlier_mname_VIU.vec)   # 1028 markers out of HWE

# How many unique HWE deviating markers?  
markers_to_drop <- unique(c(hwe_outlier_mname_BC.vec, hwe_outlier_mname_JPN.vec, hwe_outlier_mname_VIU.vec))
length(markers_to_drop) # 2595 ( the sum total of each is 3513, so there are ~1000 markers that are seen twice)

markers_to_keep <- setdiff(x = locNames(obj), y = markers_to_drop)
length(markers_to_keep) # 6,804 markers to keep

obj <- obj[, loc=markers_to_keep]
obj

## Remove Hobs > 0.5 markers 
# Which markers are greater than 0.5 heterozygosity?
hobs.outliers <- per_loc_stats.df[per_loc_stats.df$Hobs > 0.5, "mname"] # 210 markers

# How many hobs outliers remain after samples were dropped for HWE deviation?
hobs.outliers.remaining <- intersect(hobs.outliers, locNames(obj))
length(hobs.outliers.remaining) # 73 remain, should drop these too

keep <- setdiff(x = locNames(obj), y = hobs.outliers)
length(keep) # 6,731 remain

# Drop Hobs > 0.5 loci from genind
obj <- obj[, loc=keep] # 
obj


# Post hwe filter, could sep pop, calc per_locus_stats per pop, then correlate Hobs
sep.obj <- seppop(x = obj)



##### 03.5 Post-all filters #####
# Save out colours to be used downstream
colours
colnames(x = colours) <- c("collection", "colour")
write.csv(x = colours, file = "00_archive/formatted_cols.csv", quote = F, row.names = F)


#### 04. Analysis ####
## Multivariate
# PCA from genind
pca_from_genind(data = obj, PCs_ret = 4, colour_file = "00_archive/formatted_cols.csv")

# DAPC from genind
dapc_from_genind(data = obj, plot_allele_loadings = TRUE, colour_file = "00_archive/formatted_cols.csv", n.pca = 10, n.da = 2)

## Genetic differentiation
calculate_FST(format = "genind", dat = obj, separated = FALSE, bootstrap = TRUE)

## Private alleles
regional_obj <- obj

# Combine related pops to query private alleles at regional level
unique(pop(regional_obj))
pop(regional_obj) <- gsub(pattern = "VIU_offspring|VIU_parent", replacement = "VIU", x = pop(regional_obj)) # combine VIU
pop(regional_obj) <- gsub(pattern = "PEN|FRA|JPN", replacement = "JPN", x = pop(regional_obj))              # combine JPN lineage
unique(pop(regional_obj))

pa <- private_alleles(gid = regional_obj)
write.csv(x = pa, file = "03_results/private_alleles.csv", quote = F)


####### Convert genepop to Rubias format #####
# Need to create a stock code file, in the form of
# in the tab-delim format of: 
#collection	repunit
#12Mile_Creek	GoA

stock_code.df <- as.data.frame(unique(pop(obj)))
colnames(stock_code.df) <- "collection"
stock_code.df$repunit <- stock_code.df$collection
stock_code.df
write_delim(x = stock_code.df, file = "00_archive/stock_code.txt", delim = "\t", col_names = T)
micro_stock_code.FN <- "00_archive/stock_code.txt"
# this is for annotate_rubias(), for an unknown reason it requires the name micro_stock_code.FN

## Convert genepop to rubias
obj # the current analysis object

## If running manually, here are the arguments needed
#sample_type <- "reference"
#data <- obj

genepop_to_rubias_SNP(data = obj, sample_type = "reference", custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN)

# Using this output, move to "01_scripts/ckmr_from_rubias.R"

## Simulations
# full_sim(rubias_base.FN = "03_results/rubias_output_SNP.txt", num_sim_indiv = 200, sim_reps = 100)

## Inbreeding
# # Estimating inbreeding (from adegenet tutorial)
# obj_PEN <- seppop(x = obj)$PEN
# obj_VIU_parent <- seppop(x = obj)$VIU_parent
# obj_VIU_offspring <- seppop(x = obj)$VIU_offspring
# obj_DPB <- seppop(x = obj)$DPB
# 
# # compute the mean inbreeding for each individual and plot
# #temp <- inbreeding(x = obj_PEN, N = 100)
# #temp <- inbreeding(x = obj_VIU_parent, N = 100)
# #temp <- inbreeding(x = obj_VIU_offspring, N = 100)
# temp <- inbreeding(x = obj_DPB, N = 100)
# 
# class(temp)
# head(names(temp))
# temp[[1]] # temp is a list of values sampled from the likelihood distribution of each individual; means values are obtained for all indiv using sapply
# Fbar <- sapply(temp, mean)
# hist(Fbar, col = "firebrick", main = "Average inbreeding in Pendrell")
# hist(Fbar, col = "firebrick", main = "Average inbreeding in VIU parents")
# hist(Fbar, col = "firebrick", main = "Average inbreeding in VIU offspring")
# hist(Fbar, col = "firebrick", main = "Average inbreeding in DPB")


## Per sample heterozygosity

# The following would need extensive coding to make happen
#rubias_to_vcf() # write out, then use instructions here to get per individual heterozygosity in vcftools
# https://github.com/bensutherland/ms_oyster_popgen/blob/master/01_scripts/heterozygosity.sh
# per population heterozygosity


# related would be good to run after here

