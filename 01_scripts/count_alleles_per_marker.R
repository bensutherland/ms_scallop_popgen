# Count alleles per locus
# uses input as described in README from the stacks output microhap VCF
# B. Sutherland, 2023-03-31
setwd(dir = "~/Documents/pyes/stacks_workflow/")

input.FN <- "05-stacks/popn_out_microhaps/populations.haps_alleles.txt"
data.df <- read.delim2(file = input.FN, header = F, sep = "\t")
dim(data.df)
head(data.df)

# Add colnames
colnames(data.df) <- c("chr", "locus", "mname", "alleles")
head(data.df)

# Calculate the number of alleles observed (AN)
data.df$allele.count <- lengths(strsplit(as.character(data.df$alleles), ","))
head(data.df, n = 25)

write.table(x = data.df
            , file = "05-stacks/popn_out_microhaps/populations.haps_alleles_counts.txt"
            , sep = "\t"
            , row.names = F
            )
