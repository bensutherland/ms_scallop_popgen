# To predict the number of cut sites using an in silico digest

# Install packages
install.packages("BiocManager")
require("BiocManager")
BiocManager::install("ShortRead")
BiocManager::install("SimRAD")
require("ShortRead")
require("SimRAD")

# for GC percent
install.packages("seqinr")
require(seqinr)

#### Import the reference genome ####
# Yesso scallop
# GenBank
rfsq <- ref.DNAseq(FASTA.file = "/Users/wayne/Documents/00_sutherland_bioinformatics/VIU_scallop/00_genome/GCA_002113885.2_ASM211388v2_genomic.fna.gz"
           , subselect.contigs = T
           , prop.contigs = 0.1
           #, prop.contigs = 1
           )

# Testing with Pacific oyster
# GenBank
rfsq <- ref.DNAseq(FASTA.file = "~/Documents/genomes/Cgig/GCA_000297895.1_oyster_v9_genomic.fa"
                   , subselect.contigs = T
                   , prop.contigs = 0.1
                   #, prop.contigs = 1
)

# What is the size of the ref genome sampled?
width(rfsq) # 10% = 109,980,884 bp
            # 10% of Pacific oyster = 52,798,009 bp

# What is the GC content:
GC(s2c(rfsq)) # 10% = 36.5 %
              # 10% of Pacific oyster is 33.4%

##### ddRAD Set your cutters: ####
# https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities

# ddRAD w/ NsiI and MspI
# NsiI = ATGCA/T
cs_5p1 <- "ATGCA"
cs_3p1 <- "T"

# MspI = C/CGG
cs_5p2 <- "C"
cs_3p2 <- "CGG"

# GO TO NEXT STEP #
# OR #

# ddRAD w/ PstI and MspI
# PstI = CTGCA/G
cs_5p1 <- "CTGCA"
cs_3p1 <- "G"

# MspI = C/CGG
cs_5p2 <- "C"
cs_3p2 <- "CGG"

# GO TO NEXT STEP #
# OR #

# ddRAD w/ SbfI and MspI
# SbfI = CCTGCA/GG
cs_5p1 <- "CCTGCA"
cs_3p1 <- "GG"

# MspI = C/CGG
cs_5p2 <- "C"
cs_3p2 <- "CGG"


#### Simulate digest ####
rfsq.dig <- insilico.digest(rfsq, cs_5p1, cs_3p1, cs_5p2, cs_3p2
                              , verbose = T
                              )

# From digest, select only those reads digested by the two enzymes
rfsq.sel <- adapt.select(rfsq.dig, type="AB+BA"
                         , cs_5p1, cs_3p1, cs_5p2, cs_3p2)

# From those selected w/ correct cut sites, select based on size
# 100-250 bp (standard at IBIS, not including barcodes etc)
wid.rfsq.sel <- size.select(rfsq.sel,  min.size = 100, max.size = 250
                            , graph=TRUE, verbose=TRUE)

