# Scallop Population Genomics
This repository is in development stage, and comes with no guarantees. 

### Requirements
`Stacks v2` http://catchenlab.life.illinois.edu/stacks/     
`cutadapt`  https://cutadapt.readthedocs.io/en/stable/index.html     
`bwa` http://bio-bwa.sourceforge.net/   
`samtools` http://www.htslib.org/    
`R` https://www.r-project.org     
`FastQC` https://www.bioinformatics.babraham.ac.uk/projects/fastqc/    
`multiQC` https://multiqc.info/      
`stacks_workflow` https://github.com/enormandeau/stacks_workflow    
`plink` https://zzz.bwh.harvard.edu/plink/plink2.shtml      
`fineRADstructure` https://www.milan-malinsky.org/fineradstructure        


_note: see R scripts for individual package requirements_

All commands are executed from `stacks_workflow`. This repo should be cloned in the same parent folder as `stacks_workflow`, as some scripts are used from here.      

## 0. Project planning
Before starting, plan the appropriate RADseq enzymes for use with the scallop, considering target number of loci and relation to sequencing depth.     
Use the following script, after updating the locations for the scallop and Pacific oyster genome (as a reference) reference genomes:       
`01_scripts/scallop_genome_digest.R`       

## 1. Prepare data
### a. Setup
1. Copy links to all raw data in `02-raw` (use cp -l)    
2. Prepare the tab-delimited `sample_information.csv` as explained in `stacks_workflow` (also see template sample_information.csv).     
3. Download reference genome: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_002113885.1/      
note: source citation is Wang et al. 2017, Nat Ecol Evol. https://pubmed.ncbi.nlm.nih.gov/28812685/     

### b. Quality control and trim
View raw data with fastqc and multiqc:    
```
mkdir 02-raw/fastqc_raw    
fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 5    
multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw   
```
Prepare lane_info.txt file with automated script:    
`./00-scripts/00_prepare_lane_info.sh`    

Trim adapters and too short reads:    
`./00-scripts/01_cutadapt.sh <numCPUs>`    

View trimmed data with fastqc and multiqc:     
```
mkdir 02-raw/trimmed/fastqc_trimmed/    
fastqc -t 5 02-raw/trimmed/*.fastq.gz -o 02-raw/trimmed/fastqc_trimmed/
multiqc -o 02-raw/trimmed/fastqc_trimmed/ 02-raw/trimmed/fastqc_trimmed       
```

### c. De-multiplex reads
Detect cut sites and barcodes to de-multiplex and truncate reads to 80 bp with process_radtags in parallel:     
`00-scripts/02_process_radtags_2_enzymes_parallel.sh 80 pstI mspI 8`    

Collect samples from multiple runs together and rename samples:    
`./00-scripts/03_rename_samples.sh`

Prepare the population map file:     
`./00-scripts/04_prepare_population_map.sh`


## 2. Align samples and quality filter on samples
### a. Align samples against the reference genome
Index the reference genome with bwa:    
`bwa index <ref_genome>`    

Point the alignment script to the reference genome (full path) using the GENOMEFOLDER and GENOME variables, and run:         
`00-scripts/bwa_mem_align_reads.sh 6`     

### b. Inspect alignment results
Compare per-sample reads and alignments, and per-sample reads and number of aligned scaffolds:      
`./../ms_scallop_popgen/01_scripts/assess_results.sh`    

Produces:      
```
04-all_samples/reads_per_sample_table.txt
04-all_samples/mappings_per_sample.txt
# A graph of number reads and aligned reads in the main directory
```

Optional other calculations:
```
# Total reads in all samples:     
awk '{ print $2 } ' 04-all_samples/reads_per_sample_table.txt | paste -sd+ - | bc
# Total reads before de-multiplexing (note: divide by 4 due to fastqc):   
for i in $(ls 02-raw/*.fastq.gz) ; do echo $i ; gunzip -c $i | wc -l ; done
```

### c. Filter out low reads/alignments individuals
Make directory to remove samples:    
`mkdir 04-all_samples/removed_samples`

1. Manually remove specific individuals with too few reads:     
`mv 04-all_samples/PIP_631.* 04-all_samples/removed_samples/`      
note: you can now go back and recalculate sample stats if you choose, but save the original outputs to avoid writing over.       

2. Manually remove these entries from the pop map

## 3. Genotype
### a. Run Stacks v.2.0 genotyper    
Update the number of cores, then run:      
`00_scripts/stacks2_gstacks_reference.sh`        

### b. Run populations to output VCF to identify outlier het samples
This step will be used to identify any problematic samples. Use single SNP per locus VCF:            
```
# Calculate inbreeding coefficient per sample
vcftools --het --vcf 05-stacks/populations.snps.vcf --out ./07-filtered_vcfs/samples

# Format data
awk '{ print $5, $1, $1 }' 07-filtered_vcfs/samples.het | cut -d "_" -f 1,2 > 07-filtered_vcfs/samples.het.data

# Plot
./00-scripts/utility_scripts/plot_heterozygozity.R 07-filtered_vcfs/samples.het.data

# Observe output and view any problematic individuals
 
```

If you are concerned about any individuals, remove them from the population map, then re-run gstacks. In our case, we will remove individual VIU002, due to a beyond high level of heterozygosity.           

### c. Re-run Stacks v.2.0 genotyper
`00_scripts/stacks2_gstacks_reference.sh`        


### d. Filter using the populations module      
Update the populations module script `00-scripts/stacks2_populations_reference.sh` as follows, run once for single-snp and once for microhaps:          
```
# single SNP per locus
populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 3 -r 0.7 \
    --ordered-export --fasta-loci --vcf \
    --min-maf 0.01 --hwe --plink --write-single-snp

# Save out all populations files appropriately
mkdir 05-stacks/popn_out_single_snp
mv 05-stacks/populations.* 05-stacks/popn_out_single_snp/

# microhaplotypes
populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 3 -r 0.7 \
    --ordered-export --fasta-loci --vcf \
    --min-maf 0.01 --hwe --plink \
    --radpainter

# Save out all populations files appropriately
mkdir 05-stacks/popn_out_microhaps
mv 05-stacks/populations.* 05-stacks/popn_out_microhaps/

``` 

### c. Convert output plink files
Convert output plink files to useable format for adegenet:        
`plink --ped 05-stacks/populations.plink.ped --map 05-stacks/populations.plink.map --maf 0.01 --recode A --allow-extra-chr --out 05-stacks/populations_single_snp`       
Note: --recode creates a new text fileset, and the A modifier causes additive (0/1/2) format to be generated, which is useful for input into R.           
Will output `05-stacks/populations_single_snp.raw`, which is the input file in the R code for the analysis.      


## 4. Genetic analysis in R
Open and run interactively:       
`ms_scallop_popgen/01_scripts/pyes_popgen_analysis.R`        
note: this will assume you have the `*.raw` file in your folder `05-stacks`, and will output results into `stacks_workflow/12-results`.     

Initial steps will allow you to view the minor allele frequency, and to convert the genlight file into a genind. After the data is in genind format, proceed to `simple_pop_stats`.         


## 5. Relatedness
Use the output radpainter file, with fineRADstructure.      
If stacks workflow is at same level as current repo:     
`cp ../stacks_workflow/05-stacks/populations.haps.radpainter ./02_inputs/`            
note: future version should create its own folder, as output of initial command outputs to input folder too.       
note: future versions should ensure populations output has RAD loci ordered according to genomic coord

fineRADstructure process:        
```
# Infer the co-ancestry matrix from RAD-seq data
RADpainter paint 02_inputs/populations.haps.radpainter

# Infer recent shared ancestry using haplotype linkage information and coalescence; assign indiv to popn
finestructure -x 100000 -y 100000 -z 1000 02_inputs/populations.haps_chunks.out 02_inputs/populations.haps_chunks.mcmc.xml

# Build tree
finestructure -m T -x 10000 02_inputs/populations.haps_chunks.out 02_inputs/populations.haps_chunks.mcmc.xml 02_inputs/populations.haps_chunks.mcmcTree.xml

# Copy R scripts from your programs folder (included with installation of fineRADstructure)
cp ~/programs/fineRADstructure/fineRADstructurePlot.R ./01_scripts/
cp ~/programs/fineRADstructure/FinestructureLibrary.R ./01_scripts/

# Edit the required sections of fineRADstructurePlot.R and run interactively
# The output figures will be in 04_fineRADstructure/ 

```








