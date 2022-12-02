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
`plink` https://zzz.bwh.harvard.edu/plink/download.shtml      

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
`00_scripts/stacks2_populations_reference.sh`        

### b. Filter using the populations module      
Update to ensure output of plink files, then run:      
`00-scripts/stacks2_populations_reference.sh`      

