Scripts for the "Turnover of strain-level diversity modulates functional traits in the honeybee gut microbiome between nurses and foragers" manuscript available as pre-print [here](https://doi.org/10.1101/2022.12.29.522137).

# Pre-requisites

The scripts and analyses rely on a manually curated database of bacterial strain genomes built by Kirsten Ellegaard ([beebiome_db](10.5281/zenodo.10182034)). These strains were mostly isolated from guts of *Apis mellifera*, *Apis cerana* and bumble bee species.

In terms of software, the scripts require:

* Python 3 (version 3.6 or higher)
* Bash
* R  (version 4.2.1, including packages: "segmented", "plyr", "ape" (version 5.6-2), "data.table" (version 1.14.2), "qvalue" (version 2.28.0), "PERMANOVA" (version 0.2.0), "vegan" (2.6-2), "ggplot2" (version 3.3.6), "reshape" (version 0.8.9), "gridExtra" (version 2.3), "ggrepel" (version 0.9.1), "RColorBrewer" (version 1.1-3) and "scales" (version 1.2.1))
* [samtools](http://www.htslib.org) (version  1.15-8-gbdc5bb8, should be added to the path)
* [bwa](https://github.com/lh3/bwa) (for the mapping) (version 0.7.15-r1142-dirty, should be added to the path)
* [SPAdes](https://github.com/ablab/spades) (for the assemblies) (version 3.10.1, should be added to the path)
* [freebayes](https://github.com/freebayes/freebayes) (for the detection of single nucleotide polymorphisms) (version v1.3.1-19-g54bf409, should be added to the path)
* [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (quality check of the input data) (should be added to the path)
* [Trimmomatic](https://github.com/usadellab/Trimmomatic) (trimming of the reads and removal of the sequencing adapters) (version 0.35, should be added to the path)
* [prodigal](https://github.com/hyattpd/Prodigal) (detection of ORFs) (version 2.6.3, should be added to the path)
* [eggnog](https://github.com/eggnogdb/eggnog-mapper) (annotation of the gene and ORF catalogue) (version 1.0.3-33-g70ff1ab, should be added to the path)
* [OrthoFinder](https://github.com/davidemms/OrthoFinder) (creation of orthologous groups of gene and orf catalogue) (version 2.3.5, should be added to the path)database

A list of packages from the conda environments used for plotting and downstream analysis can be found in `envs/python_env.yml` and `envs/r_env.yml` and the files be used to create those environments.

--------

# Running the pipeline

## Data preparation: mapping and filtering 

First, download all the necessary data. These are the raw sequencing reads (which should be put in a `raw_reads` directory).

Then you can run the script `00_data_preparation.sh`.

The script requires the file `NexteraPE-PE.fa`, containing the Nextera PE adapter sequences obtained from the trimmomatic github [repo](https://github.com/timflutre/trimmomatic/blob/master/adapters/NexteraPE-PE.fa). This script will:

1. make a `QC` directory removing any pre-existing directory by that name
2. make a `filtered_reads` directory removing any pre-existing directory by that name
3. for each sample in `raw_reads` directory the expected format is `{sample_name}_L{sequencing_lane_number}_R{orientation_1_or_2}_{unique_tag}.fastq.gz`
4. runs [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the raw reads
5. trims and filters the raw reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic)
6. runs [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the filtered samples
7. merges the filtered reads


## Community composition analysis

### The reference databases

The reference database is to be downloaded as the (`beebiome_db`) directory. The reference database can be downloaded from [here](10.5281/zenodo.10182034) **TODO: add zenodo link**. This link contains a zipped file. After unzipping, it should contain the following files and directories:

* `honeybee_genome.fasta` : fasta file containing the host (_Apis mellifera_) genome sequence
* `beebiome_db` : fasta file of 198 concatenated genomes with one genome per entry (multi-line fasta) where the headers represent the genome identifier
* `beebiome_red_db` : fasta file of 39 species representative genomes with one genome per entry (multi-line fasta) where the headers represent the genome identifier to be used for the analysis of intra-specific variation
* `fna_files` : directory containing genome sequence files and concatenated files where the concatenated files contain one fasta entry renamed to the genome identifier and all contigs concatenated into one entry
* `ffn_files` : directory containing one file per genome listing the nucleotide sequence of all the predicted genes
* `faa_files` : directory containing one file per genome listing the amino acid sequence of all the predicted genes
* `bed_files` : directory containing bed files where the location of each of the predicted genes are indicated based on their position in the concatenated genome file
* `single_ortho` : directory containing one file per phylotype listing all the single-copy orthogroups (OGs) identified by orthofinder where each line represents an OG id followed by a list of genes from each of the genomes of that phylotype that belong to that OG and the corresponding sequences of these genes can be found in the ffn file belonging to the respective genome
* `red_bed_files` : directory containing bed files for species representative genomes that only list the positions genes that belong to the core orthogroups of their phylotype

### Running the community composition analysis

The script `./01_community_composition_analysis.sh` will: 

1. map the filtered reads against the core genes of the beebiome bacterial database
   * It creates (after deleting any pre-existing directory of the same name) `mapping_full_db`
   * For each sample, it maps the reads to the beebiome_db and filters the resulting sam file using the script `filter_sam_aln_length.pl` to remove alignments with less than 50 matches. It also removes chimeric alignments and filters by edit distance \< 5 and the deleted all the intermediate files.
   * 
   * It extracted the reads that were unmapped and mapped with fewer than 50 matches and saves them as `mapping_full_db/${sample_tag}_vs_db_unmapped_R${1_or_2}.fastq`.
2. use the mapping data to infer terminus coverage of each of the bacterial speciess
   * It creates (after deleting any pre-existing directory of the same name) `community_composition`
   * It first makes symlinks of relavent scripts and files from the earlier directory inside this directory (./core_cov.pl, ./core_cov.R, ../mapping_full_db/\*_filt2x_sorted.bam ../beebiome_db/bed_files/\*.bed ../beebiome_db/single_ortho_files/*_single_ortho.txt)
   * For each phylotype for which a \*_single_ortho.txt file is available, it runs the ./core_cov.pl script which writes the output files `*_corecov.txt` per phylotype containing in each tab-separated column, the species, sample, gene, reference position (position of the corresponding gene of the same orthogroup in the reference genome of the species) and coverage (of the gene in this sample). It summarises the coverage of each gene per line.
   * Next, the core_cov.R script used the `*_corecov.txt` file to create the `*_corecov_coord.txt` file summarising in each line per species per sample, in tab-separated columns the species, sample, terminus coverage, peak-trough-ratio and difference in max and min ori coverage (represented by the left and right extremes in the plot of coverage vs ref position). This script also makes the plots in `*_corecov.pdf`.
   * The resulting `*_corecov_coord.txt` files are used by the script `figures_community_composition.R` for plotting the figures.

NOTE: this script is not optimized to run more easily. You can easily make it faster by parallelizing the mapping and filtering of mapping data operations. Most of the scripts here were written by Kirsten Ellegaard. Note that if you are using a subset of the samples or a subset of the database (or adding samples or strain genomes to the database), you will need to edit the perl "core_cov.pl" script in which the list of samples and genomes is hardcoded.

## Strain-level analysis

The script `./02_strain-level_analysis.sh` runs the first part of the analysis of strain-level variation.

1. maps reads to the reduced database `beebiome_red_db` and filters bam files as done before
2. makes softlinks to the required files and scripts and repeats community analysis using the scripts `core_cov_red.pl`, `core_cov.R`
3. using these results, filters the bed_red files to *filt.bed files which only include core genes that are also present in at least 10x coverage in each sample that the species is present at at least 20x coverage and the concatenates the them to `all_filt.bed` which demarcates the regions to be used for snv analysis
4. calculates the resulting core genome length using the script `calc_core_length.sh`
5. using freebayes, carries out variant calling across all samples with thresholds of 0.1 for min allele frequency and 10 for minimum coverage
6. using `vcfintersect` it subsets the resulting vcf files to inclucde only the regions specified in `all_filt.bed` and `vcfbreakmulti` to split lines specifying more that 2 alternate alleles in the same position, into separate line
7. further filters the vcf file to only keep samples with suffecient coverage of each species using `filter_vcf_samples.sh` and `filt_vcf_samples.pl` and the scripts also split the combined vcf files into multiple files each containing only the information relevant to each species

`parse_vcf_snvs.py` is then used to read the vcf files and pickle dictionaries containing information in each position as custom class (GenomePosition) objects. These dictionaries were then used to summarize results of various allele frequency thresholds and to generate the final set of intermediate files that were then parsed by `figures_strain-level_analyses.R` to plot figures

## Functional gene content analysis

The script `./03_functional_gene_content.sh` will:

1. combine reads unmapped to the host database with reads that were mapped to the beebiome_db and assemble them using [SPAdes](https://github.com/ablab/spades) with the meta option specified.
2. predict ORFs on the assembled contigs using [prodigal](https://github.com/hyattpd/Prodigal) and keep only ORFs that are at least 300 bp in length and marked as having a normal start and stop codon using the script `filt_ffn.pl`.
3. concatenates all the genes from genomes in the beebiome_db database and does a blast search of each gene predicted from the assemblies against it to assign each contig to a species if at least 80% of the ORFs on it were matched against genes of that species and at least 10% of the ORFs were matched against genes of that species and no other species and the results summarised in `functional_gene_content/all_filt_ORFs_list_assigned.txt`.
4. computes OGs using amino acid sequences of the ORFs and genes from genomes of the database with [OrthoFinder](https://github.com/davidemms/OrthoFinder)
5. maps reads against a catalog of nucleotide sequences of filtered orfs collected from across all samples and genes of genomes from the beebiome_db database and and calculates the coverage of each orf in each sample using [samtools](http://www.htslib.org) `coverage`
6. The catalog is annotated using [eggnog](https://github.com/eggnogdb/eggnog-mapper).
7. The rest of the analysis is carried out in the script `figures_functional_gene_content.R` which is run by the shell script as the last step and also makes the figures.

## Other analyses

### Core genome length summaries

The shell script can be used to calculate core lengths of a given set of bed files. The python script reads the lengths from the resulting file and then calculates the length of the total genomes from the fna files and then calculates the % that comprises the core genome before and after filtering.

### Gene length plots

The python script here was used to make coverage vs gene length plots of genes in orthogroups that were detected to be significant and non-significant to check that this was not baised by the length of the genes. A few example plots can be seen in the figures directory. If the mapping is done and the coverage files made by samtools from the bam file are available, the code can be used (after appropriately modifying the path) to make many summary plots to visualize coverage across genes. The intermediate files are large and can take a while to be made.

--------

# Acknowledgement

Most of the scripts and analyses from this manuscript have been adapted, developed or taken directly from the work of Kirsten Ellegaard, in particular the following publications:

> Kirsten Maren Ellegaard & Philipp Engel. **Genomic diversity landscape of the honey bee gut microbiota**; _Nature Communications_ **10**, Article number: 446 (2019).
> PMID: 30683856;
> doi:[10.1038/s41467-019-08303-0](https://www.nature.com/articles/s41467-019-08303-0)

> Kirsten Maren Ellegaard, Shota Suenami, Ryo Miyasaki, Philipp Engel. **Vast differences in strain-level diversity in the gut microbiota of two closely related honey bee species**; _Current Biology_ **10**, Epub 2020 Jun 11.
> PMID: 32531278;
> doi: [10.1016/j.cub.2020.04.070](https://www.cell.com/current-biology/fulltext/S0960-9822(20)30586-8)

If you are using the scripts found here, please cite them as well.

The scripts are authored and updated by multiple individuals. The perl scripts and some of the other scripts were authored by Kirsten Ellegaard and used in the analyses for the above mentioned publications. Other scripts were added and put together in this pipeline by Gilles Baud. Python scripts and updated R scripts for plotting were authored by Aiswarya Prasad.

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg