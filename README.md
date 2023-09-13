Scripts for the "Ecological turnover of strain-level diversity in the honeybee gut microbiota between nurses and foragers" manuscript.
=======

This repository contains the necessary scripts and metadata required to reproduce the analyses and figures from the "Ecological turnover of strain-level diversity in the honeybee gut microbiota between nurses and foragers" manuscript. They use metagenomic data from honeybee gut DNA sequencing with Illumina.

# Acknowledgement

Most of the scripts and analyses from this manuscript have been adapted, developed or taken directly from the work of Kirsten Ellegaard, in particular the following publications:

> Kirsten Maren Ellegaard & Philipp Engel. **Genomic diversity landscape of the honey bee gut microbiota**; _Nature Communications_ **10**, Article number: 446 (2019).
> PMID: 30683856;
> doi:[10.1038/s41467-019-08303-0](https://www.nature.com/articles/s41467-019-08303-0)

> Kirsten Maren Ellegaard, Shota Suenami, Ryo Miyasaki, Philipp Engel. **Vast differences in strain-level diversity in the gut microbiota of two closely related honey bee species**; _Current Biology_ **10**, Epub 2020 Jun 11.
> PMID: 32531278;
> doi: [10.1016/j.cub.2020.04.070](https://www.cell.com/current-biology/fulltext/S0960-9822(20)30586-8)

If you are using the scripts found here, please cite them as well.

 
# Pre-requisites
----------

The scripts and analyses rely on a manually curated database of bacterial strain genomes built by Kirsten Ellegaard ([beebiome_db](https://drive.switch.ch/index.php/s/vHmM6aIIFyDGQwm)). These strains were mostly isolated from guts of *Apis mellifera*, *Apis cerana* and bumble bee species.

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
* (optional) [OrthoFinder](https://github.com/davidemms/OrthoFinder) (creation of orthologous groups of gene and orf catalogue) (version 2.3.5, should be added to the path) - results of this step are packaged with the beebiome database

A list of packages from the conda environments used for plotting and downstream analysis can be found in envs/python_env.yml and envs/r_env.yml and the files be used to create those environments.

# Running the pipeline
--------

## Data preparation: mapping and filtering 

First, download all he necessary data. These are the raw sequencing reads (which should be put in a `raw_reads` directory).

Then you can run `00_data_preparation.sh` script.
```./00_data_preparation.sh```

## Community composition analysis

### The reference database - beebiome_db

The reference database is to be downloaded as the (`beebiome_db`) directory. The reference database can be downloaded from [here](https://drive.switch.ch/index.php/s/vHmM6aIIFyDGQwm) **TODO: add zenodo link**. This link contains a zipped file. After unzipping it should contian the following files and directories:

* `honeybee_genome.fasta` : fasta file containing the host (_Apis mellifera_) genome sequence
* `beebiome_db` : fasta file of 198 concatenated genomes with one genome per entry (multi-line fasta) where the headers represent the genome identifier
* `beebiome_red_db` : fasta file of 39 species representative genomes with one genome per entry (multi-line fasta) where the headers represent the genome identifier to be used for analysis of intra-specific variation
* `fna_files` : directory containing genome sequence files and concatenated files where the concatenated files contain one fasta entry renamed to the genome identifier and all contigs concatenated into one entry
* `ffn_files` : directory containing one file per genome listing the nucleotide sequence of all the predicted genes
* `faa_files` : directory containing one file per genome listing the amino acid sequence of all the predicted genes
* `bed_files` : directory containing bed files where the location of each of the predicted genes are indicated based on their position in the concatenated genome file
* `single_ortho` : directory containing one file per phylotype listing all the single-copy orthogroups (OGs) identified by orthofinder where each line represents an OG id followed by a list of genes from each of the genomes of that phylotype that belong to that OG and the corresponding sequences of these genes can be found in the ffn file belonging to the respective genomecat fna
* `red_bed_files` : directory containing bed files for species representative genomes that only list the positions genes that belong to the core orthogroups of their phylotype

```./01_community_composition_analysis.sh```

## Strain-level analysis

Same here.  

```./02_strain-level_analysis.sh```

## Functional gene content analysis

And here too.  

```./03_functional_gene_content.sh```

