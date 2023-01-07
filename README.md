
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CVrefDB R package

## Overview

Cross Validation (CV) of DNA barcodes reference databases (refDB) using
blast.

The main aim of this package is to repeatedly blast sequences of a DNA
barcode reference database against itself to evaluate the quality of the
taxonomic assignments at various taxonomic levels and using different
strategies.

There are two main functions :

- `CV_blastn()` : This function will repeatedly extract a random sample
  of sequences from a DNA barcode reference database in fasta format,
  then blast these query sequences against the remainder of the database
  (k-fold cross-validation) or the whole database (leaked
  cross-validation) and return the output of the top blast hits with the
  true and predicted taxonomies and blastn statistics (bit score,
  identity,…).
- `assign_taxonomy()` : This function will assign a taxonomy to
  sequences based on their top blastn hits. It can use different methods
  for the assignment based on the best bit score, the consensus score
  and it allows optionally to specify the minimum identity, length and
  E-value to take into consideration for each taxonomic level. The
  output contains a taxonomic assignment for each taxonomic level
  (species, genus, family,…) and identity and consensus scores that can
  be used to decide to dismiss untrustworthy assignments

The vignette provides a simple example on how these functions interact
with a few smaller functions useful to exploit and visualize the
outputs. You can access the vignette after package installation with the
following R command : `browseVignettes("CVrefDB")` or dirrectly on
[github](https://raw.githubusercontent.com/GillesSanMartin/CVrefDB/master/vignettes/CVrefDB.pdf).

## Installation

**1. Install blast**

`blastn` must be installed locally on your machine through you favorite
package manager or by downloading installers from [the NCBI
website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).  
Under Windows you need to install for example the latest win64.exe file
available on the [download
page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).  
Under debian based GNU/Linux you can for example use
`sudo apt install ncbi-blast+`.

To check that the installation works fine and that blast is visible from
R you can type the folowing command in the R console :
`system("blastn -version")` which should print your `blatn` version (and
no error message).

**2. Install `Biostrings` from (Bioconductor)**

This package is used to read/write fasta files. It is not on CRAN but on
another package repository dedicated to bioinformatics and its
installation is slightly different than CRAN packages :

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("Biostrings")

**3. Install `CVrefDB`**

This package is a companion to a publication and is not intended to be
published on the CRAN.

To install the latest version form github :

    if (!require("devtools", quietly = TRUE))
        install.packages("devtools")

    devtools::install_github("GillesSanMartin/CVrefDB")

To install it from a bundled `tar.gz` file in the R console :

    # install first dependencies from CRAN
    install.packages(c("data.tables", "dplyr", "tidyr", "ggplot2"))
    # then install from source
    install.packages("/<insert your path here>/CVrefDB_0.0.1.tar.gz", 
                     type = "source", repos = NULL)

## Usage

Look at the [pdf
vignette](https://raw.githubusercontent.com/GillesSanMartin/CVrefDB/master/vignettes/CVrefDB.pdf)
for more details, graphs, …

Type `help(package = "CVrefDB")` in the R console to get the package
help index. Detailed examples are provided for each function.

``` r
library(CVrefDB)
library(Biostrings) # from Bioconductor

# Change printing options for tibbles
options(tibble.print_max = 100, tibble.print_min = 30, tibble.width = Inf)

# path to example files from the package
fasta_path <- system.file("extdata/ITS2_Rosales_Restricted.fasta",
                          package = "CVrefDB")
taxonomy_path <- system.file("extdata/ITS2_Rosales_Restricted.tsv",
                             package = "CVrefDB")

# Read the fasta and tsv files
fasta <- Biostrings::readDNAStringSet(fasta_path)
reftaxo <- read.table(taxonomy_path, sep = "\t", header = TRUE)

# inspect the content of these files :
fasta
head(reftaxo)

# 10 fold cross validation
output_10FoldCV <- CV_blastn(fasta_db = fasta, taxo = reftaxo, seed = 12, verbose = TRUE)
head(output_10FoldCV)

# assign taxonomy
assigned_long <- assign_taxonomy(output_10FoldCV, taxo = reftaxo,  Order = NA )
head(assigned_long)

# Reorganize the output in wide format
assigned_wide <- pivot_assign_taxonomy(assigned_long)
head(assigned_wide)

# compute the % of correct assignments for each family, at family, genus and species level
scores <- 
score_per_taxon(assigned_long[assigned_long$Method == "TopHitPlus",], 
                 grouping_tax_level = "Family", 
                 predicted_NA_wrong = TRUE)
scores[order(scores$Tax_level, -scores$Pct),] 
```
