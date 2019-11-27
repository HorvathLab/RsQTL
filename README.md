# RsQTL: correlation of expressed SNVs with splicing patterns using RNA-sequencing data

This toolkit contains the required scripts to transform sequencing files into RsQTL input files and run the MatrixEQTL R package to identify significant variation-splicing relationships.


&nbsp;


## Getting Started

These instructions will get you a copy of the scripts up and running on your machine for development and testing purposes. See *Running the scripts* for notes on how to use the project on a live system. We have provided sample data that can be used to test the pipeline. It also serves as an example of the data format the pipeline expects.

This package was developed on R version 3.6.1 on macOS High Sierra.

### Prerequisites

* The following R packages installed on your machine:

  ```
  tidyverse
  MatrixEQTL
  data.table
  GenomicRanges
  factoextra
  ggpubr
  ```
  
  You may install the packages using the following commands:
  ```
  install.packages(c("tidyverse", "MatrixEQTL", "data.table", "factoextra", "ggpubr", "BiocManager"))
  BiocManager::install("GenomicRanges")
  ```

* Each of the following scripts copied to a working directory on your machine

	```
	build_junction_input_RsQTL.R
	build_VAF_matrix_RsQTL.R
	harmonize_matrices_RsQTL.R
	build_pca_covariate_matrix_RsQTL.R
	run_matrix_RsQTL.R
	annotate_cis_trans_RsQTL.R
	plot_RsQTL.R
	```
	You can obtain the full toolkit [here.](https://github.com/HorvathLab/RsQTL/archive/master.zip)
* Output *.csv* files from our ReadCounts tool (https://github.com/HorvathLab/NGS/tree/master/readCounts) containing the read counts extracted per SNV for each sample

* A *.counts* matrix from *LeafCutter*

## Running the scripts

The scripts are designed to be run from the *Unix command line* (Terminal on macOS) from the root directory (by default RsQTL-master if the toolkit was downloaded from the link above). Make sure to *cd* to this directory before beginning.

***

&nbsp;

### build_splicing_matrix_RsQTL.R

Transforms the raw LeafCutter output into a matrix suitable for RsQTL analysis

#### Input
* The path to a raw LeafCutter output matrix (-c)
* The desired prefix of the output junction matrices (-o)
* *OPTIONAL:* the minimum number of reads spanning an intron junction for inclusion (--read-thresh, default = 30)
* *OPTIONAL:* maximum fraction of samples allowed to have a cluster read ratio of NA for a given intron (--na-thresh, default = 0.8)
* *OPTIONAL:* maximum fraction of samples allowed to have a cluster read ratio of 1 for a given intron  (--one-thresh, default = 0.8)
* *OPTIONAL:* maximum fraction of samples allowed to have a cluster read ratio of 0 for a given intron  (--zero-thresh, default = 0.8)
* *OPTIONAL:* logical specifying whether to remove introns above the X%ile of NA values (--excl-na, default = F)
* *OPTIONAL:* %ile threshold above which to exclude introns (--perc-thresh, default = 0.9)
* *OPTIONAL:* logical specifying whether to exclude introns which represent the only member of a cluster (--excl-sing, default = F)



#### Output
* One file (in the output directory) with the cluster read ratios for each intron in each sample
* One file (in the output directory) with the genomic locations for each intron

#### Sample Command
```
Rscript build_splicing_matrix_RsQTL.R -c data/sample_junction_matrix.counts -o RsQTL_test
```

&nbsp;

***

&nbsp;

### build\_VAF_matrix_RsQTL.R

Transforms the read counts into a variant fraction matrix with information from all provided samples

#### Input
* A directory containing the *.csv* files from the output of Readcounts	(-r)
* The desired prefix of the output SNV matrix and SNV location files (-o)


#### Output
* One file (in the script’s directory) with the SNV locations for MatrixEQTL 
* One file (in the script’s directory) with the SNV variant allele fraction matrix for MatrixEQTL


#### Sample command
```
Rscript build_VAF_matrix_RsQTL.R -r data/ -o RsQTL_test
```
&nbsp;

***

&nbsp;

### build\_pca_covariate_matrix_RsQTL.R

Creates a covariate matrix of principal components from the VAF and junction matrices. Additional covariates can be included (see sample data for expected format).

#### Input
* The path to the VAF matrix created by build_VAF_matrix_RsQTL.R (-r)
* The path to the splicing matrix created by build_junction_matrix_RsQTL.R (-s)
* *OPTIONAL:* The path to any additional covariates (-c)
* *OPTIONAL:* The number of PCs to include as covariates (-n, default = 10)
  * *NOTE:* You may want to run this script once with the default value then adjust the number of PCs used based on the variance plots. Read more about selecting an appropriate number of principal components [here.](http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/)
* The desired prefix of the output covariate matrix (-o)


#### Output
* A covariate matrix containing the top *n* VAF and cluster read ratio principal components as well as any additional supplied covariates
* Scree plots showing the percentage of variance explained by principal components of the VAF and splicing matrices


#### Sample command
```
Rscript build_pca_covariate_matrix_RsQTL.R -r output/RsQTL_test_VAF_matrix.txt -s output/RsQTL_test_splicing_matrix.txt -c data/additional_covariates_matrix.txt -n 10 -o RsQTL_test
```

&nbsp;

***

&nbsp;

### harmonize\_matrices_RsQTL.R

Harmonizes matrices so that all inputs for run_matrix_RsQTL.R contain the same samples (optional, but helps to avoid errors)

#### Input

* The path to the VAF matrix created by build_VAF_matrix_RsQTL.R
* The path to the splicing matrix created by build_splicing_matrix_RsQTL.R
* The path to the covariate matrix created by build_pca_covariate_matrix_RsQTL.R (if no covariate information is used, you will need to modify this script to remove the references to the covariate matrix)

#### Output
* Three matrices (in the output directory) corresponding to the three input files, but including only samples that were present in all three input files


#### Sample command
```
Rscript harmonize_matrices_RsQTL.R -r output/RsQTL_test_VAF_matrix.txt -s output/RsQTL_test_splicing_matrix.txt -c output/RsQTL_test_covariate_pca_matrix.txt
```

&nbsp;

***

&nbsp;

### run\_matrix_RsQTL.R

*This script is based off of the sample code from Shabalin, et al (2012)*

Runs the RsQTL analysis using MatrixEQTL

#### Input

* Names of the SNV matrix (-s), SNV location (-sl), splicing intron matrix (-i), and splicing intron location (-il) files
* *OPTIONAL:* Name of the covariates matrix (-c); we include an example in the "data" folder
* Logical (T or F, -ct) specifying whether to split the output into *cis* and *trans*
* Name of the output file for the qq plot (-qq)
* Prefix of the output file(s) (-o)
* P-value thresholds for the *cis* (-pcis) and *trans* (-ptr) output files or the unified output file (-p) depending on the logical specified above


#### Output
* One file (in the output directory) with the *cis* RsQTLs
* One file (in the output directory) with the *trans* RsQTLs
OR
* One file (in the output directory) with all of the unified RsQTLs depending on the logical specified above


#### Sample commands

Splitting *cis* and *trans*
```
Rscript run_matrix_RsQTL.R -s output/RsQTL_test_VAF_matrix_harmonized.txt -sl output/RsQTL_test_VAF-loc_matrix.txt -i output/RsQTL_test_splicing_matrix_harmonized.txt -il output/RsQTL_test_splicing-loc_matrix.txt -c output/RsQTL_test_pca_covariate_matrix_harmonized.txt -ct T -qq RsQTL_test_qqplot -pcis 0.1 -ptr 0.1 -o RsQTL_test
```

Unified *cis* and *trans*
```
Rscript run_matrix_RsQTL.R -s output/RsQTL_test_VAF_matrix_harmonized.txt -sl output/RsQTL_test_VAF-loc_matrix.txt -i output/RsQTL_test_splicing_matrix_harmonized.txt -il output/RsQTL_test_splicing-loc_matrix.txt -c output/RsQTL_test_pca_covariate_matrix_harmonized.txt -ct F -qq RsQTL_test_qqplot -p 0.1 -o RsQTL_test
```
&nbsp;

### annotate\_cis_trans_RsQTL.R

Annotates the output of RsQTL as cis/trans based on whether the SNV and paired intron lie within the same gene

#### Input

* The path to a RsQTL results file (-r)
* The path to a list of gene locations (-g, see sample data for a suggested file)
* The desired prefix of the output annotated results file (-o)

#### Output
* One file (in the output directory) with the RsQTLs annotated as *cis* or *trans*


#### Sample command
```
Rscript annotate_cis_trans_RsQTL.R -r output/RsQTL_test_all_RsQTLs.txt -g data/gene_locations_hg38.txt -o RsQTL_test
```
&nbsp;

### plot_RsQTL.R

Plots either the top n most significant RsQTLs in the input file or a specific SNV-intron pair

#### Input

* The path to the input VAF matrix (-r)
* The path to the input splicing matrix (-s)
* Plotting mode (-m, "bulk" or "single")
* (If plotting in bulk mode) the path to a results file (-res)
* (If plotting in bulk mode) the number of correlations to plot (-n)
* (If plotting in single mode) the SNV to plot (-snv)
* (If plotting in single mode) the intron to plot (-intron)
* The desired prefix of the output annotated results file (-o)

#### Output
* One file (in the output directory) containing either the bulk or individual plot


#### Sample command

Bulk mode
```
Rscript plot_RsQTL.R -r output/RsQTL_test_VAF_matrix_harmonized.txt -s output/RsQTL_test_splicing_matrix_harmonized.txt -res output/RsQTL_test_RsQTLs_cistrans_ann.txt -m bulk -n 200 -o RsQTL_test
```

Single mode
```
Rscript plot_RsQTL.R -r output/RsQTL_test_VAF_matrix_harmonized.txt -s output/RsQTL_test_splicing_matrix_harmonized.txt -res output/RsQTL_test_RsQTLs_cistrans_ann.txt -m single -snv "1:952657_T>C" -intron "chr12:56160320_56161387" -o RsQTL_test
```

&nbsp;

## Authors and Acknowledgements

*Justin Sein and Liam F. Spurr*, Pavlos Bousounis, Nawaf Alomran, Dacian Reece-Stremtan, Prashant N M, Hongyu Liu, and Anelia Horvath

We would like to thank the Matrix EQTL team (Shabalin, et al. 2012) for their sample code and R package upon which *run\_matrix_RsQTL.R* is based.
