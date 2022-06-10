## Compositional stability of sediment microbial communities during a seagrass meadow decline
This is the repository for the manuscript "Compositional stability of sediment microbial communities during a seagrass meadow decline" written by Marsej Markovski, Mirjana Najdek, Gerhard J. Herndl and Marino Korlević. The raw sequencing data with the exception of negative controls have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession numbers SAMEA11293274 - SAMEA11293412 and SAMEA6648825. Negative control samples are part of this repository and are located in data/raw/. To be able to reproduce the results the mothur compatible SILVA reference file (Release 138.1, available under a [CC-BY license](https://www.arb-silva.de/silva-license-information/)) must be created according to instruction given on [the mothur blog](https://mothur.org/blog/2021/SILVA-v138_1-reference-files/) and in the Makefile. This README file contains an overview of the repository structure, information on software dependencies and instructions how to reproduce and rerun the analysis.

### Overview

	project
	|- README                       # the top level description of content (this doc)
	|- LICENSE                      # the license for this project
	|
	|- submission/                  # files necessary for manuscript or supplementary information rendering, e.g. executable Rmarkdown
	| |- manuscript.Rmd             # executable Rmarkdown for the manuscript of this study
	| |- manuscript.md              # Markdown (GitHub) version of the manuscript.Rmd file
	| |- manuscript.tex             # TeX version of manuscript.Rmd file
	| |- manuscript.pdf             # PDF version of manuscript.Rmd file
	| |- manuscript.aux             # auxiliary file of the manuscript.tex file, used for cross-referencing
	| |- header.tex                 # LaTeX header file to format the PDF version of manuscript
	| |- supplementary.Rmd          # executable Rmarkdown for the supplementary information of this study
	| |- supplementary.md           # Markdown (GitHub) version of the supplementary.Rmd file
	| |- supplementary.tex          # TeX version of supplementary.Rmd file
	| |- supplementary.pdf          # PDF version of supplementary.Rmd file
	| |- supplementary.aux          # auxiliary file of the supplementary.tex file, used for cross-referencing
	| |- header_supplementary.tex   # LaTeX header file to format the PDF version of supplementary information
	| |- references.bib             # BibTeX formatted references
	| +- citation_style.csl         # csl file to format references
	|
	|- data/                        # raw, reference, and primary data are not changed once created
	| |- references/                # reference files to be used in analysis
	| |- raw/                       # raw data, will not be altered
	| +- mothur/                    # mothur processed data
	|
	|- code/                        # any programmatic code
	|
	|- results/                     # all output from workflows and analyses
	| |- figures/                   # manuscript or supplementary information figures
	| +- numerical/                 # results of the statistics or other numerical results for manuscript or supplementary information
	|
	|-. gitignore                   # gitignore file for this study
	|-. Rprofile                    # Rprofile file containing information on which R libraries to load,
	|                               # rendering options for knitr and Rmarkdown, functions etc.
	+- Makefile                     # executable Makefile for this study

### How to regenerate this repository

#### Dependencies
* GNU Bash (v. 4.2.46(2))
* GNU Make (v. 3.82); should be located in user's PATH
* mothur (v. 1.45.2)
* R (v. 3.6.3); should be located in user's PATH
* R packages:
  * `stats (v. 3.6.3)`
  * `knitr (v. 1.33)`
  * `rmarkdown (v. 2.11)`
  * `tinytex (v. 0.31)`
  * `tidyverse (v. 1.3.1)`
  * `vegan (v. 2.5.7)`
  * `RColorBrewer (v. 1.1.2)`
  * `kableExtra (v. 1.3.4)`
  * `cowplot (v. 1.1.1)`
  * `lemon (v. 0.4.5)`

#### Running analysis
Before running the analysis be sure to generate the mothur compatible SILVA reference file and indicate in the Makefile its location. The manuscript and supplementary information can be generated on a Linux computer by running the following commands:
```
git clone https://github.com/MicrobesRovinj/Markovski_SalineSediment16S_x_2022.git
cd Markovski_SalineSediment16S_x_2021/
make submission/manuscript.pdf
```
If something goes wrong and the analysis needs to be restarted run the following command from the project home directory before rerunning the analysis:
```
make clean
```
