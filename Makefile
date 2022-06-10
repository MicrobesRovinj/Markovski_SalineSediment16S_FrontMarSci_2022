export MOTHUR = code/mothur/mothur
export RAW = data/raw
export MOTH = data/mothur
export REFS = data/references
BASIC_STEM = data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster
export FIGS = results/figures
export NUM = results/numerical
export FINAL = submission

# Obtaining the Linux version of mothur from the mothur GitHub repository
$(MOTHUR) :
	wget --no-check-certificate https://github.com/mothur/mothur/releases/download/v1.45.2/Mothur.Centos7.zip
	unzip Mothur.Centos7.zip
	mv mothur code/
	rm Mothur.Centos7.zip
	rm -rf __MACOSX
	# The binary distribution of vsearch (part of mothur) is not working on cluster Isabella,
	# so vsearch was build from source.
	wget https://github.com/torognes/vsearch/archive/v2.17.0.tar.gz
	tar xzf v2.17.0.tar.gz
	mkdir code/vsearch/
	cd vsearch-2.17.0/; ./autogen.sh
	cd vsearch-2.17.0/; ./configure --prefix=$$(dirname $$PWD)/code/vsearch/
	make -C vsearch-2.17.0/
	make install -C vsearch-2.17.0/
	mv code/vsearch/bin/vsearch code/mothur/
	rm -rf vsearch-2.17.0/ v2.17.0.tar.gz code/vsearch/

#########################################################################################
#
# Part 1: Creating the reference files
#
# 	We will need several reference files to complete the analyses including the
# SILVA reference alignment and taxonomy. As we are analyzing both Bacteria and
# Archaea we need to optimize the procedure described on the mothur blog
# (https://mothur.org/blog/2021/SILVA-v138_1-reference-files/).
#
#########################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. We will use the
# Release 138.1. The curation of the reference files to make them compatible with
# mothur is described at https://mothur.org/blog/2021/SILVA-v138_1-reference-files/.
# As we are using primers from the Earth Microbiome Project that are targeting
# both Bacteria and Archaea (http://www.earthmicrobiome.org/protocols-and-standards/16s/)
# we need to modify the procedure described at
# https://mothur.org/blog/2021/SILVA-v138_1-reference-files/
# as this approach is removing shorter archeal sequences.
#
# The SILVA Release 138.1 was downloaded from
# https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/ARB_files/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
# opened in ARB and exported to silva.full_v138_1.fasta file as described at
# https://mothur.org/blog/2021/SILVA-v138_1-reference-files/ uder the
# section Getting the data in and out of the ARB database. A total of 446,881
# sequences were exported.

# Defining environmental variables
export REFS_START_SCREEN = 11894
export REFS_END_SCREEN = 25319
export REFS_MAXAMBIG = 5
# Defining the SILVA taxa mapping file
export SILVA_MAP_ADDRESS = https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz
export SILVA_MAP = tax_slv_ssu_138.1.txt.gz

# Creating the reference files
$(REFS)/silva%tax\
$(REFS)/silva.pcr%align\
$(REFS)/silva.pcr.unique%align : code/get_references.sh\
                                 $(MOTHUR)\
                                 ~/references/silva.full_v138_1/silva.full_v138_1.fasta
	# Copying the silva.full_v138_1.fasta file
	cp ~/references/silva.full_v138_1/silva.full_v138_1.fasta $(REFS)/silva.full.fasta
	# Running the script
	bash code/get_references.sh

#########################################################################################
#
# Part 2: Running data through mothur and calculating the sequencing error
#
# 	Processing fastq files through mothur and generating files that will be used
# in the overall analysis.
#
#########################################################################################

# Generating raw.files for mothur's function make.contigs
$(RAW)/raw.files : $(RAW)/metadata.csv
	cut -f 1,2,3 data/raw/metadata.csv | tail -n +2 > $(RAW)/raw.files

# Downloading project fastq.gz files from the European Nucleotide Archive (ENA)
$(RAW)/18118-*%fastq : $(RAW)/NC_*.fastq
	wget "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB48899&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true&limit=0" -O $(RAW)deposited_list.txt
	cut -f $$(head -1 $(RAW)deposited_list.txt | tr "\t" "\n" | cat -n | grep "submitted_ftp" | cut -f 1) $(RAW)deposited_list.txt > $(RAW)fastq.gz.txt
	tail -n +2 $(RAW)fastq.gz.txt | tr ";" "\n" > $(RAW)fastq.gz_list.txt
	sed -e "s/^/ftp:\/\//" -i $(RAW)fastq.gz_list.txt
	wget -i $(RAW)fastq.gz_list.txt -P $(RAW)
	wget "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR399/ERR3991476/18118-0422_R1_all.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR399/ERR3991476/18118-0422_R2_all.fastq.gz" -P $(RAW)
	gzip -d $(RAW)/*.gz

# The primer.oligos file containing the sequences of gene specific primers must be
# specified. Also, the atcc_v4.fasta file containing sequences of the mock community
# members must be added.

# Defining environmental variables
export MAXAMBIG = 0
export MINLENGTH = 225
export MAXLENGTH = 275
export ALIGNREF = $(REFS)/silva.pcr.unique.align
export STARTSCREEN = 1968
export ENDSCREEN = 11550
export MAXHOMOP = 19
export CLASSIFYREF = $(REFS)/silva.pcr.align
export CLASSIFYTAX = $(REFS)/silva.tax
export CONTAMINENTS = Chloroplast-Mitochondria-unknown-Eukaryota
export REMOVEGROUPS = NC_3-NC_4-ATCC_2-ATCC_3-ATCC_4-ATCC_5
export MOCKGROUPS = ATCC_2-ATCC_3-ATCC_4-ATCC_5
export MOCKREF = $(REFS)/atcc_v4.fasta

# Running fastq files through mothur, generating the summary.txt file to check if
# all went alright throughout the analysis and calculating the sequencing error.
$(BASIC_STEM).pick.silva.wang.tax%summary\
$(BASIC_STEM).pick.pick.pick.opti_mcc%shared\
$(MOTH)/summary%txt\
$(BASIC_STEM).pick.pick.pick.error%summary : code/get_processed.batch\
                                             code/get_summary.batch\
                                             $(MOTHUR)\
                                             $(RAW)/primer.oligos\
                                             $(RAW)/raw.files\
                                             $(RAW)/NC_*.fastq\
                                             $(RAW)/18118-*.fastq\
                                             $(REFS)/silva.pcr.align\
                                             $(REFS)/silva.pcr.unique.align\
                                             $(REFS)/silva.tax\
                                             $(REFS)/atcc_v4.fasta
	$(MOTHUR) code/get_processed.batch
	rm $(MOTH)/*.map
	$(MOTHUR) code/get_summary.batch

#########################################################################################
#
# Part 3: Figures and tables generation
#
#       Run scripts to generate figures and tables.
#
#########################################################################################

# Constructig the rarefaction plot
$(FIGS)/rarefaction_a%jpg\
$(FIGS)/rarefaction_b%jpg : $(BASIC_STEM).pick.pick.pick.opti_mcc.shared\
                            code/get_rarefaction.batch\
                            $(MOTHUR)\
                            code/plot_rarefaction.R\
                            $(RAW)/metadata.csv
	$(MOTHUR) code/get_rarefaction.batch
	R -e "source('code/plot_rarefaction.R')"

# Generating randomly rarefied community data
$(NUM)/rarefied.Rdata : $(BASIC_STEM).pick.pick.pick.opti_mcc.shared\
                        code/get_rarefied.R
	R -e "source('code/get_rarefied.R')"

# Calculating richness estimators and diversity indices
$(NUM)/estimators_indices_metadata.Rdata : $(NUM)/rarefied.Rdata\
                                           code/get_estimators_indices_metadata.R\
                                           $(RAW)/metadata.csv
	R -e "source('code/get_estimators_indices_metadata.R')"

# Plotting richness and diversity calculators
$(FIGS)/calculators.jpg : $(NUM)/estimators_indices_metadata.Rdata\
                          code/plot_calculators.R
	R -e "source('code/plot_calculators.R')"

# Generating temporal plots of estimators and indices
$(FIGS)/estimators_months%jpg\
$(FIGS)/diversity_indices_month%jpg : $(NUM)/estimators_indices_metadata.Rdata\
                                      code/plot_estimators_per_month.R
	R -e "source('code/plot_estimators_per_month.R')"

# Constructing PCoA plots
$(FIGS)/pcoa_figure%jpg\
$(FIGS)/pcoa_figure_areas_layers%jpg : $(NUM)/rarefied.Rdata\
                                       $(RAW)/metadata.csv\
                                       code/plot_pcoa.R\
                                       code/plot_pcoa_areas_layers.R
	R -e "source('code/plot_pcoa.R')"
	R -e "source('code/plot_pcoa_areas_layers.R')"

# Plotting similarity coefficients
$(FIGS)/matrix.jpg : $(NUM)/rarefied.Rdata\
                     $(RAW)/metadata.csv\
                     code/plot_matrix.R
	R -e "source('code/plot_matrix.R')"

# Formating data used for community composition bar plots
$(NUM)/community.Rdata : $(BASIC_STEM).pick.silva.wang.tax.summary\
                         code/get_community.R
	R -e "source('code/get_community.R')"

# Generating community composition bar plots
$(FIGS)/community_barplot%jpg\
$(FIGS)/community_barplot_major_1%jpg\
$(FIGS)/community_barplot_major_2%jpg\
$(FIGS)/community_barplot_month%jpg : $(NUM)/community.Rdata\
                                      $(RAW)/metadata.csv\
                                      $(RAW)/group_colors.csv\
                                      code/plot_community_barplot.R\
                                      code/plot_community_barplot_major.R\
                                      code/plot_community_barplot_month.R
	R -e "source('code/plot_community_barplot.R')"
	R -e "source('code/plot_community_barplot_major.R')"
	R -e "source('code/plot_community_barplot_month.R')"

#########################################################################################
#
# Part 4: Combining all together
#
# 	Render the manuscript and supplementary information.
#
#########################################################################################

$(FINAL)/manuscript.pdf\
$(FINAL)/supplementary.pdf : $(BASIC_STEM).pick.pick.pick.error.summary\
                             $(BASIC_STEM).pick.silva.wang.tax.summary\
                             $(RAW)/metadata.csv\
                             $(BASIC_STEM).pick.pick.pick.opti_mcc.shared\
                             $(NUM)/rarefied.Rdata\
                             $(NUM)/estimators_indices_metadata.Rdata\
                             $(NUM)/community.Rdata\
                             $(FIGS)/map.jpg\
                             $(FIGS)/calculators.jpg\
                             $(FIGS)/pcoa_figure.jpg\
                             $(FIGS)/pcoa_figure_areas_layers.jpg\
                             $(FIGS)/community_barplot.jpg\
                             $(FIGS)/community_barplot_major_1.jpg\
                             $(FIGS)/community_barplot_major_2.jpg\
                             $(FIGS)/rarefaction_a.jpg\
                             $(FIGS)/rarefaction_b.jpg\
                             $(FIGS)/estimators_months.jpg\
                             $(FIGS)/diversity_indices_month.jpg\
                             $(FIGS)/matrix.jpg\
                             $(FIGS)/community_barplot_month.jpg\
                             $(FINAL)/citation_style.csl\
                             $(FINAL)/header.tex\
                             $(FINAL)/header_supplementary.tex\
                             $(FINAL)/references.bib\
                             $(FINAL)/manuscript.Rmd\
                             $(FINAL)/supplementary.Rmd
	R -e 'render("$(FINAL)/supplementary.Rmd", clean=FALSE)'
	rm $(FINAL)/supplementary.knit.md
	R -e 'render("$(FINAL)/manuscript.Rmd", clean=FALSE)'
	rm $(FINAL)/manuscript.knit.md
	rm $(FINAL)/*.log $(FINAL)/*.out

# Cleaning
.PHONY: clean
clean :
	rm -f my_job.sge.* || true
	rm -f $(REFS)/tax* || true
	rm -f $(REFS)/silva* || true
	rm -f $(MOTH)/raw.* || true
	rm -f $(MOTH)/current_files.summary || true
	rm -f $(MOTH)/summary.txt || true
	rm -f $(RAW)/18118-*.fastq || true
	rm -f $(RAW)/raw.files || true
	rm -rf code/mothur/ || true
	rm -f $(FIGS)/*.jpg || true
	rm -f $(NUM)*.Rdata || true
	rm -f $(NUM)*.tex || true
	rm -f mothur*logfile || true
	rm -f $(FINAL)/manuscript.pdf || true
	rm -f $(FINAL)/supplementary.pdf || true
