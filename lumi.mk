# Run Illumina microarray pipeline
#
# The primary targets in this file are:
#
# split_anno		Split an annotatio file from illumina into main and control
#                       probes.
# read_lumi		Read array data with Lumi
# preprocess_lumi	Run Lumi's preprocessing routines to background correct, 
#                       normalise, variance stabilise.
# extract_matrices	Extracts text-format matrices from serialised R objects.
# run_limma		Differential expression analysis with Limma.
# run_roast		Differential gene set analysis, also with Limma.
# make_shiny_object	Compile data for into an ExploratorySummarizedExperiment object as 
# 			defined by the shinyngs package.

# Default target
.PHONY: all
all: split_anno read_lumi preprocess_lumi extract_matrices run_limma run_roast make_shiny_object 

ROOT_DIR			:= /my/root/dir
ANALYSIS_DIR			:= $(ROOT_DIR)/analysis
OUT_DIR				:= $(ANALYSIS_DIR)/pipeline

# Input files

EXPERIMENT	 		:= <EXPERIMENT FILE HERE>
CONTRASTS			:= <CONTRAST FILE HERE>
BATCH_VARIABLE			:= <BATCH VARIABLE FROM EXPERIMENT FILE, or NA>
ILLUMINA_ANNO_FILE		:= <ILLUMINA ANNOTATION FILE HERE>
DATA_DIR			:= <ARRAY DATA LOCATION HERE>
SAMPLE_PROBE_PROFILE 		:= $(DATA_DIR)/<SAMPLE PROBE PROFILE HERE>
CONTROL_PROBE_PROFILE		:= $(DATA_DIR)/<CONTROL PROBE PROFILE HERE>
SAMPLES_TABLE			:= $(DATA_DIR)/<SAMPLES TABLE>
MSIGDB_DIR			:= <DIRECTORY WITH .gmt FILES FROM MSIGDB OR SIMILAR>
SPECIES				:= <DEFINE AN ENSEMBL-TYPE SPECIES, e.g. hsapiens>
THREADS				:= 8

# Illumina's probe annotation for control probes etc

ANNO_DIR			:= $(OUT_DIR)/annotation
ILLUMINA_ANNO_MAIN		:= $(patsubst $(ROOT_DIR)/%.txt,$(OUT_DIR)/annotation/%_main.txt,$(ILLUMINA_ANNO_FILE))
ILLUMINA_ANNO_CONTROL		:= $(patsubst $(ROOT_DIR)/%.txt,$(OUT_DIR)/annotation/%_control.txt,$(ILLUMINA_ANNO_FILE))
ILLUMINA_ANNO_SPLIT_SCRIPT	:= parseIlluminaAnnotation.pl

# Read the Illumina data files

RAW_LUMI			:= $(OUT_DIR)/lumi/raw.rds
PARSE_LUMI_SCRIPT		:= readIllumina.R

# Normalise the lumiBatch to make an ExpressionSet

LUMI_BG_METHOD			:= bgAdjust
LUMI_NORM_METHOD		:= quantile
LUMI_VS_METHOD			:= vst
LUMI_PREPROCESS_SCRIPT		:= lumiExpresso.R

PROCESSED_ESET			:= $(OUT_DIR)/lumi/normalised.rds

# Get the matrix from the eset

EXTRACT_SCRIPT			:= extractMatrix.R
EXPRESSION_DIR			:= $(OUT_DIR)/expression
EXPRESSION_FILES		:= $(EXPRESSION_DIR)/probe/raw.csv $(EXPRESSION_DIR)/control/raw.csv $(EXPRESSION_DIR)/probe/normalised.csv

# Differential expression

DE_DIR				:= $(OUT_DIR)/diffexp
DIFFEXP_FILES			:= $(DE_DIR)/probe/pvals.csv $(DE_DIR)/probe/qvals.csv
LIMMA_SCRIPT			:= arrayLimma.R

# Annotation

ANNOTATION_DIR			:= $(OUT_DIR)/annotation
ANNOTATION_FILES		:= $(ANNOTATION_DIR)/probe.csv $(ANNOTATION_DIR)/control.csv
ANNOTATION_SCRIPT		:= getAnnotation.R
BIOMART_ID_FIELD		:= illumina_humanht_12_v4
BIOMART_ENTREZ_GENE		:= entrezgene
BIOMART_SYMBOL			:= external_gene_name
BIOMART_GENE_FIELDS		:= ensembl_gene_id,entrezgene,chromosome_name,strand,start_position,end_position,external_gene_name,description,gene_biotype
ANNOTATE_CONTROL_SCRIPT		:= annotateControls.R

# Gene set analysis with mroast()

ROAST_SCRIPT			:= arrayLimmaRoast.R
ROAST_NROT			:= 9999
GENE_SET_RESULTS		:= $(DE_DIR)/mroast/results.csv

# Shiny interctive output

SHINY_CONFIG		:= config.txt
SHINY_OBJECT		:= $(OUT_DIR)/data.rds
SHINY_OBJECT_SCRIPT	:= makeShiny.R

.PHONY:split_anno
split_anno:$(ILLUMINA_ANNO_CONTROL) $(ILLUMINA_ANNO_MAIN)

$(ANNO_DIR)/%_main.txt $(ANNO_DIR)/%_control.txt:$(ANALYSIS_DIR)/%.txt
	mkdir -p $(dir $@) && perl $(ILLUMINA_ANNO_SPLIT_SCRIPT) $< $(ANNO_DIR)/$*_main.txt $(ANNO_DIR)/$*_control.txt

.PHONY:read_lumi
read_lumi:$(RAW_LUMI)

$(RAW_LUMI):$(SAMPLE_PROBE_PROFILE) $(CONTROL_PROBE_PROFILE) $(SAMPLES_TABLE) $(ILLUMINA_ANNO_MAIN) $(ILLUMINA_ANNO_CONTROL) $(EXPERIMENT)
	mkdir -p $(dir $@) && Rscript $(PARSE_LUMI_SCRIPT) $^ $@ > $(subst .rds,.log,$@) 2>&1

.PHONY:preprocess_lumi
preprocess_lumi:$(PROCESSED_ESET)

$(PROCESSED_ESET):$(RAW_LUMI)
	Rscript $(LUMI_PREPROCESS_SCRIPT) $< $(LUMI_BG_METHOD) $(LUMI_NORM_METHOD) $(LUMI_VS_METHOD) $@ > $(subst .rds,.log,$@) 2>&1

.PHONY:extract_matrices
extract_matrices:$(EXPRESSION_FILES)

%/expression/probe/raw.csv %/expression/control/raw.csv:%/lumi/raw.rds
	mkdir -p $(EXPRESSION_DIR)/probe && mkdir -p $(EXPRESSION_DIR)/control && Rscript $(EXTRACT_SCRIPT) $< $*/expression/probe/raw.csv NA $*/expression/control/raw.csv > $(subst .csv,.log,$@) 2>&1

%/expression/probe/normalised.csv:%/lumi/normalised.rds
	mkdir -p $(dir $@) && Rscript $(EXTRACT_SCRIPT) $< $*/expression/probe/normalised.csv $(LUMI_VS_METHOD) NA > $(subst .csv,.log,$@) 2>&1

.PHONY:run_limma
run_limma: $(DIFFEXP_FILES) 

%/diffexp/probe/pvals.csv %/diffexp/probe/qvals.csv: %/lumi/normalised.rds
	mkdir -p $(dir $@) && Rscript $(LIMMA_SCRIPT) $< $(EXPERIMENT) $(CONTRASTS) $*/diffexp/probe/pvals.csv $*/diffexp/probe/qvals.csv $(BATCH_VARIABLE) > $*/diffexp/probe.log 2>&1

.PHONY:get_annotation
get_annotation: $(ANNOTATION_FILES)

$(ANNOTATION_DIR)/control.csv:$(EXPRESSION_DIR)/control/raw.csv $(ILLUMINA_ANNO_CONTROL)
	Rscript $(ANNOTATE_CONTROL_SCRIPT) $^ $@

$(ANNOTATION_DIR)/probe.csv:$(EXPRESSION_DIR)/probe/raw.csv
	Rscript $(ANNOTATION_SCRIPT) $(SPECIES) $< $(BIOMART_GENE_FIELDS) $(BIOMART_ID_FIELD) $@ > $(subst .csv,.log,$@) 2>&1

# Run a gene set enrichment analysis with Limma's Roast

.PHONY: run_roast
run_roast: $(GENE_SET_RESULTS)

$(GENE_SET_RESULTS): $(PROCESSED_ESET) $(ANNOTATION_DIR)/probe.csv $(DIFFEXP_FILES)
	mkdir -p $(dir $@) && Rscript $(ROAST_SCRIPT) $< $(EXPERIMENT) $(CONTRASTS) $(word 2,$^) $(BIOMART_ENTREZ_GENE) $(BIOMART_ID_FIELD) $(ROAST_NROT) $(THREADS) $(MSIGDB_DIR) $@ $(BATCH_VARIABLE) 

#> $(subst .csv,.log,$@) 2>&1

# Make an object to use in a shiny

.PHONY: make_shiny_object
make_shiny_object: $(SHINY_OBJECT)

$(SHINY_OBJECT): $(ANNOTATION_FILES) $(CONTRASTS) $(EXPRESSION_FILES) $(DIFFEXP_FILES) $(GENE_SET_RESULTS) $(SHINY_CONFIG)
	mkdir -p $(dir $@) && Rscript $(SHINY_OBJECT_SCRIPT) $(CONTRASTS) $(EXPRESSION_DIR) $(DE_DIR) $(ANNO_DIR) $(GENE_SET_RESULTS) $(SHINY_CONFIG) $(SHINY_OBJECT) > $(subst .R,.log,$(SHINY_OBJECT_SCRIPT)) 2>&1

