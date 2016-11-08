cl <- commandArgs(trailingOnly = TRUE)

library(shiny)
library(shinyngs)
library(SummarizedExperiment)

contrastsfile <- cl[1]
expressiondir <- cl[2]
testsdir <- cl[3]
annotationdir <- cl[4]
gene_set_table_file <- cl[5]
configfile <- cl[6]
shinyobject <- cl[7]

print("Reading inputs")

config_lines <- readLines(configfile)
config <- list(experiments = c())

# Parse the config file

section <- NA
for (line in config_lines){
  l <- unlist(strsplit(line, '\t'))
  field <- l[1]
  value <- l[2]
  
  if (grepl('^#', field)){
    section <- sub('#', '', field)
    config[[section]] <- list()
    config$experiments <- c(config$experiments, section)
    next
  }else{
    field_names_and_values <- do.call(rbind, unlist(lapply(strsplit(value,','), function(x) strsplit(x, '\\|')), recursive=F))
    if (ncol(field_names_and_values) == 1){
      v <- field_names_and_values[,1]
    }else{
      v <- field_names_and_values[,2]
      v <- setNames(v, field_names_and_values[,1])
    }
    
    if (is.na(section)){
      config[[field]] <- v
    }else{
      config[[section]][[field]] <- v
    }
  }
}

description <- readChar(config$description, file.info(config$description)$size)

contrasts <- read.delim(contrastsfile, stringsAsFactors=FALSE)

# Read in the expression matrices

experiments <- structure(config$experiments, names = config$experiments)
experiment_data <- lapply(experiments, function(x) read.delim(config[[x]]$experiment))

# Create assays with matching rows, filling with emtpy values where necessary
# NOTE: we must match assay columns to experiment rows- it will not be done
# automatically

print("Reading expression...")

expression <- lapply(experiments, function(exp){
  assaydetails <- file.info(list.files(file.path(expressiondir, exp), full.names = TRUE, pattern = '.csv'))
  assaydetails = assaydetails[with(assaydetails, order(as.POSIXct(mtime), decreasing = TRUE)), ]
  assayfiles <- rownames(assaydetails)
  
  assays <- lapply(assayfiles, function(x){
	as.matrix(read.csv(x, row.names=1, check.names = FALSE))[,rownames(experiment_data[[exp]])]
  })
  names(assays) <- basename(tools::file_path_sans_ext(assayfiles))
  
  # The assays slot of a summarised experiment needs the same dimensions for every matrix
  
  all_rows <- Reduce(union, lapply(assays, rownames))
  
  lapply(assays, function(as){
    
    missing_rows <- all_rows[! all_rows %in% rownames(as)]
    empty_rows <- data.frame(matrix(NA, nrow = length(missing_rows), ncol = ncol(as)), row.names = missing_rows)
    colnames(empty_rows) <- colnames(as)
    as <- rbind(as, empty_rows)
    round(as.matrix(as[all_rows,]), 2)
  })
})

# Derive annotation matrices to match the expression

print ("Reading annotation")

experiment_annotation <- lapply(experiments, function(exp){
  anno <- read.csv(file.path(annotationdir, paste(exp, 'csv', sep='.')), stringsAsFactors=FALSE)
  anno[match(rownames(expression[[exp]][[1]]), anno[[config[[exp]]$idtype]]), ]
})

# Contrasts

contrast_names <- apply(contrasts, 1, function(x) paste(x[1], paste(x[-1], collapse = '-'), sep=':'))
contrast_names <- structure(contrast_names, names = contrast_names)

# Make the ExploratorySummarisedExperiments

print("Constructing ExploratorySummarizedExperiments")

expsumexps <- lapply(experiments, function(exp){
  
  ese <- ExploratorySummarizedExperiment(
    assays = SimpleList(
      expression[[exp]]
    ),
    colData = DataFrame(experiment_data[[exp]]),
    annotation <- experiment_annotation[[exp]],
    idfield = config[[exp]]$idtype,
    entrezgenefield = "entrezgene",
    labelfield = "external_gene_name",
    assay_measures = as.list(config[[exp]]$assay_measures)
  )
  
  # Rename the samples if specified
  
  if ('name_by' %in% names(config)){
    colnames(ese) <- experiment_data[[exp]][colnames(ese), config$name_by]
  }
  ese
})


# Tests: a list of statistics (e.g. p/q values). List names must be the same
# (or a subset of) those used in the assays slot.

print("Loading derived statstics")

if (file.exists(testsdir)){
  test_stats <- structure(c('pvals', 'qvals'), names = c('pvals', 'qvals'))
  tests <- lapply(experiments, function(exp){
    stats <- list(lapply(test_stats, function(ts){
      statsfile <- file.path(testsdir, exp, paste(ts, 'csv', sep='.'))
      if (file.exists(statsfile)){
        test <- read.csv(statsfile, row.names = 1)
        test[rownames(expression[[exp]][[1]]),,drop = FALSE]
      }else{
        NULL
      }
      
    }))
    names(stats) <- names(expression[[exp]])[1]
    stats
  })
  
  # Add to the objects
  
  for (exp in names(tests)){
    expsumexps[[exp]]@tests = tests[[exp]]
  }
}

# Gene sets

gene_sets_for_object <- list()

if ('gene_sets' %in% names(config)){

  print("Loading gene sets and associated analysis")
  
  gene_sets <- config$gene_sets
  
  gene_set_ids <- sub('.entrez', '', basename(tools::file_path_sans_ext(gene_sets)))
  names(gene_set_ids) <- names(gene_sets)
  
  gene_sets_for_object <- lapply(gene_sets, GSEABase::getGmt)
  
  if (file.exists(gene_set_table_file)){
    gene_set_table <- read.csv(gene_set_table_file, stringsAsFactors = FALSE)
    
    gene_set_table[, c('PropDown', 'PropUp')] <- round(gene_set_table[, c('PropDown', 'PropUp')], 2)
    
    gene_set_names <- structure(names(gene_sets), names = names(gene_sets))
    gst <- lapply(gene_set_names, function(gsn){
      lapply(contrast_names, function(cn){
        tab <- subset(gene_set_table, contrast == cn & gene_set_type == gene_set_ids[[gsn]])
        rownames(tab) <- tab$gene_set
        tab[,c('NGenes', 'PropDown', 'PropUp', 'Direction', 'PValue', 'FDR')]
      })
    })
    
    # Make a list for the GSA
    
    mainexperiment <- names(expsumexps)[1]
    mainassay <- names(assays(expsumexps[[mainexperiment]]))[1]
    newgsa <- list()
    newgsa[[mainassay]] <- gst
    
    expsumexps[[mainexperiment]]@gene_set_analyses <- newgsa
  
  }
}

print("Creating ExploratorySummarizedExperimentList")

expsumexpslist <- ExploratorySummarizedExperimentList(
  expsumexps,
  title = config$title,
  author = config$author,
  description = as.character(includeMarkdown(config$description)),
  group_vars = config$group_vars,
  default_groupvar = config$default_groupvar,
  contrasts = lapply(1:nrow(contrasts), function(x) as.character(contrasts[x,])),
  url_roots = as.list(config$url_roots),
  gene_sets = gene_sets_for_object
)

# Species 

if ('ensembl_species' %in% names(config)){
  expsumexpslist@ensembl_species <- config$ensembl_species
}

# Read attrition report etc

for (exp in names(expsumexpslist)){
  if ('read_reports' %in% names(config[[exp]])){
    rrfs <- lapply(config[[exp]]$read_reports, function(x) unlist(strsplit(x, ';')))
    names(rrfs) <- basename(tools::file_path_sans_ext(rrfs))
      
      slot(expsumexpslist[[exp]], 'read_reports') <- lapply(rrfs, function(rrf){
        readrep <- as.matrix(read.csv(rrf, row.names = 1, check.names = F))[rownames(experiment_data[[exp]]),]
        if ('name_by' %in% names(config)){
          rownames(readrep) <- experiment_data[[exp]][rownames(readrep), config$name_by]
        }
        readrep
      })
    }
}

if ('read_reports' %in% names(config)){
 
  print("Loading read reports")
 
  read_report_files <- lapply(parse_multifield(config$read_reports), function(x) unlist(strsplit(x, ';')))
  
  for (exp in names(read_report_files)){
    rrfs <- read_report_files[[exp]] 
    names(rrfs) <- basename(tools::file_path_sans_ext(rrfs))
    
    slot(expsumexpslist[[exp]], 'read_reports') <- lapply(rrfs, function(rrf){
      readrep <- as.matrix(read.csv(rrf, row.names = 1, check.names = F))[rownames(experiment_data[[exp]]),]
      if ('name_by' %in% names(config)){
        rownames(readrep) <- experiment_data[[exp]][rownames(readrep), config$name_by]
      }
      readrep
    })
  }
}

saveRDS(expsumexpslist, file = shinyobject)
