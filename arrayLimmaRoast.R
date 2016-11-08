library(limma)
library(GSEABase)
library(lumi)

cl <- commandArgs(trailingOnly=TRUE)

eset_file <- cl[1]
experimentfile <- cl[2]
contrastsfile <- cl[3]
annotationfile <- cl[4]
entrezgenefield <- cl[5]
idfield <- cl[6]
nrot <- as.numeric(cl[7])
threads <- as.numeric(cl[8])
msigdbdir <- cl[9]
outfile <- cl[10]
batchvars <- cl[11]

experiment <- read.delim(experimentfile, sep='\t', stringsAsFactors=F)
allcontrasts <- read.delim(contrastsfile, sep='\t', stringsAsFactors=F)
batchvars <- unlist(strsplit(batchvars, ','))
batchvars <- batchvars[batchvars %in% colnames(experiment)]

eset <- readRDS(eset_file)

annotation <- read.csv(annotationfile, row.names=NULL, stringsAsFactors = FALSE)
annotation[[entrezgenefield]] <- as.character(annotation[[entrezgenefield]])

# Get the gene sets

msigdb_files <- list.files(msigdbdir, full.names=TRUE, pattern='.gmt')

gene_sets <- lapply(msigdb_files, getGmt)
names(gene_sets) <-  sub('.entrez.gmt', '', basename(msigdb_files))

# Remove anything with a tiny gene set

gene_sets <- lapply(gene_sets, function(pgss){
  pgss[ unlist(lapply(pgss, function(x) length(geneIds(x)))) >= 5 ]
})

# Convert to probe IDs

print("Converting gene sets to probes")

gene_sets <- lapply(gene_sets, function(gene_set_collection) {
  
  # gene_set_collection doesn't behave exactly like a list (it's a GSEABase object), so we have to make sure the result
  # gets named properly
  
  gsc <- lapply(gene_set_collection, function(gene_set) {
    set_gene_ids <- GSEABase::geneIds(gene_set)
    gs <- annotation[[idfield]][annotation[[entrezgenefield]] %in% set_gene_ids]
    gs[!is.na(gs)]
  })
  names(gsc) <- names(gene_set_collection)
  gsc
})

print("Done gene set conversion")

# This one isn't very useful

gene_sets <- gene_sets[names(gene_sets) != "c2.all.v5.0"]

roastres <- lapply(unique(allcontrasts$variable), function(contrast_variable){
  
  print(contrast_variable)
  
  # Exclude samples with blank values for this variable
  
  subexperiment <- experiment[! is.na(experiment[[contrast_variable]]),, drop = FALSE]
  subexperiment[[contrast_variable]] <- factor(subexperiment[[contrast_variable]])
  
  subeset <- eset[,rownames(subexperiment), drop = FALSE]

  contrasts <- allcontrasts[allcontrasts$variable == contrast_variable,]
  contrast_names_for_output <- apply(contrasts, 1, function(contrast) paste(contrast[1], paste(contrast[-1], collapse='-'), sep=':'))
  
  if (length(batchvars) > 0){
    contmodel <- paste('~0', contrast_variable, paste(batchvars, collapse = '+'), sep='+') 
    for (bv in batchvars){
      subexperiment[[bv]] <- factor(subexperiment[[bv]])
    }
  }else{
    print("No batch variable provided")
    contmodel <- paste('~0', contrast_variable, sep='+') 
  }
  
  print(paste('Model:', contmodel))
  
  design <- model.matrix( as.formula(contmodel), data=subexperiment)
  colnames(design) <- sub(contrast_variable, paste0(contrast_variable, '.'), colnames(design))
 
  fit <- lmFit(subeset, design)
 
  # Contrasts bit

  contrast_names <- paste(paste(contrast_variable, make.names(contrasts$group1), sep="."), paste(contrast_variable, make.names(contrasts$group2), sep="."), sep="-")[which(contrasts$variable == contrast_variable)] # for limma et al
  contrast_names <- gsub(".X", ".", contrast_names)
  contrast.matrix <- makeContrasts(contrasts=contrast_names, levels=design)
  
  mroast <- lapply(1:nrow(contrasts), function(n){
    
    print(contrast_names[n])
    
    contrast <- as.character(contrasts[n,])

    gsnames <- names(gene_sets)
    names(gsnames) <- gsnames
    
    #clust <- makeCluster(getOption("cl.cores", threads))
    #clusterExport(clust, c('gsnames', 'contrast.matrix', 'gene_sets', 'subeset', 'design', 'nrot', 'contrast_names_for_output', 'n'), envir=environment())
    #clusterEvalQ(clust, {library(limma)})
    #clusterEvalQ(clust, {library(lumi)})
    
    #gsres <- parLapply(clust, gsnames, function(pgss){
    gsres <- lapply(gsnames, function(pgss){
      print(paste0('...', pgss))
      
      # Reverse the sense of the contrast so it makes sense
      
      mrcont <- contrast.matrix[,n]
      mrcont[contrast.matrix[,n] == -1] <- 1
      mrcont[contrast.matrix[,n] == 1] <- -1
      
      res <- mroast(
        y=subeset, 
        index=ids2indices(gene_sets[[pgss]], rownames(subeset)),
        design=design,
        contrast=mrcont,
        nrot=nrot)
      
      data.frame(
        cbind(
          contrast_names_for_output[n], 
          gene_set_type=pgss, 
          gene_set=rownames(res), 
          res,
          row.names = NULL
        ), 
        row.names=NULL
      )  
      
    })
    #stopCluster(clust)
    
    do.call(rbind, gsres)
    
  })
  do.call(rbind, mroast)
})

final <- do.call(rbind, roastres)
colnames(final)[1] <- 'contrast'
colnames(final) <- gsub('\\.', '_', colnames(final))

write.csv(final, file=outfile, row.names=F)
