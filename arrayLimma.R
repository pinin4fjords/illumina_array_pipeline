library(limma)
library(edgeR)
library(multtest)
library(lumi)

cl <- commandArgs(trailingOnly=TRUE)

eset_file <- cl[1]
experimentfile <- cl[2]
contrastsfile <- cl[3]
pvalsfile <- cl[4]
qvalsfile <- cl[5]
batchvars <- cl[6]

print(paste('pvals file', pvalsfile))
print(paste('qvals file', qvalsfile))

experiment <- read.delim(experimentfile, sep='\t', stringsAsFactors=F)
contrasts <- read.delim(contrastsfile, sep='\t', stringsAsFactors=F)
eset <- readRDS(eset_file)

batchvars <- unlist(strsplit(batchvars, ','))
batchvars <- batchvars[batchvars %in% colnames(experiment)]

results <- do.call(cbind, lapply(unique(contrasts$variable), function(contrast_variable){

  # Exclude samples with blank values for this variable

  subexperiment <- experiment[! is.na(experiment[[contrast_variable]]),, drop = FALSE]
  subexperiment[[contrast_variable]] <- factor(subexperiment[[contrast_variable]])
 
  subeset <- eset[,rownames(subexperiment),drop = FALSE]
  
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
  
  contrast_names_for_output <- apply(contrasts, 1, function(contrast) paste(contrast[1], paste(contrast[-1], collapse='-'), sep=':'))[which(contrasts$variable == contrast_variable)]
 
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  out <- do.call(cbind, lapply(1:ncol(contrast.matrix), function(k){
	multadjust <- mt.rawp2adjp(fit2$p.value[,k], proc=c("BH"))
    	data.frame(
	    'PValue' = fit2$p.value[,k], 
	    'FDR' = (multadjust$adjp[order(multadjust$index),])[,2]	
	)
  }))
  colnames(out) <- paste(colnames(out), unlist(lapply(contrast_names_for_output, function(x) rep(x, 2))), sep=':')

  out
}))

print("Done analysis")

pvals <- results[,grep('PValue', colnames(results)), drop = FALSE]
qvals <- results[,grep('FDR', colnames(results)), drop = FALSE]

colnames(pvals) <- sub('PValue:', '', colnames(pvals))
colnames(qvals) <- sub('FDR:', '', colnames(qvals))

print(paste('Minimum p value:', min(pvals)))
print(paste('Minimum q value:', min(qvals)))

# Output the table of p values

write.csv(pvals, file=pvalsfile)
write.csv(qvals, file=qvalsfile)

