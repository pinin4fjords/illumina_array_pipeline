library(biomaRt)

cl <- commandArgs(trailingOnly=TRUE)

species <- cl[1]
expr_matrix <- cl[2]
fields <- unlist(strsplit(cl[3],','))
idfield <- cl[4]
output <- cl[5]

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = paste0(species, '_gene_ensembl'), host='www.ensembl.org')
t2g <- getBM(attributes = unique(c(idfield, fields)), filters = paste0("with_", idfield), values = TRUE, mart = mart)

expr <- as.matrix(read.csv(expr_matrix, row.names=1, check.names = FALSE))
t2g <- t2g[t2g[[idfield]] %in% rownames(expr), ]

# Put genes with Entrez IDs first, so that when we pick the first of a set it's
# one with an Entrez ID where available.

t2g <- rbind(subset(t2g, ! is.na(entrezgene)), subset(t2g, is.na(entrezgene)))

genes_by_probe <- split(t2g, t2g[[idfield]])

multis <- which(unlist(lapply(genes_by_probe, function(x) nrow(x) > 1)))

genes_by_probe[multis]  <- lapply(genes_by_probe[multis], function(x){
  
  if (length(unique(x$ensembl_gene_id)) == 1){
    x[1,, drop = FALSE]
  }else if (sum(x$gene_biotype == 'protein_coding') == 1){
    subset(x, gene_biotype == 'protein_coding')
  }else if (sum(! is.na(x$entrezgene)) == 1){
    subset(x, ! is.na(entrezgene))
  }else if (length(unique(x$family)) == 1){
    x[1, ,drop = FALSE]
  }else{
    NULL
  }
})

genes_by_probe <- genes_by_probe[unlist(lapply(genes_by_probe, function(x) ! is.null(x)))]

anno <- data.frame(rownames(expr), stringsAsFactors = FALSE)
colnames(anno) <- idfield

anno <- merge(anno, data.table::rbindlist(genes_by_probe), all.x = TRUE, sort = FALSE)
anno <- anno[match(rownames(expr), anno[[idfield]]), ]

write.csv(anno, file=output, row.names=FALSE)
