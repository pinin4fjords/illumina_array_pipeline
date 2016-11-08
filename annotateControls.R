##############################################################################
#
# Annotate control probes based on the Illumina file. Control probes can
# fall in multiple categories
#
##############################################################################

cl <- commandArgs(trailingOnly = TRUE)

library(data.table)

expression_file <- cl[1]
control_annotation_file <- cl[2]
outfile <- cl[3]

expression <- read.csv(expression_file, row.names = 1, check.names = FALSE)
control_annotation <- read.delim(control_annotation_file, stringsAsFactors = FALSE)

control_annotation <- rbindlist(lapply(split(control_annotation, control_annotation$Array_Address_Id), function(x){
  if (nrow(x) > 1){
    z <- x[1,,drop=F]
    z[1,] <- apply(x,2,function(y) paste(unique(y), collapse = ','))
    z
  }else{
    x
  }
}))

control_annotation <- control_annotation[match(rownames(expression), control_annotation$Array_Address_Id),]

write.csv(control_annotation, file = outfile, row.names = FALSE)
