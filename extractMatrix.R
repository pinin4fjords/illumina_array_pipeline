###############################################################################
#
# Derive an untransformed expresion matrix from an object like a LumiBatch
# or an ExpressionSet
#
# May require un-transformation
#
###############################################################################

cl <- commandArgs(trailingOnly = TRUE)

library(affy)
library(lumi)

objectfile <- cl[1]
mainfile <- cl[2]
transform_method <- cl[3]
controlfile <- cl[4]

object <- readRDS(objectfile)

if (transform_method == 'vst'){
  exprs <- exprs(lumiB(inverseVST(object), method="forcePositive"))
}else if (transform_method == 'log2'){
  exprs <- 2^exprs(object)
}else{
  exprs <- exprs(object)
}

write.csv(exprs, file = mainfile)

# Now the control data

if ('controlData' %in%  slotNames(object) && controlfile != 'NA'){

  control_data <- object@controlData
  
  # Chuck out duplicates where probe in two categories
  
  control_data <- control_data[match(unique(control_data$ProbeID), control_data$ProbeID),]
  
  rownames(control_data) <- control_data$ProbeID

  write.csv(control_data[,c(-1,-2)], file = controlfile)
}