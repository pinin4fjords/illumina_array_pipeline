#############################################
#
# Take a set of .CEL files, output affyBatch
#
#############################################

# Galaxy is sensitive to output messages, so make sure we control what's output

library(affy)

cl = commandArgs(trailingOnly=T)

sample_probe_profile = cl[1]
control_probe_profile = cl[2]
sample_info_file = cl[3]
illumina_annotation_main = cl[4]
illumina_annotation_control = cl[5]
experimentfile = cl[6]
outfile = cl[7]

res = ''

# First extract the tar file to get the CEL files

# Then read in the experiment

print (paste("Reading experiment file ", experimentfile))
experiment <- read.delim(file=experimentfile, sep="\t", row.names=1, stringsAsFactors=FALSE)
print ("Done reading experiment")

# Now read the expression files

library(lumi)

# Now read in, dending on the optional ones

if (sample_info_file != 0){
    result <- lumiR.batch(sample_probe_profile, sampleInfoFile=sample_info_file, convertNuID=FALSE)
}else{
    result <- lumiR.batch(sample_probe_profile, convertNuID=FALSE)
}

############
#
# Check IDs 
#
############

# if these are not ILMN... ids, then the annotation packages don't work. Can be converted with the appropriate information extracted from the files from illumina at http://www.switchtoi.com/annotationfiles.ilmn
    
if (length(grep("ILMN", featureNames(result))) == 0){
    
    if (illumina_annotation_main != 0){
    		
        illumina_anno <- read.delim(file=illumina_annotation_main, sep="\t", header=TRUE)
    	featureNames(result) <- illumina_anno$Probe_Id[match(featureNames(result), illumina_anno$Array_Address_Id)]
    
    }else{
        stop("This dataset has non-ILMN IDs, almost certainly 'Arrray Address Ids'. Go to http://www.switchtoi.com/annotationfiles.ilmn and download the annotation file from your chip, upload it with the upload tool, and parse it with 'Parse Illumina annoation...'")
    }
}

#########

# Re-arrange the data to agree with the experiment. Do this BEFORE adding the controls, because it messes up the controlData slot

result <- result[,match(rownames(experiment), rownames(pData(result)))]

# Add in the controls

if (control_probe_profile != 0){
    result <- addControlData2lumi(control_probe_profile, result)
}

# Rename to the 'Name' column in two places

if (is.null(experiment$Name)){
    experiment$Name <- rownames(experiment)
}

pData(phenoData(result))$sampleID <- experiment$Name[match(pData(phenoData(result))$sampleID, rownames(experiment))] 
sampleNames(result) <- experiment$Name[match(sampleNames(result), rownames(experiment))] 

# Now we do the output

saveRDS(result, file = outfile)
