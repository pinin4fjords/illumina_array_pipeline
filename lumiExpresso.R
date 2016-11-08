#############################################
#
# Take a lumiBatch, process with expresso, output eset
#
#############################################

cl = commandArgs(trailingOnly=T)
lumibatch = cl[1]
bg.method = cl[2]
norm.method = cl[3]
vs.method = cl[4]
outfile = cl[5]

library(affy)
library(lumi)
library(limma)
library(multtest)
library(vsn)

# Read the expressionset object

result <- readRDS(lumibatch)

# Use the options to determine how to run expresso

if (bg.method == 0){
    if (norm.method == 0){
        if (vs.method == 0){
            result <- lumiExpresso(result, bg.correct = FALSE, normalize = FALSE, variance.stabilize=FALSE, verbose = FALSE)
        }
        else{
            result <- lumiExpresso(result, bg.correct = FALSE, normalize = FALSE, variance.stabilize=TRUE, varianceStabilize.param = list(method=vs.method), verbose = FALSE)
        } 
    }else{
        if (vs.method == 0){
            result <- lumiExpresso(result, bg.correct = FALSE, normalize = TRUE, normalize.param = list(method=norm.method), variance.stabilize=FALSE, verbose = FALSE)
        }
        else{
            result <- lumiExpresso(result, bg.correct = FALSE, normalize = TRUE, normalize.param = list(method=norm.method), variance.stabilize=TRUE, varianceStabilize.param = list(method=vs.method), verbose = FALSE)
        } 
    }
}else{
    if (norm.method == 0){
        if (vs.method == 0){
            result <- lumiExpresso(result, bg.correct = TRUE, bgcorrect.param = list(method=bg.method), normalize = FALSE, variance.stabilize=FALSE, verbose = FALSE)
        }
        else{
            result <- lumiExpresso(result, bg.correct = TRUE, bgcorrect.param = list(method=bg.method), normalize = FALSE, variance.stabilize=TRUE, varianceStabilize.param = list(method=vs.method), verbose = FALSE)
        } 
    }else{
        if (vs.method == 0){
            result <- lumiExpresso(result, bg.correct = TRUE, bgcorrect.param = list(method=bg.method), normalize = TRUE, normalize.param = list(method=norm.method), variance.stabilize=FALSE, verbose = FALSE)
        }
        else{

	    print("Background correcting, normalising and variance stabilising")

            result <- lumiExpresso(result, bg.correct = TRUE, bgcorrect.param = list(method=bg.method), normalize = TRUE, normalize.param = list(method=norm.method), variance.stabilize=TRUE, varianceStabilize.param = list(method=vs.method), verbose = FALSE)
        } 
    }
}

# Now output the eset

saveRDS(result, file = outfile)

