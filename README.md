# Basic Illumina array pipeline

This is a simple pipeline for analysing Illumina array data, using the Lumi [1] and Limma [2] packages of Bioconductor. To use it, edit the indicated variables at the top of the lumi.mk makefile, and execute it with 'make -f lumi.mk'. A directory called 'pipeline' will be created and will contain the outputs.

# Inputs

This is not a bead-level analysis, it assumes you have the following:

* Sample probe profile, control probe profile and samples table, all usually exported from BeadStudio.
* An illumina annotation file, e.g. "HumanHT-12_V4_0_R2_15002873_B.txt" 
* A tab-delimted experiment file describing your samples and with rows matching columns of the input data.
* A tab-delimited file defining contrasts.
* A set of .gmt format gene set files for differential gene set analysis. 

## Experiment file

The experiment file is tab-delimited without a column name for sample IDs, like:

```
age	gender
Sample1	25	M
Sample2	30	F
Sample3 22	F
Sample4 12	M
Sample5	50	M
Sample6	70	F
```

## Contrasts file

The contrasts file defines contrasts in terms of the variables found in the experiment, like:

```
variable	group1	group2
gender	F	M
```

# Analysis

Makefile targtets are:

## all

Run all of the following. The default.

## split_anno

Split the Illumina annotation file into main- and control- probes

## read_lumi

Run readIllumina.R to Make a valid lumiBatch object from the inputs (will be serialised to .RDS)

## preprocess_lumi 

Run lumiExpresso.R to Call lumiExpresso() to perform background correction, normalisation and variance stabilisation (check for LUMI parameters in the makefile to tweak options). 

## extract_matrices

Use extractMatrix.R to derive csv-formatted matrices we can use later for exploratory purposes.

## run_limma

Run arrayLimma.R to look at the specified contrasts using limma and produce matrices of uncorrected and corrected p values.

## run_roast

Employ limma's mroast() method to perform differential gene set analysis.

## make_shiny_object

Using makeShiny.R, take the text-format outputs and make a data structure for use with [shinyngs](https://github.com/pinin4fjords/shinyngs). This will be serialised to data.rds, and can the be loaded for visualisation:

```R
eselist <- readRDS('data.rds')
app <- prepareApp('illuminaaarray', eselist)
shiny::shinyApp(ui = app$ui, server = app$server)
```

# References

* [1] Du, P., Kibbe, W.A., Lin and S.M. (2008). &ldquo;lumi: a pipeline for processing Illumina microarray.&rdquo; <em>Bioinformatics</em>.  </p>  <p>P D, X Z, CC H, N J, WA K, L H and SM L (2010). &ldquo;Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis.&rdquo; <em>BMC Bioinformatics</em>.  </p>  <p>Lin, S.M., Du, P., Kibbe and W.A. (2008). &ldquo;Model-based Variance-stabilizing Transformation for Illumina Microarray Data.&rdquo; <em>Nucleic Acids Res</em>.  </p>  <p>Du, P., Kibbe, W.A., Lin and S.M. (2007). &ldquo;nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays.&rdquo; <em>Biology Direct</em>.
* [2] Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). &ldquo;limma powers differential expression analyses for RNA-sequencing and microarray studies.&rdquo; <em>Nucleic Acids Research</em>, <b>43</b>(7), pp. e47.
* [3] >Huber, W., Carey, J. V, Gentleman, R., Anders, S., Carlson, M., Carvalho, S. B, Bravo, C. H, Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, D. K, Irizarry, A. R, Lawrence, M., Love, I. M, MacDonald, J., Obenchain, V., Ole's, K. A, Pag'es, H., Reyes, A., Shannon, P., Smyth, K. G, Tenenbaum, D., Waldron, L., Morgan and M. (2015). &ldquo;Orchestrating high-throughput genomic analysis with Bioconductor.&rdquo; <em>Nature Methods</em>, <b>12</b>(2), pp. 115&ndash;121. <a href=\"http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html\">http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html</a>. 
