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


# References

* [1] Du, P., Kibbe, W.A., Lin and S.M. (2008). &ldquo;lumi: a pipeline for processing Illumina microarray.&rdquo; <em>Bioinformatics</em>.  </p>  <p>P D, X Z, CC H, N J, WA K, L H and SM L (2010). &ldquo;Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis.&rdquo; <em>BMC Bioinformatics</em>.  </p>  <p>Lin, S.M., Du, P., Kibbe and W.A. (2008). &ldquo;Model-based Variance-stabilizing Transformation for Illumina Microarray Data.&rdquo; <em>Nucleic Acids Res</em>.  </p>  <p>Du, P., Kibbe, W.A., Lin and S.M. (2007). &ldquo;nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays.&rdquo; <em>Biology Direct</em>.
* [2] Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). &ldquo;limma powers differential expression analyses for RNA-sequencing and microarray studies.&rdquo; <em>Nucleic Acids Research</em>, <b>43</b>(7), pp. e47.
* [3] >Huber, W., Carey, J. V, Gentleman, R., Anders, S., Carlson, M., Carvalho, S. B, Bravo, C. H, Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, D. K, Irizarry, A. R, Lawrence, M., Love, I. M, MacDonald, J., Obenchain, V., Ole's, K. A, Pag'es, H., Reyes, A., Shannon, P., Smyth, K. G, Tenenbaum, D., Waldron, L., Morgan and M. (2015). &ldquo;Orchestrating high-throughput genomic analysis with Bioconductor.&rdquo; <em>Nature Methods</em>, <b>12</b>(2), pp. 115&ndash;121. <a href=\"http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html\">http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html</a>. 
