# BiFET

Bias-free Footprint Enrichment Test (BiFET) is a statistical method 
to identify transcription factors (TFs) whose
footprints are over-represented in target regions compared to background 
regions after correcting for differences in read counts and GC contents between
target and background regions, where regions represent ATAC-seq or DNAse-seq peaks.

__Author__

Ahrim Youn

__Installation:__

To install this package, start R and enter:

`if (!requireNamespace("BiocManager", quietly=TRUE))`
    `install.packages("BiocManager")`

`BiocManager::install("BiFET")`

Alternatively, the most updated version of this package may be installed via the commands in R:

`library(devtools)`

`install_github('UcarLab/BiFET')`

__Publication__:

For more details about BiFET, please read our publication in _Nucleic Acids Reserarch_: https://www.ncbi.nlm.nih.gov/pubmed/30428075

__Tutorial__:

For instructions on how to use BiFET, please see the package vignette: 
https://bioconductor.org/packages/release/bioc/vignettes/BiFET/inst/doc/BiFET.html

__Disclaimer of Warranties and Liabilities__

The Jackson Laboratory provides the software “as is” without warranty of any kind, implied or expressed. You assume full responsibility and risk of loss resulting from your downloading and use of the content of the software. We expressly disclaim any warranty of merchantability, title, security, accuracy and non-infringement. In no event shall The Jackson Laboratory be liable for any claim, damages or other liability arising from the software or the use of the software. You may only use our content in academic research but not for commercial purposes. The software is provided as an information resource only, and should not be used or relied on for any diagnostic or treatment purposes.
