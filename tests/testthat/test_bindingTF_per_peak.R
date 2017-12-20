# Test file for bindingTF_per_peak function
library(BiFET)
context(desc = "Binding TF per peak")

# load peak file
peak_file <- system.file("extdata", "islet_PBMC_consensus_peak.Rdata",
                         package = "BiFET")
load(peak_file)
peaks_sel <- peaks[1:1000,]
footprints <- system.file("extdata", "PBMC_PIQ.bed", package = "BiFET")
PBMCmotif <- read.delim(footprints, header = FALSE)
PBMCmotif <- motif_processing(PBMCmotif)
tf_res <- bindingTF_per_peak(peaks = peaks_sel, motif = PBMCmotif)

# test that output matrix has same number of rows as peak matrix
test_that("input matrix row number is same as output matrix column number", {
  expect_equal(object = ncol(tf_res), expected = nrow(peaks_sel))
})

