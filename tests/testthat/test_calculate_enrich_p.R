# Test file for calculate_enrich_p function
library(BiFET)
context(desc = "Calculate enrich p")

# load peak file
peak_file <- system.file("extdata", "islet_PBMC_consensus_peak.Rdata",
                         package = "BiFET")
load(peak_file)
peaks_sel <- peaks[1:5000,]
footprints <- system.file("extdata", "PBMC_PIQ.bed", package = "BiFET")
PBMCmotif <- read.delim(footprints, header = FALSE)
PBMCmotif <- motif_processing(PBMCmotif)
TFbinding.mat <- bindingTF_per_peak(peaks_sel, PBMCmotif)
targetpeak <- which(peaks_sel[, "peaktype"] == "target")
backgroundpeak <- which(peaks_sel[, "peaktype"] == "background")
reads <- peaks_sel[, "reads"]
GCcontent <- peaks_sel[, "GCcontent"]

test_that("correct object types",  {
  expect_type(object = reads, type = "double")
  expect_type(object = GCcontent, type = "double")
  expect_type(object = targetpeak, type = "integer")
  expect_type(object = backgroundpeak, type = "integer")
})


result <- calculate_enrich_p(TFbinding.mat, reads,
                             GCcontent, targetpeak, backgroundpeak)

test_that("correct output results", {
  expect_type(object = result, type = "list")
  expect_equal(object = length(result), expected = 3)
})


