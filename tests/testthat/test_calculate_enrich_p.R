# Test file for calculate_enrich_p function
library(BiFET)
context(desc = "Calculate enrich p")

# load peak file
peak_file <- system.file("extdata", "input_peak_motif.Rdata",
                         package = "BiFET")
load(peak_file)
TFbinding.mat <- bindingTF_per_peak(GRpeaks, GRmotif)
targetpeak <- which(GRpeaks$peaktype == "target")
backgroundpeak <- which(GRpeaks$peaktype == "background")
reads <- GRpeaks$reads
GCcontent <- GRpeaks$GC

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


