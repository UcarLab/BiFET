# Test file for bindingTF_per_peak function
library(BiFET)
context(desc = "Binding TF per peak")

# load peak file
peak_file <- system.file("extdata", "input_peak_motif.Rdata",
                         package = "BiFET")
load(peak_file)
tf_res <- bindingTF_per_peak(GRpeaks = GRpeaks, GRmotif = GRmotif)

# test that output matrix has same number of rows as peak matrix
test_that("input matrix row number is same as output matrix column number", {
  expect_equal(object = ncol(tf_res), expected = length(GRpeaks))
})

