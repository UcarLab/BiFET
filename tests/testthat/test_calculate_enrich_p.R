# Test file for calculate_enrich_p function
library(BiFET)
context(desc = "Calculate enrich p")

# load peak file
peak_file <- system.file("extdata", "input_peak_motif.Rdata",
                         package = "BiFET")
load(peak_file)
result <- calculate_enrich_p(GRpeaks, GRmotif)
test_that("correct output results", {
  expect_type(object = result, type = "list")
  expect_equal(object = length(result), expected = 3)
})


