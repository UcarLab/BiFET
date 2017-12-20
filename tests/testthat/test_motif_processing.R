# Test file for motif_processing function
library(BiFET)
context(desc = "Motif processing")

# load in sample matrix
footprints <- system.file("extdata", "PBMC_PIQ.bed", package = "BiFET")
PBMCmotif <- read.delim(footprints, header = FALSE, stringsAsFactors = F)
res <- motif_processing(motif = PBMCmotif)

test_that("input matrix or data frame has 4 or more columns", {
  expect_gte(object = ncol(res), expected = 4)
  expect_type(object = res[,1], type = "character")
  expect_type(object = res[,2], type = "integer")
  expect_type(object = res[,3], type = "integer")
  expect_type(object = res[,4], type = "character")
})

# test that resulting matrix has fewer columns than input
test_that("output processed matrix has fewer columns that input matrix", {
  expect_lt(object = ncol(res), expected = ncol(PBMCmotif))
})
