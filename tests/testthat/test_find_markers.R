# Test file for find.markers.R function
library(iasva)
context(desc = "find_markers")

# load test data
counts_file <- system.file("extdata", "iasva_counts_test.Rds", package = "iasva")
counts <- readRDS(counts_file)
anns_file <- system.file("extdata", "iasva_anns_test.Rds", package = "iasva")
anns <- readRDS(anns_file)
Geo_Lib_Size <- colSums(log(counts + 1))
Patient_ID <- anns$Patient_ID
mod <- model.matrix(~Patient_ID + Geo_Lib_Size)

# test that read counts input is matrix
test_that("correct input matrix format", {
  expect_gt(object = nrow(counts), expected = 1)
  expect_gt(object = ncol(counts), expected = 1)
})

iasva.res <- iasva(t(counts), mod[, -1], num.sv = 5, permute = FALSE) 
# test that output is list
test_that("correct output results", {
  expect_type(object = iasva.res, type = "list")
  expect_equal(object = length(iasva.res), expected = 4)
})

markers <- find_markers(t(counts), iasva.res$sv)
# test output of markers is matrix
test_that("correct output marker matrix format", {
  expect_gt(object = nrow(markers), expected = 1)
  expect_gt(object = ncol(markers), expected = 1)
})



