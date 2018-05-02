# Test file for fast_iasva.R function
library(iasva)
library(SummarizedExperiment)
context(desc = "fast_iasva")

# load test data
counts_file <- system.file("extdata", "iasva_counts_test.Rds", package = "iasva")
counts <- readRDS(counts_file)
anns_file <- system.file("extdata", "iasva_anns_test.Rds", package = "iasva")
anns <- readRDS(anns_file)
Geo_Lib_Size <- colSums(log(counts + 1))
Patient_ID <- anns$Patient_ID
mod <- model.matrix(~Patient_ID + Geo_Lib_Size)
# create summarized experiment object
summ_exp <- SummarizedExperiment(assays = counts)

# test that input is summarized experiment
test_that("correct input format", {
  expect_equal(object = class(summ_exp)[1], expected = "SummarizedExperiment")
  expect_gt(object = nrow(assay(summ_exp)), expected = 1)
  expect_gt(object = ncol(assay(summ_exp)), expected = 1)
})

iasva.res <- fast_iasva(summ_exp, mod[, -1], num.sv = 5)
# test that output is list
test_that("correct output results", {
  expect_type(object = iasva.res, type = "list")
  expect_equal(object = length(iasva.res), expected = 3)
})

