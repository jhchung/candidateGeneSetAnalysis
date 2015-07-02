library(testthat)
library(tbx1GeneNetwork)
library(plyr)
library(dplyr)

test_check("tbx1GeneNetwork")

test_that("calculate_empirical_pvalue gives the right result", {
  
  orig_test <- 5
  resampled_test <- c(1, 4, 5, 0, 2, 4, 7)
  
  expect_equal(calculate_empirical_pvalue(orig_test, resampled_test),
               0.375)
})

test_that("remove_prevously_selected_genes filters correctly", {
  empty_subset <- data.frame()
  test_random_subset <- data.frame(
    chrom = c(1),
    start = c(34610),
    end = c(36081),
    Gene = c("FAM138F")
  )
  input_size_subset <- data.frame(
    chrom = c(1, 1, 1, 1),
    start = c(11873, 14361, 34610, 69090),
    end = c(14409, 16765, 36081, 70008),
    Gene = c("DDX11L1", "WASH7P", "FAM138F", "OR4F5")
  )
  expected_results <- data.frame(
    chrom = c(1, 1, 1),
    start = c(11873, 14361, 69090),
    end = c(14409, 16765, 70008),
    Gene = c("DDX11L1", "WASH7P", "OR4F5")
  )
  
  initial_results <- remove_prevously_selected_genes(
    input = input_size_subset,
    previous = empty_subset
  )
  test_one_gene <- remove_prevously_selected_genes(
    input = input_size_subset,
    previous = test_random_subset
  )
  test_all_genes <- remove_prevously_selected_genes(
    input = input_size_subset,
    previous = input_size_subset
  )
  expect_identical(input_size_subset, initial_results)
  expect_identical(test_one_gene, expected_results)
  expect_equal(nrow(test_all_genes), 0)
})