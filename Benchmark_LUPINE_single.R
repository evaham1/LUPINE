#### Script to compare original LUPINE_single with modified LUPINE_single_timepoint on microbiome data
#### run using the microbiome HFHSdata from vignette, only on day 1

# Sometimes just check that the result you get with LUPINE_single is the same as the original unmodified LUPINE package:
# devtools::install_github("https://github.com/SarithaKodikara/LUPINE")

## Run original LUPINE_single function
res_orig <- LUPINE_single(HFHSdata$OTUdata_Normal,
              day = 1,
              excluded_taxa = HFHSdata$low_Normal_var,
              is.transformed = FALSE,
              lib_size = HFHSdata$Lib_Normal,
              ncomp = 1)
# Estimates
res_orig$Estimate
dim(res_orig[[1]]) # 212, 212
res_orig[[1]][1:4, 1:4] # note the below output is from the original LUPINE installed from github, should match exactly
# OTU_3       OTU_7      OTU_8      OTU_5
# OTU_3         NA  0.04243286  0.2052954  0.0152339
# OTU_7 0.04243286          NA  0.3509201 -0.2096235
# OTU_8 0.20529540  0.35092010         NA -0.1492258
# OTU_5 0.01523390 -0.20962352 -0.1492258         NA

# P-values
res_orig$pvalue
dim(res_orig[[2]]) # 212, 212
res_orig[[2]][1:4, 1:4] # note the below output is from the original LUPINE installed from github, should match exactly
# OTU_3     OTU_7     OTU_8     OTU_5
# OTU_3        NA 0.8475554 0.3473617 0.9449986
# OTU_7 0.8475554        NA 0.1006321 0.3370656
# OTU_8 0.3473617 0.1006321        NA 0.4967753
# OTU_5 0.9449986 0.3370656 0.4967753        NA


## Run modified LUPINE_single_timepoint and make sure get the same results

# extract data for first timepoint
data_first_timepoint <- HFHSdata$OTUdata_Normal[, , 1]
dim(data_first_timepoint) # 23 x 212

# filter var in that data
var_names <- colnames(data_first_timepoint)[!(colnames(data_first_timepoint) %in% HFHSdata$low_Normal_var[[1]])]
data_first_timepoint <- data_first_timepoint[, var_names]
dim(data_first_timepoint) # 23, 102

# run modified LUPINE after filtering
library(BiocParallel)
res_new <- LUPINE_single_timepoint(data_first_timepoint,
                               lib_size = HFHSdata$Lib_Normal[, 1],
                               ncomp = 1)
# Estimates
res_new$Estimate
dim(res_new[[1]]) # 212, 212
res_new[[1]][1:4, 1:4]

# P-values
res_new$pvalue
dim(res_new[[2]]) # 212, 212
res_new[[2]][1:4, 1:4]

## Results look the same, quick check that they are identical:
library(testthat)

test_that("LUPINE_single_timepoint output is identical to LUPINE_single output for day 1 data", {

  res_orig <- LUPINE_single(HFHSdata$OTUdata_Normal,
                            day = 1,
                            excluded_taxa = HFHSdata$low_Normal_var,
                            is.transformed = FALSE,
                            lib_size = HFHSdata$Lib_Normal,
                            ncomp = 1)

  data_first_timepoint <- HFHSdata$OTUdata_Normal[, , 1]
  var_names <- colnames(data_first_timepoint)[!(colnames(data_first_timepoint) %in% HFHSdata$low_Normal_var[[1]])]
  data_first_timepoint <- data_first_timepoint[, var_names]
  res_new <- LUPINE_single_timepoint(data_first_timepoint,
                                     lib_size = HFHSdata$Lib_Normal[, 1],
                                     ncomp = 1)

  expect_identical(res_orig[[1]], res_new[[1]])
  expect_identical(res_orig[[2]], res_new[[2]])
  expect_true(all(dim(res_orig) == dim(res_new)), "Matrix dimensions do not match")
  expect_true(all(rownames(res_orig) == rownames(res_new)), "Row names do not match")
  expect_true(all(colnames(res_orig) == colnames(res_new)), "Column names do not match")
})

test_that("LUPINE_single_timepoint output is identical to LUPINE_single output for day 2 data", {

  res_orig <- LUPINE_single(HFHSdata$OTUdata_Normal,
                            day = 2,
                            excluded_taxa = HFHSdata$low_Normal_var,
                            is.transformed = FALSE,
                            lib_size = HFHSdata$Lib_Normal,
                            ncomp = 1)

  data_first_timepoint <- HFHSdata$OTUdata_Normal[, , 2]
  var_names <- colnames(data_first_timepoint)[!(colnames(data_first_timepoint) %in% HFHSdata$low_Normal_var[[2]])]
  data_first_timepoint <- data_first_timepoint[, var_names]
  res_new <- LUPINE_single_timepoint(data_first_timepoint,
                                     lib_size = HFHSdata$Lib_Normal[, 2],
                                     ncomp = 1)

  expect_identical(res_orig[[1]], res_new[[1]])
  expect_identical(res_orig[[2]], res_new[[2]])
  expect_true(all(dim(res_orig) == dim(res_new)), "Matrix dimensions do not match")
  expect_true(all(rownames(res_orig) == rownames(res_new)), "Row names do not match")
  expect_true(all(colnames(res_orig) == colnames(res_new)), "Column names do not match")
})

### Now test if LUPINE_single_timepoint is any faster than LUPINE_single as removed loop and can run apply in parallel
library(microbenchmark)

# benchmark with default BPPARAM ie not in parallel
benchmark_results <- microbenchmark(
  LUPINE_single(HFHSdata$OTUdata_Normal,
                day = 1,
                excluded_taxa = HFHSdata$low_Normal_var,
                is.transformed = FALSE,
                lib_size = HFHSdata$Lib_Normal,
                ncomp = 1),
  LUPINE_single_timepoint(data_first_timepoint,
                          lib_size = HFHSdata$Lib_Normal[, 1],
                          ncomp = 1),
  times = 10
)

print(benchmark_results)
benchmark_results$expr <- factor(
  benchmark_results$expr,
  levels = unique(benchmark_results$expr),
  labels = c("LUPINE_single", "LUPINE_single_timepoint")
)
boxplot(benchmark_results,
        ylab = "Execution Time (s)") # modified function is faster

# benchmark with two cores
benchmark_results <- microbenchmark(
  LUPINE_single(HFHSdata$OTUdata_Normal,
                day = 1,
                excluded_taxa = HFHSdata$low_Normal_var,
                is.transformed = FALSE,
                lib_size = HFHSdata$Lib_Normal,
                ncomp = 1),
  LUPINE_single_timepoint(data_first_timepoint,
                          lib_size = HFHSdata$Lib_Normal[, 1],
                          ncomp = 1,
                          BPPARAM = SnowParam(workers = 2)),
  times = 10
)

print(benchmark_results)
benchmark_results$expr <- factor(
  benchmark_results$expr,
  levels = unique(benchmark_results$expr),
  labels = c("LUPINE_single", "LUPINE_single_timepoint")
)
boxplot(benchmark_results,
        ylab = "Execution Time (s)") # modified function is even faster when run on 2 cores



