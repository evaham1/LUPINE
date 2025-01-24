### R script to examine how LUPINE function works and make adjustments to generalise to all data types

# load libraries and test data
library(LUPINE)
data("HFHSdata")

## initial testing with LUPINE on single timepoint
debug(LUPINE)
debug(LUPINE_single)

# run LUPINE_single through LUPINE wrapper
net_Normal <- LUPINE(HFHSdata$OTUdata_Normal,
                     is.transformed = FALSE,
                     lib_size = HFHSdata$Lib_Normal,
                     ncomp = 1, single = TRUE,
                     excluded_taxa = HFHSdata$low_Normal_taxa,
                     cutoff = 0.05
)

# check outputs
length(net_Normal) # 4, each timepoints
dim(net_Normal[[1]]) # 212 x 212 variable interaction table
table(net_Normal[[1]])
# 0     1
# 43084  1860

# run LUPINE_single outside LUPINE wrapper
cutoff = 0.05
res <- sapply(1:4, function(d) {
  net <- LUPINE_single(HFHSdata$OTUdata_Normal,
                       timepoint = d,
                       excluded_taxa = HFHSdata$low_Normal_taxa,
                       is.transformed = FALSE,
                       lib_size = HFHSdata$Lib_Normal,
                       ncomp = 1
  )$pvalue
  net <- apply(net<cutoff, c(1, 2), function(x) {
    ifelse(is.na(x), 0, x)
  })
  return(net)
}, simplify = FALSE)

# check outputs
length(res)
dim(res[[1]])
table(res[[1]])
# all same as above


## simplify further - remove the cutoff filtering step for now
res <- sapply(1:4, function(d) {
  net <- LUPINE_single(HFHSdata$OTUdata_Normal,
                     timepoint = d,
                     excluded_taxa = HFHSdata$low_Normal_taxa,
                     is.transformed = FALSE,
                     lib_size = HFHSdata$Lib_Normal,
                     ncomp = 1)
})
dim(res) # 2, 4
dim(res[[1]]) # 212, 212
res[[1]][1:4, 1:4]
# OTU_3       OTU_7      OTU_8       OTU_5
# OTU_3        NA  0.30955285  0.2736252  0.15273790
# OTU_7 0.3095528          NA  0.4122583 -0.01348296
# OTU_8 0.2736252  0.41225832         NA -0.21531472
# OTU_5 0.1527379 -0.01348296 -0.2153147          NA

## run LUPINE_single_timepoint on just the first timepoint data

# extract data for first timepoint
data_first_timepoint <- HFHSdata$OTUdata_Normal[, , 1]
dim(data_first_timepoint) # 23 x 212

# filter taxa in that data
taxa_names <- colnames(data_first_timepoint)[!(colnames(data_first_timepoint) %in% HFHSdata$low_Normal_taxa[[1]])]
data_first_timepoint <- data_first_timepoint[, taxa_names]
dim(data_first_timepoint) # 23, 102

# run modified LUPINE after filtering
net <- LUPINE_single_timepoint(data_first_timepoint,
                       is.transformed = FALSE,
                       lib_size = HFHSdata$Lib_Normal,
                       ncomp = 1)



# check outputs
length(res)
dim(res[[1]])
res[[1]][1:4, 1:4]





