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

# run LUPINE_single_timepoint in exact same way
cutoff = 0.05
res <- sapply(1:4, function(d) {
  net <- LUPINE_single_timepoint(HFHSdata$OTUdata_Normal,
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







