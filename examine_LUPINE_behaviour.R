### R script to examine how LUPINE function works and make adjustments to generalise to all data types

# load libraries and test data
library(LUPINE)
data("HFHSdata")

## initial testing with LUPINE on single timepoint
debug(LUPINE)
debug(LUPINE_single)

# try what is the difference for is.transformed TRUE/FALSE
data <- HFHSdata$OTUdata_Normal
net_Normal <- LUPINE(HFHSdata$OTUdata_Normal,
                     is.transformed = FALSE,
                     lib_size = HFHSdata$Lib_Normal,
                     ncomp = 1, single = TRUE,
                     excluded_taxa = HFHSdata$low_Normal_taxa,
                     cutoff = 0.05
)

# quick visualisation to check for consistency
Day0 <- HFHSdata$OTUdata_Normal[, , 1]
taxa_info <- HFHSdata$filtered_taxonomy[colnames(Day0), ]
taxa_info$X5 <- factor(taxa_info$X5)
col_vec <- rep(
  c(
    "green", "gray", "darkgreen", "darkred", "firebrick2", "pink",
    "tomato", "orange", "blue", "purple", "hotpink", "lightblue"
  ),
  summary(taxa_info$X5)
)
p1 <- netPlot_HFHS(net_Normal[[3]], col_vec, "Normal Day7")
p1

net_Normal <- LUPINE(HFHSdata$OTUdata_Normal,
                     is.transformed = FALSE,
                     lib_size = HFHSdata$Lib_Normal,
                     ncomp = 1, single = TRUE,
                     excluded_taxa = HFHSdata$low_Normal_taxa,
                     cutoff = 0.05
)

# quick visualisation to check for consistency
Day0 <- HFHSdata$OTUdata_Normal[, , 1]
taxa_info <- HFHSdata$filtered_taxonomy[colnames(Day0), ]
taxa_info$X5 <- factor(taxa_info$X5)
col_vec <- rep(
  c(
    "green", "gray", "darkgreen", "darkred", "firebrick2", "pink",
    "tomato", "orange", "blue", "purple", "hotpink", "lightblue"
  ),
  summary(taxa_info$X5)
)
p1 <- netPlot_HFHS(net_Normal[[3]], col_vec, "Normal Day7")
p1
