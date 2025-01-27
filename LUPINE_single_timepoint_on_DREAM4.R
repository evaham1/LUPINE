### Try edited LUPINE_single_timepoint on DREAM4 data to see if can predicted network
### running on just first timepoint

library(LUPINE)
library(BiocParallel)

## insilico_size100_1_timepoint
# n=10, t=21, p=100 normalised gex data, ground-truth network boolean for each p-p interaction
network_1_data <- read.csv("./data/DREAM4/insilico_size100_1_timeseries.tsv", sep = "\t")
network_1_data <- network_1_data %>%
  mutate(Sample = rep(1:10, each = 21)) %>%
  select(Sample, everything())

# just use first timepoint
network_1_data_timepoint1 <- network_1_data %>%
  filter(Time == 0) %>%
  select(-c(Sample, Time))

# run LUPINE single timepoint
res <- LUPINE_single_timepoint(network_1_data_timepoint1, ncomp = 1)
length(res)
dim(res[[1]]) # 100, 100
res[[1]][1:4, 1:4] # estimates
# G1          G2          G3          G4
# G1          NA -0.01335283 -0.15371531 -0.45057124
# G2 -0.01335283          NA -0.17274200 -0.07775505
# G3 -0.15371531 -0.17274200          NA  0.06738056
# G4 -0.45057124 -0.07775505  0.06738056          NA

dim(res[[2]]) # 100, 100
res[[2]][1:4, 1:4] # p-values

# apply cutoff on p-values to get boolean network
cutoff = 0.05
network <- apply(res[[2]], c(1, 2), function(x) {
  ifelse(is.na(x), NA, ifelse(x > cutoff, 0, 1))
})

### Read in ground-truth network
goldstandard <- readRDS("~/dev/repos/LUPINE/data/DREAM4/DREAM4_network_1_goldstandard.RDS")

# total number of interactions 346 - way more
sum(network == 1, na.rm = TRUE) / 2 # 346
sum(goldstandard == 1, na.rm = TRUE) / 2 # 168

### Compare network to ground truth
p1 <- netPlot_HFHS(net_Normal[[3]], col_vec, "Normal Day7")



