### Try edited LUPINE_single_timepoint on DREAM4 data to see if can predicted network
### running on just first timepoint

library(LUPINE)
library(BiocParallel)
library(igraph)

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
class(network)
class(goldstandard)
dim(network)
dim(goldstandard)

# total number of interactions 346 - way more
sum(network == 1, na.rm = TRUE) / 2 # 346
sum(goldstandard == 1, na.rm = TRUE) / 2 # 168
diag(goldstandard) <- 0
any(is.na(network)) # TRUE - why is this??
network[is.na(network)] <- 0
any(is.na(goldstandard)) # FALSE
any(is.na(network)) # FALSE

### Compare network to ground truth - network distance
dst <- distance_matrix(list(goldstandard, network))
dst # 4.493636, comparing to Saritha's vignette this is not too bad (4-5 for different days of same condition)

### Compare network to ground truth - Mantel test
mantel_res <- MantelTest_matrix(list(goldstandard, network))
mantel_res # 0.78, pretty high comparing to vignette

### Compare network to ground truth - IVI values
net_gold <- graph_from_adjacency_matrix(goldstandard, mode = "undirected")
ivi_gold <- influential::ivi(net_gold)
net_lupine <- graph_from_adjacency_matrix(network, mode = "undirected")
ivi_lupine <- influential::ivi(net_lupine)

compare_ivi <- data.frame(Goldstandard = ivi_gold,
                          LUPINE_output = ivi_lupine)
compare_ivi <- compare_ivi %>%
  mutate(Difference = abs(Goldstandard - LUPINE_output))
mean(compare_ivi$Difference)

tail(ivi_gold[order(ivi_gold)])
tail(ivi_lupine[order(ivi_lupine)])

## try plotting to see difference
plot(net_gold)
plot(net_lupine)

# Venn diagram of interactions?
sum(network == 1, na.rm = TRUE) / 2 # 346
sum(goldstandard == 1, na.rm = TRUE) / 2 # 168

extract_interactions <- function(mat) {
  interactions <- which(mat == 1, arr.ind = TRUE)
  data.frame(
    Var1 = rownames(mat)[interactions[, 1]],
    Var2 = colnames(mat)[interactions[, 2]]
  )
}

LUPINE_interactions <- extract_interactions(network)
Gold_interactions <- extract_interactions(goldstandard)
nrow(LUPINE_interactions)
nrow(Gold_interactions)

LUPINE_interactions <- extract_interactions(network)
nrow(LUPINE_interactions) # 692
LUPINE_interactions <- LUPINE_interactions %>%
  rowwise() %>%
  mutate(
    OrderedVar1 = min(Var1, Var2), # Assign smaller name to OrderedVar1
    OrderedVar2 = max(Var1, Var2)  # Assign larger name to OrderedVar2
  ) %>%
  ungroup() %>%
  select(OrderedVar1, OrderedVar2) %>% # Keep only ordered pairs
  distinct() %>%
  mutate(Interaction_name = paste0(OrderedVar1, "-", OrderedVar2))
nrow(LUPINE_interactions) # 346

Gold_interactions <- extract_interactions(goldstandard)
nrow(Gold_interactions) # 336
Gold_interactions <- Gold_interactions %>%
  rowwise() %>%
  mutate(
    OrderedVar1 = min(Var1, Var2), # Assign smaller name to OrderedVar1
    OrderedVar2 = max(Var1, Var2)  # Assign larger name to OrderedVar2
  ) %>%
  ungroup() %>%
  select(OrderedVar1, OrderedVar2) %>% # Keep only ordered pairs
  distinct() %>%
  mutate(Interaction_name = paste0(OrderedVar1, "-", OrderedVar2))
nrow(Gold_interactions) # 168

intersect(LUPINE_interactions$Interaction_name, Gold_interactions$Interaction_name)

library(VennDiagram)

# Create a Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(LUPINE_interactions$Interaction_name),
  area2 = length(Gold_interactions$Interaction_name),
  cross.area = length(intersect(LUPINE_interactions$Interaction_name, Gold_interactions$Interaction_name)),
  category = c("LUPINE predicted network", "Gold standard network"),
  fill = c("blue", "green"),
  alpha = 0.5
)





