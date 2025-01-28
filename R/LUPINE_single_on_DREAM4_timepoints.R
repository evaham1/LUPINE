### Try edited LUPINE_single_timepoint on DREAM4 data to see if can predicted network
### running on each timepoint to see which are more similar

library(LUPINE)
library(BiocParallel)
library(igraph)
library(pheatmap)

### Read in ground-truth network
goldstandard <- readRDS("~/dev/repos/LUPINE/data/DREAM4/DREAM4_network_1_goldstandard.RDS")
diag(goldstandard) <- 0
sum(is.na(goldstandard)) # 0

### Read in in silico gex data
network_1_data <- read.csv("./data/DREAM4/insilico_size100_1_timeseries.tsv", sep = "\t")
network_1_data <- network_1_data %>%
  mutate(Sample = rep(1:10, each = 21)) %>%
  select(Sample, everything())
cutoff = 0.05

# Run one timepoint 1
network_1_data_timepoint <- network_1_data %>%
  filter(Time == 0) %>%
  select(-c(Sample, Time))
res <- LUPINE_single_timepoint(network_1_data_timepoint, ncomp = 1)
network_t1 <- apply(res[[2]], c(1, 2), function(x) {
  ifelse(is.na(x), NA, ifelse(x > cutoff, 0, 1))
})
diag(network_t1) <- 0
sum(is.na(network_t1)) # 0

# Run one timepoint 5
network_1_data_timepoint <- network_1_data %>%
  filter(Time == 200) %>%
  select(-c(Sample, Time))
res <- LUPINE_single_timepoint(network_1_data_timepoint, ncomp = 1)
network_t5 <- apply(res[[2]], c(1, 2), function(x) {
  ifelse(is.na(x), NA, ifelse(x > cutoff, 0, 1))
})
diag(network_t5) <- 0
sum(is.na(network_t5)) # 0

# Run one timepoint 10
network_1_data_timepoint <- network_1_data %>%
  filter(Time == 450) %>%
  select(-c(Sample, Time))
res <- LUPINE_single_timepoint(network_1_data_timepoint, ncomp = 1)
network_t10 <- apply(res[[2]], c(1, 2), function(x) {
  ifelse(is.na(x), NA, ifelse(x > cutoff, 0, 1))
})
diag(network_t10) <- 0
sum(is.na(network_t10)) # 0

# Run one timepoint 15
network_1_data_timepoint <- network_1_data %>%
  filter(Time == 700) %>%
  select(-c(Sample, Time))
res <- LUPINE_single_timepoint(network_1_data_timepoint, ncomp = 1)
network_t15 <- apply(res[[2]], c(1, 2), function(x) {
  ifelse(is.na(x), NA, ifelse(x > cutoff, 0, 1))
})
diag(network_t15) <- 0
sum(is.na(network_t15)) # 0

# Run one timepoint 21
network_1_data_timepoint <- network_1_data %>%
  filter(Time == 1000) %>%
  select(-c(Sample, Time))
res <- LUPINE_single_timepoint(network_1_data_timepoint, ncomp = 1)
network_t21 <- apply(res[[2]], c(1, 2), function(x) {
  ifelse(is.na(x), NA, ifelse(x > cutoff, 0, 1))
})
diag(network_t21) <- 0
sum(is.na(network_t21)) # 0


### Compare networks - network distance
dst <- distance_matrix(list(goldstandard, network_t1, network_t5, network_t10, network_t15, network_t21))
rownames(dst) <- c("Goldstandard", "T1", "T5", "T10", "T15", "T21")
colnames(dst) <- c("Goldstandard", "T1", "T5", "T10", "T15", "T21")
pheatmap(dst, cluster_rows = FALSE, cluster_cols = FALSE)

dst <- distance_matrix(list(network_t1, network_t5, network_t10, network_t15, network_t21))
rownames(dst) <- c("T1", "T5", "T10", "T15", "T21")
colnames(dst) <- c("T1", "T5", "T10", "T15", "T21")
pheatmap(dst, cluster_rows = FALSE, cluster_cols = FALSE)

### Compare networks - Mantel test
mantel_res <- MantelTest_matrix(list(goldstandard, network_t1, network_t5, network_t10, network_t15, network_t21))
rownames(mantel_res) <- c("Goldstandard", "T1", "T5", "T10", "T15", "T21")
colnames(mantel_res) <- c("Goldstandard", "T1", "T5", "T10", "T15", "T21")
pheatmap(mantel_res, cluster_rows = FALSE, cluster_cols = FALSE)

mantel_res <- MantelTest_matrix(list(network_t1, network_t5, network_t10, network_t15, network_t21))
rownames(mantel_res) <- c("T1", "T5", "T10", "T15", "T21")
colnames(mantel_res) <- c("T1", "T5", "T10", "T15", "T21")
pheatmap(mantel_res, cluster_rows = FALSE, cluster_cols = FALSE)


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

LUPINE_interactions <- extract_interactions(network) %>%
  mutate(Pair = ifelse(Var1 < Var2, paste(Var1, Var2, sep = "_"), paste(Var2, Var1, sep = "_"))) %>%
  distinct(Pair)

Gold_interactions <- extract_interactions(goldstandard) %>%
  mutate(Pair = ifelse(Var1 < Var2, paste(Var1, Var2, sep = "_"), paste(Var2, Var1, sep = "_"))) %>%
  distinct(Pair)

nrow(LUPINE_interactions)
nrow(Gold_interactions)
intersect_interactions <- length(intersect(LUPINE_interactions$Pair, Gold_interactions$Pair)) # 11

library(VennDiagram)

# Create a Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = nrow(LUPINE_interactions),
  area2 = nrow(Gold_interactions),
  cross.area = intersect_interactions,
  category = c("LUPINE predicted network", "Gold standard network"),
  fill = c("blue", "green"),
  alpha = 0.5
)





