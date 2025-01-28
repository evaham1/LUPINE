### Try edited LUPINE_single_timepoint on DREAM4 data to see if can predicted network
### running on each timepoint to see which are more similar

library(LUPINE)
library(BiocParallel)
library(igraph)
library(pheatmap)
library(ggplot2)

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

fit <- data.frame(cmdscale(dst, k = 3))
p2 <- ggplot(fit, aes(x = X1, y = X2)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(fit)),
            parse = TRUE, hjust = 0.5, vjust = 1.5,
            show.legend = FALSE)
p2

# no goldstandard
# dst <- distance_matrix(list(network_t1, network_t5, network_t10, network_t15, network_t21))
# rownames(dst) <- c("T1", "T5", "T10", "T15", "T21")
# colnames(dst) <- c("T1", "T5", "T10", "T15", "T21")
# pheatmap(dst, cluster_rows = FALSE, cluster_cols = FALSE)

### Compare networks - Mantel test
mantel_res <- MantelTest_matrix(list(goldstandard, network_t1, network_t5, network_t10, network_t15, network_t21))
rownames(mantel_res) <- c("Goldstandard", "T1", "T5", "T10", "T15", "T21")
colnames(mantel_res) <- c("Goldstandard", "T1", "T5", "T10", "T15", "T21")
pheatmap(mantel_res, cluster_rows = FALSE, cluster_cols = FALSE)

# no goldstandard
mantel_res <- MantelTest_matrix(list(network_t1, network_t5, network_t10, network_t15, network_t21))
rownames(mantel_res) <- c("T1", "T5", "T10", "T15", "T21")
colnames(mantel_res) <- c("T1", "T5", "T10", "T15", "T21")
pheatmap(mantel_res, cluster_rows = FALSE, cluster_cols = FALSE)


### Compare network to ground truth - IVI values
net_gold <- graph_from_adjacency_matrix(goldstandard, mode = "undirected")
ivi_gold <- influential::ivi(net_gold)
net_t1 <- graph_from_adjacency_matrix(network_t1, mode = "undirected")
ivi_t1 <- influential::ivi(net_t1)
net_t5 <- graph_from_adjacency_matrix(network_t5, mode = "undirected")
ivi_t5 <- influential::ivi(net_t5)
net_t10 <- graph_from_adjacency_matrix(network_t10, mode = "undirected")
ivi_t10 <- influential::ivi(net_t10)
net_t15 <- graph_from_adjacency_matrix(network_t15, mode = "undirected")
ivi_t15 <- influential::ivi(net_t15)
net_t21 <- graph_from_adjacency_matrix(network_t21, mode = "undirected")
ivi_t21 <- influential::ivi(net_t21)
IVI_comb <- rbind(ivi_gold, ivi_t1, ivi_t5, ivi_t10, ivi_t15, ivi_t21)
pca_ivi <- mixOmics::pca(IVI_comb)
fit1 <- data.frame(pca_ivi$variates$X) %>% cbind(name = pca_ivi$names$sample)

p1 <- ggplot(fit1, aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text(aes(label = rownames(fit1)),
            parse = TRUE, hjust = 0.5, vjust = 1.5,
            show.legend = FALSE)
p1

