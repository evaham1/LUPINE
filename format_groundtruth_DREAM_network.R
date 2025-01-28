### Script to take the tsv file of ground-truth network data from DREAM4 and convert it to matrix format like that outputted by LUPINE

library(tidyverse)

## read in gold-standard network data, table with one gene in each column and then boolean to indicate interaction or not
network_1_data <- read.csv("./data/DREAM4/insilico_size100_1_goldstandard.tsv", sep = "\t")
colnames(network_1_data) <- c("Var1", "Var2", "Value")

head(network_1_data)
# Var1 Var2 Value
# 1   G1   G3     1
# 2   G1   G4     1

# how many variables and therefore possible pairwise interactions are there
length(unique(c(network_1_data$Var1, network_1_data$Var2))) # 100 unique variables
# for 100 variables, the total number of pairwise interactions (undirected) possible excluding self-interactions is:
(100*99) / 2 # 4950
nrow(network_1_data) # 9899 - this is about double, because this network is DIRECTED

## need to turn the directed network into an undirected one
# i.e. if A interacts with B, there is an interaction between A and B

# add interaction names ordered so that we can identify duplicates
network_1_data_sorted <- network_1_data %>%
  mutate(Pair = ifelse(Var1 < Var2, paste(Var1, Var2, sep = "_"), paste(Var2, Var1, sep = "_")))
length(unique(network_1_data_sorted$Pair)) # 4950 - now we have covered all pairwise undirected interactions

network_interactions <- network_1_data_sorted %>%
  group_by(Pair) %>%
  summarize(Interactions = sum(Value), .groups = "drop")
table(network_interactions$Interactions)
# 0    1    2
# 4782  161    7
two_way_interactions <- network_interactions %>%
  filter(Interactions == 2)  # Two-way means both directions (Value == 1) are present
two_way_interactions$Pair
# [1] "G10_G37" "G10_G44" "G44_G62" "G54_G55" "G54_G82" "G65_G66" "G67_G69"


# aggregate by Pair and set 'Value' to 1 if there is any interaction in either direction
network_undirected <- network_1_data_sorted %>%
  group_by(Pair) %>%
  summarize(Value = ifelse(any(Value == 1), 1, 0), .groups = "drop") %>%
  separate(Pair, into = c("Var1", "Var2"), sep = "_")

# checks
nrow(network_undirected) # 4950 pairwise interactions covered
table(network_undirected$Value)
# 0    1
# 4782  168

## Change format from dataframe into a symmetrical matrix with variables on rows and cols in order

# ## making test data to transform to a matrix
# target_mat <- matrix(c(0, 0, 1, 1, 0,
#          0, 0, 0, 0, 1,
#          1, 0, 0, 0, 1,
#          1, 0, 0, 0, 1,
#          0, 1, 1, 1, 0), nrow = 5, ncol = 5, byrow = TRUE)
# diag(target_mat) <- NA
# rownames(target_mat) <- c("G1", "G2", "G3", "G4", "G5")
# colnames(target_mat) <- c("G1", "G2", "G3", "G4", "G5")
# target_mat

# variables <- paste0("G", seq(1, 5))
# network_1 <- matrix(0, nrow = 5, ncol = 5)
# data <- network_1_data[1:5, ]

# Convert tsv of network from DREAM4 into matrix like LUPINE output
variables <- paste0("G", seq(1, 100))
network_1 <- matrix(0, nrow = 100, ncol = 100)
data <- network_undirected

diag(network_1) <- NA
rownames(network_1) <- variables
colnames(network_1) <- variables

# Loop through the input data to update the matrix
for (i in 1:nrow(data)) {
  var1 <- data$Var1[i]
  var2 <- data$Var2[i]
  if (data$Value[i] == 1) {
    # Update both (Var1, Var2) and (Var2, Var1) to 1
    network_1[var1, var2] <- 1
    network_1[var2, var1] <- 1
  }
}

network_1

# check total number of interactions is the same
sum(network_1 == 1, na.rm = TRUE) / 2 # 168
sum(data$Value) # 168

# check matrix is symmetrical
all(is.na(network_1) | network_1 == t(network_1)) # TRUE

# save reformatted network
saveRDS(network_1, "./data/DREAM4/DREAM4_network_1_goldstandard.RDS")
write.csv(network_1, "./data/DREAM4/DREAM4_network_1_goldstandard.csv")

