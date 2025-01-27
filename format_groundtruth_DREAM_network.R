### Script to take the tsv file of ground-truth network data from DREAM4 and convert it to matrix format like that outputted by LUPINE

library(tidyverse)

# read in gold-standard network data, boolean to indicate interaction or not
network_1_data <- read.csv("./data/DREAM4/insilico_size100_1_goldstandard.tsv", sep = "\t")
colnames(network_1_data) <- c("Var1", "Var2", "Value") # 9899 rows - does not include every possible interaction

# filter to only include rows displaying an interaction (i.e. Value = 1)
network_1_data <- network_1_data %>% filter(Value == 1) # only keep rows which have an interaction, 175 rows left

# filter to remove duplicate rows (order of Var1 and Var2 doesn't matter)
network_1_data[c("Var1", "Var2")] <- t(apply(network_1_data[c("Var1", "Var2")], 1, sort))
duplicates <- network_1_data[duplicated(network_1_data[c("Var1", "Var2")]) | duplicated(network_1_data[c("Var1", "Var2")], fromLast = TRUE), ]
nrow(duplicates) # 14 duplicates i.e. 7 times 2
network_1_data <- network_1_data %>%
  filter(!(Var1 %in% duplicates$Var1 & Var2 %in% duplicates$Var2) |
           (duplicated(network_1_data[c("Var1", "Var2")]) == FALSE))
nrow(network_1_data) # 168 true interactions left in data

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
data <- network_1_data

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
sum(data$Value) # 175

# save reformatted network
saveRDS(network_1, "./data/DREAM4/DREAM4_network_1_goldstandard.RDS")
write.csv(network_1, "./data/DREAM4/DREAM4_network_1_goldstandard.csv")

