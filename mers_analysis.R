library(tidyverse)
library(phylopomp)
library(pomp)
library(tictoc)
library(ape)
library(foreach)
library(doFuture)
plan(multisession)

theme_set(theme_bw())

setwd("~/projects/yang-phylopomp")
pfs_mif2 <- readRDS(file = "mif2_reduced_pop_results16.rds")


invalid_groups <- pfs_mif2 |>
  traces() |>
  melt() |>
  filter(name == "loglik", value == -Inf) |>
  pull(.L1) |>  # Extract the .L1 values
  unique()  # Keep only unique values

# Filter out entire trajectories belonging to these groups
pfs_mif2 |>
  traces() |>
  melt() |>
  filter(name %in% c("Beta11","Beta21","Beta22","psi1", "psi2", "loglik")) |>
  filter(!.L1 %in% invalid_groups) |>  # Remove groups that had -Inf in loglik
  ggplot(aes(x=iteration, y=value, group=.L1, color=factor(.L1))) +
  geom_line() +
  guides(color="none") +
  facet_wrap(~name, scales="free_y")

pfs_mif2 |>
  traces() |>
  melt() |>
  filter(name == "loglik", iteration >= 50, iteration <= 200) |>
  ggplot(aes(x = iteration, y = value, group = .L1, color = factor(.L1))) +
  geom_line() +
  guides(color = "none") 
  
treeio::read.newick("MERS_274_sCoal.combinedTyped.mcc.newick") -> x

# Suppose pfs_mif2 is your mif2List object
loglik_values <- sapply(pfs_mif2, logLik)  # Extract log-likelihoods

# Get the indices of the top 5 log-likelihoods
top5_indices <- order(loglik_values, decreasing = TRUE)[1:5]

# Extract the top 5 mif2 objects
top5_mif2 <- pfs_mif2[top5_indices]


# Set up parallel backend
registerDoFuture()
plan(multisession)  # Use all available cores

# Number of filtering runs per mif2 object
num_runs <- 10
np_particles <- 5000

tic("pfilter: ")
# Run pfilter in parallel
pfilter_results <- foreach(i = 1:5, .combine = 'c', .packages = c('pomp', 'phylopomp')) %:%
  foreach(j = 1:num_runs, .combine = 'c') %dopar% {
    logLik(pfilter(top5_mif2[[i]], Np = np_particles))
  }
toc()
  # Convert results to a data frame
pfilter_df <- data.frame(
  Model = rep(1:5, each = num_runs),
  Run = rep(1:num_runs, times = 5),
  LogLik = pfilter_results
)

# Calculate summary statistics for each model
summary_stats <- aggregate(LogLik ~ Model, data = pfilter_df, 
                           FUN = function(x) c(mean = mean(x), sd = sd(x), range = range(x)))

# Display the results
print(pfilter_df)
print(summary_stats)

# Close parallel backend
plan(sequential)

