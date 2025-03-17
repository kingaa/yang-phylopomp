library(tidyverse)
library(phylopomp)
library(pomp)
library(tictoc)
library(ape)
theme_set(theme_bw())

setwd("~/projects/yang-phylopomp")




treeio::read.newick("MERS_274_sCoal.combinedTyped.mcc.newick") -> x

plot(x, show.node.label = TRUE, root.edge = TRUE)

x |>
  treeio::write.tree() |>
  parse_newick() |>
  plot(points=TRUE,ladderize=FALSE)

max_time <- max(node.depth.edgelength(x))
mers_pomp <- parse_newick(treeio::write.tree(x))

pomp_obj <- twospecies_pomp(
  mers_pomp,
  c1 = 1, c2 = 1,
  Beta11 = 18.25, Beta12 = 0, # Species 1 is camel 18.25 from paper # Beta12 is transmission from species 2 TO 1
  Beta21 = 36.5, Beta22 = 0.01 * 365, # Species 2 is human
  gamma1 = 365/14, gamma2 = 365/10, # recovery time for human: 10 days, camel: 14 days (paper)
  b1 = 0, b2 = 0,
  d1 = 0, d2 = 0,
  psi1 = 0.1, psi2 = 0.1, # SAMPLING RATE, ESTIMATE
  omega1 = 0, omega2 = 0,
  S1_0 = 475, I1_0 = 25, R1_0 = 0, # ESTIMATE population, leave it fixed for now
  S2_0 = 5000, I2_0 = 0, R2_0 = 0 # ESTIMATE population
)


# take a better look at interpreting model parameters, time scale matters.
twospecies_params <- data.frame(
  c1 = 1, c2 = 1,
  Beta11 = 18.25, Beta12 = 0, # Species 1 is camel 18.25 from paper # Beta12 is transmission from species 2 TO 1
  Beta21 = 36.5, Beta22 = 0.01 * 365, # Species 2 is human
  gamma1 = 365/14, gamma2 = 365/10, # recovery time for human: 10 days, camel: 14 days (paper)
  b1 = 0, b2 = 0,
  d1 = 0, d2 = 0,
  psi1 = 0.1, psi2 = 0.1, # SAMPLING RATE, ESTIMATE
  omega1 = 0, omega2 = 0,
  S1_0 = 475, I1_0 = 25, R1_0 = 0, # ESTIMATE population, leave it fixed for now
  S2_0 = 5000, I2_0 = 0, R2_0 = 0, # ESTIMATE population
  #S1_0 = 1500000, I1_0 = 0, R1_0 = 0, # ESTIMATE population, leave it fixed for now
  #S2_0 = 63364000, I2_0 = 0, R2_0 = 0 # ESTIMATE population
  t0 = 0,
  time = max_time # max time seen in the tree
)



twospecies_params |>
  select(-time, -Beta22) |>
  expand_grid(
    Beta22 = seq(1, 73, length.out = 36)
  ) |>
  mutate(
    N1 = S1_0 + I1_0 + R1_0,
    N2 = S2_0 + I2_0 + R2_0
  ) |>
  uncount(5) |>  
  collect() -> params


tic("mif2")
np_particles = 5000

{ # mif2
  
  library(iterators)
  library(doFuture)
  plan(multisession)
  
  foreach(
    p = iter(params, "row")
  ) %dofuture% {
    library(phylopomp)
    
    pomp_obj |>
      mif2(
        params = p,
        Np = np_particles,
        Nmif = 50,
        cooling.fraction.50 = 0.5,
        rw.sd = rw_sd(
          Beta11 = 0.005,  
          Beta21 = 0.005,
          #omega1 = 0.005, omega2 = 0.005,
          psi1 = 0.005, psi2 = 0.005
          #S1_0 = ivp(0.2), I1_0 = ivp(0.2), R1_0 = ivp(0.2),
          #S2_0 = ivp(0.2), I2_0 = ivp(0.2), R2_0 = ivp(0.2)
        ),
        partrans = parameter_trans(
          log = c(
            "Beta11", "Beta21"
            #"omega1", "omega2"
          ),
          logit = c("psi1", "psi2"),
          #barycentric = c("S1_0", "I1_0", "R1_0", "S2_0", "I2_0", "R2_0")
        ),
        paramnames = c(
          "Beta11",  "Beta21",
          #"omega1", "omega2",
          "psi1", "psi2"
          #"S1_0", "I1_0", "R1_0",
          #"S2_0", "I2_0", "R2_0"
        )
      )
  } %seed% TRUE |>
    concat()
  
} -> mif2 

saveRDS(mif2, file = "pfs_results_profile_reduced_pop.rds")

num_runs <- 1
# Set up parallel backend

registerDoFuture()
plan(multisession)  # Use all available cores

pfilter_results <- foreach(i = 1:180, .combine = 'rbind', .packages = c('pomp', 'phylopomp')) %:%
  foreach(j = 1:num_runs, .combine = 'rbind') %dopar% {
    pf <- pfilter(mif2[[i]], Np = np_particles)
    c(logLik = logLik(pf), Beta22 = coef(pf)["Beta22"])
  }
toc()
# Load ggplot2 library
library(ggplot2)

# Convert pfilter_results matrix to a data frame
results_df <- data.frame(Beta22 = pfilter_results[, 2], logLik = pfilter_results[, 1])

# Plot using ggplot2
ggplot(data = results_df, aes(x = Beta22, y = logLik)) +
  geom_point(color = "black") +
  labs(
       x = "Beta22",
       y = "logLik") +
  theme_minimal()

ggplot(data = results_df, aes(x = Beta22, y = logLik)) +
  geom_point(color = "black") +
  geom_smooth(method = "loess", color = "green", se = FALSE, span = 1) +
  labs(
    x = "Beta22",
    y = "logLik") +
  theme_minimal()


# Remove rows where logLik is -Inf or NA
filtered_results_df <- results_df[!is.infinite(results_df$logLik) & !is.na(results_df$logLik), ]

# Fit a LOESS model to the filtered data
loess_fit <- loess(logLik ~ Beta22, data = filtered_results_df, span = 1)

# Generate a sequence of Beta22 values for smooth predictions
beta_seq <- seq(min(filtered_results_df$Beta22), max(filtered_results_df$Beta22), length.out = 1000)

# Predict logLik values using the LOESS model
predicted_logLik <- predict(loess_fit, newdata = data.frame(Beta22 = beta_seq))

# Find the maximum predicted logLik and corresponding Beta22
max_idx <- which.max(predicted_logLik)
max_Beta22 <- beta_seq[max_idx]
max_logLik <- predicted_logLik[max_idx]

# Print the results
cat("Maximum logLik:", max_logLik, "at Beta22 =", max_Beta22)

# Add the max point to the plot
ggplot(data = filtered_results_df, aes(x = Beta22, y = logLik)) +
  geom_point(color = "black") +
  geom_smooth(method = "loess", color = "green", se = FALSE, span=1) +
  geom_point(aes(x = max_Beta22, y = max_logLik), color = "red", size = 2) + # Highlight max
  labs(
    x = "Beta22",
    y = "logLik",
    title = paste("LOESS Smooth of Profile likelihood, with max LL at beta22 = ", max_Beta22)
  ) +
  theme_minimal()

write.csv(results_df, "profile_likelihood_data.csv", row.names = FALSE)

# Save the plot
ggsave("profile_likelihood.png", width = 8, height = 6)


