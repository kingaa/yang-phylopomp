library(tidyverse)
library(phylopomp)
library(pomp)
theme_set(theme_bw())

setwd('/home/pyang/projects/yang-phylopomp')



treeio::read.nexus("MERS_274_sCoal.combinedTyped.mcc.tree") -> x


x$root.edge <- 1 # edit root edge
x$tip.label[] <- "b_"
x$node.label[] <- "g_"

plot(x, show.node.label = TRUE, root.edge = TRUE)

# It seems that write.tree() returns a newick object, and from what i've seen online,
# Newick trees don't support root edges. Not sure how to proceed here...
# Probably need to find a way to edit the raw gpgen in mers_pomp...

# TODO: if more than one thing coming out of a root, insert the appropriate number of roots
# Wrap whole thing in parenthesis
#x |>
#  treeio::write.tree() |>
#  parse_newick() |>
#  curtail(time=2, prune=FALSE) |>
#  diagram(prune=FALSE)

x |>
  treeio::write.tree() |>
  parse_newick() |>
  plot(points=TRUE,ladderize=FALSE)

mers_pomp <- parse_newick(treeio::write.tree(x))
plot(mers_pomp, 5, 0)

#diagram(curtail(mers_pomp,0.5),prune=F,obscure=F)

pomp_obj <- twospecies_pomp(
  mers_pomp,
  c1 = 0, c2 = 0,
  Beta11 = 4, Beta12 = 1,
  Beta21 = 1, Beta22 = 4,
  gamma1 = 1, gamma2 = 1,
  b1 = 0, b2 = 0,
  d1 = 0, d2 = 0,
  psi1 = 1, psi2 = 1,
  omega1 = 0.5, omega2 = 0.5,
  S1_0 = 300, I1_0 = 10, R1_0 = 10,
  S2_0 = 300, I2_0 = 10, R2_0 = 10
)

# Modify the tree to include root not at time 0, have some time before !!! IMPORTANT
# ballpark 270 infections - > more than model expects 
# Exhausting S pool -> Increase to be a lot larger  ~ 300 each or so

# Can get ballpark estimate for parameters 
# Work out number of infections you see -> adjust to have a similar amount in the parameters. i.e r_0 = 4, 

# Suppose for this data -> Case detected 
# Try increasing I


twospecies_params <- data.frame(
  c1 = 0, c2 = 0,
  Beta11 = 4, Beta12 = 1,
  Beta21 = 1, Beta22 = 4,
  gamma1 = 1, gamma2 = 1,
  b1 = 0, b2 = 0,
  d1 = 0, d2 = 0,
  psi1 = 1, psi2 = 1,
  omega1 = 0.5, omega2 = 0.5,
  S1_0 = 300, I1_0 = 10, R1_0 = 10,
  S2_0 = 300, I2_0 = 10, R2_0 = 10,
  t0 = 0,
  time = 5.6512919 # max time seen in the tree
)


twospecies_params |>
  select(-time) |>
  expand_grid(
    rep=seq_len(200)
  ) |>
  mutate(
    
    N1 =S1_0+I1_0+R1_0,
    N2 =S2_0+I2_0+R2_0
  ) |>
  collect() -> params

if (TRUE) { # pfilter
  {
    library(iterators)
    library(doFuture)
    plan(multisession)
    ## cl <- makeClusterMPI(250,autostop=TRUE,verbose=FALSE)
    ## plan(cluster,workers=cl)
    foreach (
      p=iter(params,"row")
    ) %dofuture% {
      library(phylopomp)
      pomp_obj |>
        pfilter(params = p, Np = 1e3)
    } %seed% TRUE |>
      concat()
  } -> pfs
}

{ # mif2
  library(iterators)
  library(doFuture)
  plan(multisession)
  foreach (
    p = iter(params, "row")
  ) %dofuture% {
    library(phylopomp)
    pomp_obj |>
      mif2(
        params = p,
        Np = 5000,
        Nmif = 10,
        cooling.fraction.50 = 0.25,
        rw.sd = rw_sd(
          c1 = 0.01, c2 = 0.01,
          Beta11 = 0.02, Beta12 = 0.02,
          Beta21 = 0.02, Beta22 = 0.02,
          gamma1 = 0.01, gamma2 = 0.01,
          psi1 = 0.01, psi2 = 0.01,
          b1 = 0.01, b2 = 0.01,
          d1 = 0.01, d2 = 0.01,
          omega1 = 0.01, omega2 = 0.01
        ),
        partrans = parameter_trans(
          log = c("c1", "c2",
                  "Beta11", "Beta12", "Beta21", "Beta22",
                  "gamma1", "gamma2",
                  "psi1", "psi2",
                  "b1", "b2",
                  "d1", "d2",
                  "omega1", "omega2")
        ),
        paramnames = c("c1", "c2",
                       "Beta11", "Beta12", "Beta21", "Beta22",
                       "gamma1", "gamma2",
                       "psi1", "psi2",
                       "b1", "b2",
                       "d1", "d2",
                       "omega1", "omega2")
      )
  } %seed% TRUE |>
    concat()
} -> mif2




saveRDS(pfs, file = "pfs_results.rds")
saveRDS(pfs2, file = "pfilterList_object.rds")


pfs_mif2 <- readRDS(file = "pfs_results.rds")

pfs_particle <- readRDS(file = "pfilterList_object.rds")

plot(pfs_particle)







