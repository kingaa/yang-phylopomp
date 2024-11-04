library(tidyverse)
library(phylopomp)
library(pomp)
theme_set(theme_bw())

setwd('/home/pyang/projects/yang-phylopomp')



treeio::read.nexus("MERS_274_sCoal.combinedTyped.mcc.tree") -> x

x$tip.label[] <- "b_"
x$node.label[] <- "g_"

x |> plot()

x |>
  treeio::write.tree() |>
  parse_newick() |>
  plot(points=TRUE,ladderize=FALSE)

mers_pomp <- parse_newick(treeio::write.tree(x))

pomp_obj <- twospecies_pomp(
  mers_pomp, # Parameters taken from example TODO: change
  iota1=0.02,iota2=0.02,
  Beta11=4,Beta12=1,
  Beta21=1,Beta22=4,
  gamma1=1,gamma2=1,
  psi1=1,psi2=1,
  S1_0=100,I1_0=2,R1_0=10,
  S2_0=100,I2_0=2,R2_0=10,
  b1=0.1,b2=0.1,
  d1=0.1,d2=0.1,
  omega1=0.5,omega2=0.5
)


twospecies_params <- data.frame(
  iota1=0.02,iota2=0.02,
  Beta11=4,Beta12=1,
  Beta21=1,Beta22=4,
  gamma1=1,gamma2=1,
  psi1=1,psi2=1,
  S1_0=100,I1_0=2,R1_0=10,
  S2_0=100,I2_0=2,R2_0=10,
  b1=0.1,b2=0.1,
  d1=0.1,d2=0.1,
  omega1=0.5,omega2=0.5,
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

if (FALSE) { # pfilter
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
        pfilter(params = p, Np = 1e4)
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
          iota1 = 0.01, iota2 = 0.01,
          Beta11 = 0.02, Beta12 = 0.02,
          Beta21 = 0.02, Beta22 = 0.02,
          gamma1 = 0.01, gamma2 = 0.01,
          psi1 = 0.01, psi2 = 0.01,
          b1 = 0.01, b2 = 0.01,
          d1 = 0.01, d2 = 0.01,
          omega1 = 0.01, omega2 = 0.01
        ),
        partrans = parameter_trans(
          log = c("iota1", "iota2",
                  "Beta11", "Beta12", "Beta21", "Beta22",
                  "gamma1", "gamma2",
                  "psi1", "psi2",
                  "b1", "b2",
                  "d1", "d2",
                  "omega1", "omega2")
        ),
        paramnames = c("iota1", "iota2",
                       "Beta11", "Beta12", "Beta21", "Beta22",
                       "gamma1", "gamma2",
                       "psi1", "psi2",
                       "b1", "b2",
                       "d1", "d2",
                       "omega1", "omega2")
      )
  } %seed% TRUE |>
    concat()
} -> pfs

attr(pfs, "system.time")
plot(pfs)