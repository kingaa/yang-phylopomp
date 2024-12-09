library(tidyverse)
library(phylopomp)
library(pomp)
library(tictoc)
library(ape)
theme_set(theme_bw())

setwd('/home/pyang/projects/yang-phylopomp')




treeio::read.newick("MERS_274_sCoal.combinedTyped.mcc.newick") -> x


#x$root.edge <- 1 # edit root edge
#x$tip.label[] <- "b_"
#x$node.label[] <- "g_"


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


# Aaron's hypothesis: length of root is ignored, need to add new root, not just edge
x |>
  treeio::write.tree() |>
  parse_newick() |>
  plot(points=TRUE,ladderize=FALSE)

max_time <- max(node.depth.edgelength(x))
mers_pomp <- parse_newick(treeio::write.tree(x))

#mers_pomp <- curtail(mers_pomp, troot = 2) # cuts the tree, lops off everything before time troot

#mers_pomp |> plot(points=TRUE, ladderize = FALSE) # 

#diagram(curtail(mers_pomp,0.5),prune=F,obscure=F)

# can try curtailing from t = 1 in the range of around t = 2.5
pomp_obj <- twospecies_pomp(
  mers_pomp,
  c1 = 0, c2 = 0,
  Beta11 = 4, Beta12 = 0, # Species 1 is camel # Beta12 is transmission from species 2 TO 1
  Beta21 = 1, Beta22 = 0.5, # Species 2 is human
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
  Beta11 = 4, Beta12 = 0, # Species 1 is camel
  Beta21 = 1, Beta22 = 0.5, # Species 2 is human
  gamma1 = 1, gamma2 = 1,
  b1 = 0, b2 = 0,
  d1 = 0, d2 = 0,
  psi1 = 1, psi2 = 1,
  omega1 = 0.5, omega2 = 0.5,
  S1_0 = 300, I1_0 = 10, R1_0 = 10,
  S2_0 = 300, I2_0 = 10, R2_0 = 10,
  t0 = 0,
  time = max_time # max time seen in the tree
)



twospecies_params |>
  select(-time, -Beta22) |>
  expand_grid(
    rep = seq_len(5),
    Beta22 = seq(0, 4, length.out = 50)
  ) |>
  mutate(
    N1 = S1_0 + I1_0 + R1_0,
    N2 = S2_0 + I2_0 + R2_0
  ) |>
  collect() -> params

if (TRUE) { # pfilter
  tic('pfilter')
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
        pfilter(params = p, Np = 5000)
    } %seed% TRUE |>
      concat()
  } -> pfs
  toc()
  
  left_join(
    pfs |> coef() |> melt() |> pivot_wider(),
    pfs |> logLik() |> melt() |> rename(logLik=value),
    by=c(".id"="name")
  ) -> params
  
  params |>
      filter(is.finite(logLik)) |>
    with(
      mcap(logLik,Beta22,span=0.5)
    ) -> mcap
  
  plot_grid(
    A=mers_pomp |>
      plot(points=TRUE,palette="#000000")+
      labs(x="time"),
    B=params |>
      ggplot()+
      geom_point(aes(x=Beta22,y=logLik))+
      geom_line(data=mcap$fit,aes(x=parameter,y=smoothed),color="blue")+
      geom_vline(xintercept=twospecies_params$Beta22,color="red")+
      geom_vline(xintercept=mcap$ci,linetype=2)+
      geom_hline(
        yintercept=with(mcap,max(fit$smoothed)-c(0,delta)),
        linetype=2
      )+
      labs(
        color=character(0),
        y="log likelihood",
        x=expression(Beta22)
      )+
      theme_classic(),
    labels="AUTO",
    ncol=1,
    rel_heights=c(1,1)
  )
}


# 100 particles: 24 sec
# 200 particles: 41 sec
# 500 particles: 89 sec
# 5000 particles: 930 sec

# roughly linear

attr(pfs,"system.time")

if (TRUE){
  tic('mif2')
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
          Nmif = 1,
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
  toc()
}
# keeping np = 5000,
# nmif = 1
# nmif = 2
# nmif = 5





saveRDS(pfs, file = "pfs_results.rds")
saveRDS(pfs2, file = "pfilterList_object.rds")


pfs_mif2 <- readRDS(file = "pfs_results.rds")

pfs_particle <- readRDS(file = "pfilterList_object.rds")

plot(pfs_particle)







