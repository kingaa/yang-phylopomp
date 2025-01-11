library(tidyverse)
library(phylopomp)
library(pomp)
library(tictoc)
library(ape)
theme_set(theme_bw())

#setwd('/home/pyang/projects/yang-phylopomp')




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

#TODO: try estimating a few, leave parameters fixed for now..
pomp_obj <- twospecies_pomp(
  mers_pomp,
  c1 = 0, c2 = 0,
  Beta11 = 4 * 365, Beta12 = 0, # Species 1 is camel # Beta12 is transmission from species 2 TO 1
  Beta21 = 1 * 365, Beta22 = 0.5 * 365, # Species 2 is human
  gamma1 = 1 * 365, gamma2 = 1 * 365,
  b1 = 0, b2 = 0,
  d1 = 0, d2 = 0,
  psi1 = 0.5, psi2 = 0.5, # SAMPLING RATE, ESTIMATE
  omega1 = 0.5 * 365, omega2 = 0.5 * 365,
  S1_0 = 1500, I1_0 = 750, R1_0 = 0, # ESTIMATE population, leave it fixed for now
  S2_0 = 63364, I2_0 = 0, R2_0 = 0, # ESTIMATE population
  #S1_0 = 1500000, I1_0 = 750000, R1_0 = 0, # ESTIMATE population, leave it fixed for now
  #S2_0 = 63364000, I2_0 = 0, R2_0 = 0 # ESTIMATE population
)

# Modify the tree to include root not at time 0, have some time before !!! IMPORTANT
# ballpark 270 infections - > more than model expects 
# Exhausting S pool -> Increase to be a lot larger  ~ 300 each or so

# Can get ballpark estimate for parameters 
# Work out number of infections you see -> adjust to have a similar amount in the parameters. i.e r_0 = 4, 

# Suppose for this data -> Case detected 
# Try increasing I


# take a better look at interpreting model parameters, time scale matters.
twospecies_params <- data.frame(
  c1 = 0, c2 = 0, # NOT ESTIMATING
  Beta11 = 4 * 365, Beta12 = 0, # Species 1 is camel, probably ESTIMATE
  Beta21 = 1 * 365, Beta22 = 0.5 * 365, # Species 2 is human, probably ESTIMATE
  gamma1 = 1 * 365, gamma2 = 1 * 365, # recovery rate, ESTIMATE
  b1 = 0, b2 = 0,
  d1 = 0, d2 = 0,
  psi1 = 0.5, psi2 = 0.5, # SAMPLING RATE, ESTIMATE
  omega1 = 0.5 * 365, omega2 = 0.5 * 365, # Waning immunity, ESIMATE
  S1_0 = 1500, I1_0 = 750, R1_0 = 0, # ESTIMATE population, leave it fixed for now
  S2_0 = 63364, I2_0 = 0, R2_0 = 0, # ESTIMATE population
  #S1_0 = 1500000, I1_0 = 0, R1_0 = 0, # ESTIMATE population, leave it fixed for now
  #S2_0 = 63364000, I2_0 = 0, R2_0 = 0 # ESTIMATE population
  t0 = 0,
  time = max_time # max time seen in the tree
)




twospecies_params |>
  select(-time) |>
  expand_grid(
    rep = seq_len(36),
  ) |>
  mutate(
    N1 = S1_0 + I1_0 + R1_0,
    N2 = S2_0 + I2_0 + R2_0
  ) |>
  collect() -> params


if (TRUE) { # pfilter
  twospecies_params |>
    select(-time, -Beta22) |>
    expand_grid(
      Beta22 = seq(1, 730, length.out = 36)
    ) |>
    mutate(
      N1 = S1_0 + I1_0 + R1_0,
      N2 = S2_0 + I2_0 + R2_0
    ) |>
    collect() -> params
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
        pfilter(params = p, Np = 100)
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
  
  saveRDS(pfs, file = "pfilter_run.rds")
}


# 100 particles: 24 sec
# 200 particles: 41 sec
# 500 particles: 89 sec
# 5000 particles: 930 sec

# Get a sense of monte carlo error 

# roughly linear

# 

#attr(pfs,"system.time")

if (TRUE) {
  tic('mif2') # Start timing
  
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
          Np = 500,
          Nmif = 6,
          cooling.fraction.50 = 0.25,
          rw.sd = rw_sd(
            Beta11 = 0.02, Beta12 = 0.02,
            Beta21 = 0.02, Beta22 = 0.02,
            gamma1 = 0.02, gamma2 = 0.02,
            omega1 = 0.02, omega2 = 0.02,
            psi1 = 0.02, psi2 = 0.02
            #S1_0 = ivp(0.2), I1_0 = ivp(0.2), R1_0 = ivp(0.2),
            #S2_0 = ivp(0.2), I2_0 = ivp(0.2), R2_0 = ivp(0.2)
          ),
          partrans = parameter_trans(
            log = c(
              "Beta11", "Beta12", "Beta21", "Beta22",
              "gamma1", "gamma2",
              "omega1", "omega2"
            ),
            logit = c("psi1", "psi2"),
            #barycentric = c("S1_0", "I1_0", "R1_0", "S2_0", "I2_0", "R2_0")
          ),
          paramnames = c(
            "Beta11", "Beta12", "Beta21", "Beta22",
            "gamma1", "gamma2",
            "omega1", "omega2",
            "psi1", "psi2"
            #"S1_0", "I1_0", "R1_0",
            #"S2_0", "I2_0", "R2_0"
          )
        )
    } %seed% TRUE |>
      concat()
    
  } -> mif2 # End mif2 block

  toc() # End timing
  
  saveRDS(mif2, file = "pfs_results_total_pop.rds")
  #saveRDS(pfs2, file = "pfilterList_object.rds")
  
  
  pfs_mif2 <- readRDS(file = "pfs_results_total_pop.rds")
  
  #pfs_particle <- readRDS(file = "pfilterList_object.rds")
  
  #plot(pfs_particle)
  
  pfs_mif2|>
    traces() |>
    melt() |>
    filter(name %in% c("Beta11","Beta12","Beta21","Beta22","gamma1","gamma2", "omega1", "omega2", "psi1", "psi2", "loglik")) |>
    ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
    geom_line()+
    guides(color="none")+
    facet_wrap(~name,scales="free_y")
}
# keeping np = 5000,
# nmif = 1
# nmif = 2
# nmif = 5

# should reduce is to around 
# Reduce the search a bit, and then 
# keep params mutliple of 36








