library(tidyverse)
library(phylopomp)
library(pomp)
theme_set(theme_bw())
#set.seed(353691151)
param_df <- data.frame(sigma = numeric(), beta = numeric(), gamma = numeric(), stringsAsFactors = FALSE)

for (k in 1:10) {
  seirs_params <- data.frame(
    Beta=3,sigma=1,gamma=0.5,psi=0.02,omega=0.08,
    S0=200,E0=10,I0=5,R0=50,
    time=100
  )
  
  #seirs_params <- data.frame(Beta=3,sigma=1,gamma=0.5,psi=0.02,omega=0.08,S0=2000,E0=100,I0=50,R0=50,time=100)
  
  
  seirs_params |>
    with(
      runSEIR(
        Beta=Beta,sigma=sigma,gamma=gamma,psi=psi,omega=omega,
        S0=S0,E0=E0,I0=I0,R0=R0,
        time=time
      )
    )-> seirs_tree
  
  plot_grid(
    seirs_tree |>
      plot(points=FALSE,palette="#000000")+
      labs(x="time"),
    seirs_tree |>
      lineages(obscure=TRUE,prune=TRUE) |>
      plot(),
    ncol=1,
    align="v",
    axis="btlr"
  )
  
  ## dev.off()
  
  seirs_params |>
    with(
      seirs_tree |>
        seirs_pomp(
          Beta=Beta,sigma=sigma,gamma=gamma,psi=psi,omega=omega,
          S0=S0,E0=E0,I0=I0,R0=R0
        )
    ) -> po
  
  seirs_params |>
    select(-time,-sigma, -Beta, -gamma) |>
    expand_grid(
      sigma=seq(0.5,1.5,by=0.5),
      Beta = seq(1, 6, by=2.5),
      gamma = seq(.2, .8, by =0.2),
      rep=seq_len(1)
    ) |>
    mutate(
      N=S0+E0+I0+R0
    ) |>
    collect() -> params
  
  
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
      po |>
        mif2(
          params=p,
          Np=5000,
          Nmif=50,
          cooling.fraction.50=0.25, # can try decreasing to lower noise
          rw.sd=rw_sd(Beta = 0.02, gamma = 0.02, sigma = 0.02),
          partrans=parameter_trans(log=c("Beta","gamma", "sigma")),
          paramnames=c("Beta","gamma", "sigma")
        )
    } %seed% TRUE |>
      concat()
  } -> pfs
  
  attr(pfs,"system.time")
  
  pfs |>
    traces() |>
    melt() |>
    filter(name %in% c("Beta","sigma","gamma","loglik")) |>
    ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
    geom_line()+
    guides(color="none")+
    facet_wrap(~name,scales="free_y")
  
  
  max_ll = pfs[[1]]@loglik
  max_ll_idx = 1
  
  for (i in 1:36){
    if(pfs[[i]]@loglik > max_ll) {
      max_ll = pfs[[i]]@loglik
      max_ll_idx = i
    }
  }
  
  best_sigma = pfs[[max_ll_idx]]@params[7]
  best_Beta = pfs[[max_ll_idx]]@params[8]
  best_gamma = pfs[[max_ll_idx]]@params[9]
  
  param_df <- rbind(param_df, data.frame(sigma= best_sigma, beta = best_Beta, gamma = best_gamma))
  
  
}

#dev.off()
