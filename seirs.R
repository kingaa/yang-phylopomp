library(tidyverse)
library(phylopomp)
library(pomp)
theme_set(theme_bw())
#set.seed(353691151)

seirs_params <- data.frame(
  Beta=3,sigma=1,gamma=0.5,psi=0.02,omega=0.08,
  S0=200,E0=10,I0=5,R0=50,
  time=100
)

#seirs_params <- data.frame(Beta=3,sigma=1,gamma=0.5,psi=0.02,omega=0.08,S0=2000,E0=100,I0=50,R0=50,time=100)

bake(
  file="seirs1.rds",
  dependson=seirs_params,
  #seed=509673338,
  ## seed=1903942427,
  ## seed=874064367,
  ## seed=857598803,
  ## seed=1293274224,
  ## seed=925687277,
  seirs_params |>
    with(
      runSEIR(
        Beta=Beta,sigma=sigma,gamma=gamma,psi=psi,omega=omega,
        S0=S0,E0=E0,I0=I0,R0=R0,
        time=time
      )
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
        Beta=6,sigma=sigma,gamma=gamma,psi=psi,omega=omega,
        S0=S0,E0=E0,I0=I0,R0=R0
      )
  ) -> po

seirs_params |>
  select(-time,-sigma, -Beta) |>
  expand_grid(
    sigma=seq(0.5,1.5,by=0.3),
    Beta = seq(1, 6, by=0.5),
    rep=seq_len(1)
  ) |>
  mutate(
    gamma=2,
    N=S0+E0+I0+R0
  ) |>
  collect() -> params

bake(
  file="seirs2_r1.rds",
  dependson=list(params,seirs_params,seirs_tree),
  seed=751601556,
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
          Nmif=1,
          cooling.fraction.50=0.25, # can try decreasing to lower noise
          rw.sd=rw_sd(Beta = 0.02, gamma = 0.02, sigma = 0.02),
          partrans=parameter_trans(log=c("Beta","gamma", "sigma")),
          paramnames=c("Beta","gamma", "sigma")
        )
    } %seed% TRUE |>
      concat()
  }
) -> pfs

attr(pfs,"system.time")

pfs |>
  traces() |>
  melt() |>
  filter(name %in% c("Beta","sigma","gamma","loglik")) |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")

x = pfs[[1]]
pfs2 <- continue(x, Np=5000,
         Nmif=10,                        
         cooling.fraction.50=0.25,        
         rw.sd=rw_sd(Beta = 2, gamma = 2), 
         partrans=parameter_trans(log=c("Beta","gamma")),
         paramnames=c("Beta","gamma"))
pfs |>
  traces() |>
  melt() |>
  filter(name %in% c("Beta","sigma","gamma","loglik")) |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")


bake(
  file="seirs2_run1.rds",
  dependson=list(pfs, params, seirs_tree), # Update to depend on previous results
  seed=751601556,
  {
    foreach (
      p=iter(params,"row")
    ) %dofuture% {
      library(phylopomp)
      po |>
        continue(
          params=p,                        
          Np=5000,
          Nmif=10,                        
          cooling.fraction.50=0.25,        
          rw.sd=rw_sd(Beta = 0.002, gamma = 0.002), 
          partrans=parameter_trans(log=c("Beta","gamma")),
          paramnames=c("Beta","gamma")
        )
    } %seed% TRUE |>
      concat()
  }
) -> pfs_1

attr(pfs,"system.time")

## IF2 trace plots
pfs_1|>
  traces() |>
  melt() |>
  filter(name %in% c("Beta","sigma","gamma","loglik")) |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")

## make a profile likelihood plot
bind_cols(
  pfs |>
    coef() |>
    melt() |>
    pivot_wider(),
  logLik=logLik(pfs)
) |>
  ggplot(aes(x=sigma,y=logLik))+
  geom_point()+
  theme_bw()

plot(pfs)

foreach(mf=pfs,.combine=rbind,
  .options.future=list(seed=900242057)
) %dofuture% {
  evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
  ll <- logmeanexp(evals,se=TRUE)
  mf |> coef() |> bind_rows() |>
    bind_cols(loglik=ll[1],loglik.se=ll[2])
} -> results


## seirs_tree |>
##   curtail(time=1.28) |>
##   plot(prune=FALSE,points=TRUE,obscure=FALSE)

## pfs |> cond_logLik()
## pfs |> eff_sample_size()

left_join(
  pfs |> coef() |> melt() |> pivot_wider(),
  pfs |> logLik() |> melt() |> rename(logLik=value),
  by=c(".id"="name")
) -> params

params |>
  ##  filter(is.finite(logLik)) |>
  with(
    mcap(logLik,sigma,span=0.5)
  ) -> mcap

plot_grid(
  A=seirs_tree |>
    plot(points=TRUE,palette="#000000")+
    labs(x="time"),
  B=params |>
    ggplot()+
    geom_point(aes(x=sigma,y=logLik))+
    geom_line(data=mcap$fit,aes(x=parameter,y=smoothed),color="blue")+
    geom_vline(xintercept=seirs_params$sigma,color="red")+
    geom_vline(xintercept=mcap$ci,linetype=2)+
    geom_hline(
      yintercept=with(mcap,max(fit$smoothed)-c(0,delta)),
      linetype=2
    )+
    labs(
      color=character(0),
      y="log likelihood",
      x=expression(sigma)
    )+
    theme_classic(),
  labels="AUTO",
  ncol=1,
  rel_heights=c(1,1)
)

dev.off()
