library(tidyverse)
library(phylopomp)
library(pomp)
theme_set(theme_bw())

setwd("~/projects/yang-phylopomp/plots")

param_df <- data.frame(sigma_if2 = numeric(), beta_if2 = numeric(), gamma_if2 = numeric(), if2_ll = numeric(), pf_ll = numeric(), stringsAsFactors = FALSE)

for (k in 2:5) {

seirs_params <- data.frame(
  Beta=3,sigma=1,gamma=0.5,psi=0.02,omega=0.08,
  S0=200,E0=10,I0=5,R0=50,
  time=100
)


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

seirs_params |>
  with(
    seirs_tree |>
      seirs_pomp(
        Beta=Beta,sigma=sigma,gamma=gamma,psi=psi,omega=omega,
        S0=S0,E0=E0,I0=I0,R0=R0
      )
) -> po

seirs_params |>
  select(-time) |>
  expand_grid(
    rep=seq_len(200)
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
          Nmif=10,
          cooling.fraction.50=0.25, # can try decreasing to lower noise
          rw.sd=rw_sd(Beta = 0.02, gamma = 0.02, sigma = 0.02),
          partrans=parameter_trans(log=c("Beta","gamma", "sigma")),
          paramnames=c("Beta","gamma", "sigma")
        )
    } %seed% TRUE |>
      concat()
  } -> mif_run

#attr(mif_run,"system.time")

plot_1 <- mif_run |>
  traces() |>
  melt() |>
  filter(name %in% c("Beta","sigma","gamma","loglik")) |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")

ggsave(paste0(k, "_if2_plot.png"), plot = plot_1, width = 10, height = 8, dpi = 300)

foreach(mf=mif_run,.combine=rbind,
        .options.future=list()
) %dofuture% {
  evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
  ll <- logmeanexp(evals,se=TRUE)
  mf |> coef() |> bind_rows() |>
    bind_cols(loglik=ll[1],loglik.se=ll[2])
} -> results

max_ll = -100000000
max_ll_idx = -1

for (i in 1:200){
  c_ll <- results[i, "loglik"]
  
  if(c_ll > max_ll) {
    max_ll = c_ll
    max_ll_idx = i
  }
}

ifs_beta = results[max_ll_idx, "Beta"]
ifs_sigma = results[max_ll_idx, "sigma"]
ifs_gamma = results[max_ll_idx, "gamma"]
ifs_mll = max_ll

#pairs(~loglik+Beta+sigma+gamma,data=results,pch=16)

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
      pfilter(params=p,Np=1e4)
  } %seed% TRUE |>
    concat()
} -> pfilter_run

max_ll = -100000000
max_ll_idx = -1

for (i in 1:200){
  c_ll <- pfilter_run[[i]]@loglik
  
  if(c_ll > max_ll) {
    max_ll = c_ll
    max_ll_idx = i
  }
}

pfilter_beta = pfilter_run[[i]]@params[1]
pfilter_sigma = pfilter_run[[i]]@params[2]
pfilter_gamma = pfilter_run[[i]]@params[3]
pfilter_mll = max_ll

## make a profile likelihood plot
plp <- bind_cols(
  mif_run |>
    coef() |>
    melt() |>
    pivot_wider(),
  logLik=logLik(mif_run)
) |>
  ggplot(aes(x=sigma,y=logLik))+
  geom_point()+
  theme_bw()

ggsave(paste0(k, "_if2_proflike.png"), plot = plp, width = 10, height = 8, dpi = 300)

ifs_sigma <- unname(ifs_sigma)
ifs_beta <- unname(ifs_beta)
ifs_gamma <- unname(ifs_gamma)
ifs_mll <- unname(ifs_mll)
pfilter_mll <- unname(pfilter_mll)

param_df <- rbind(param_df, data.frame(sigma_if2 = ifs_sigma[[1]], beta_if2 = ifs_beta[[1]], gamma_if2 = ifs_gamma[[1]], if2_ll = ifs_mll[[1]], pf_ll = pfilter_mll[[1]]))

}

