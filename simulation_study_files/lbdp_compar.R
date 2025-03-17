library(tidyverse)
library(pomp)
library(phylopomp)
library(cowplot)
library(broom)
library(scales)
stopifnot(getRversion() >= "4.2")
stopifnot(packageVersion("pomp")>="5.2")
stopifnot(packageVersion("phylopomp")>="0.9.2")
theme_set(theme_bw(base_family="serif"))
options(
  dplyr.summarise.inform=FALSE,
  pomp_archive_dir="results/lbdp_compar"
)
set.seed(505641575)

#'
#' The following codes simulate a realization of the genealogy process induced by a linear birth-death-sampling process with birth rate $\lambda$, death rate $\mu$ and sampling rate $\psi$.
#'
#'
## ----lbdp1--------------------------------------------------------------------
true.params <- data.frame(lambda=1.2,mu=0.8,psi=1,n0=5)
with(
  true.params,
  runLBDP(lambda=lambda,mu=mu,psi=psi,n0=n0,time=5)
) -> x
plot_grid(
  x |> plot(points=TRUE),
  x |> lineages() |> plot(),
  ncol=1,
  align="hv",
  rel_heights=c(2,1)
)

expand_grid(
  rep=1:10,
  lambda=1.2,
  mu=0.8,
  psi=seq(0.2,2.5,by=0.1),
  time=3,
  n0=5,
  Np=c(1000,2000,4000,8000,16000,32000,64000,128000)
) -> params

library(doParallel)
library(doRNG)
library(iterators)
registerDoParallel()
registerDoRNG(715604304)
foreach (
  p=iter(params,"row"),
  .inorder=FALSE,
  .combine=bind_rows
) %dopar% {
  with(
    p,
    x |>
      lbdp_exact(lambda=lambda,mu=mu,psi=psi,n0=n0)
  ) -> ll1
  with(
    p,
    x |>
      lbdp_pomp(lambda=lambda,mu=mu,psi=psi,n0=n0) |>
      pfilter(Np=Np) |>
      logLik()
  ) -> ll2
  bind_cols(p,exact=ll1,pf=ll2)
}-> params

params |>
  mutate(
    diff=(pf-exact)/exact
  ) |>
  group_by(Np) |>
  summarize(
    rmse=sqrt(mean(diff*diff)),
    sd=sd(diff),
    bias=abs(mean(exact-pf))
  ) |>
  ungroup() -> stats

stats |>
  pivot_longer(-Np) |>
  group_by(name) |>
  nest() |>
  reframe(tidy(lm(log(value)~log(Np),data=data[[1]]))) |>
  select(name,term,estimate,std.error) |>
  arrange(term,name)

pal <- c(viridis_pal(option="H",begin=0.1,end=0.8)(8),"#000000")
names(pal) <- c("1k","2k","4k","8k","16k","32k","64k","128k","exact")

plot_grid(
  params |>
    pivot_longer(c(exact,pf)) |>
    unite(name,name,Np) |>
    mutate(
      name=if_else(grepl("exact",name),"exact",name),
      name=gsub("pf_","",name),
      name=gsub("000","k",name),
      name=ordered(name,levels=names(pal))
    ) |>
    group_by(lambda,mu,psi,time,n0,name) |>
    reframe(
      type=c("logLik","logLik_se"),
      value=logmeanexp(value,se=TRUE)
    ) |>
    ungroup() |>
    pivot_wider(names_from=type) |>
    mutate(
      y=logLik,
      ymax=logLik+2*logLik_se,
      ymin=logLik-2*logLik_se
    ) |>
    filter(logLik>max(logLik)-16) |>
    ggplot(aes(x=psi,group=name,color=name,
      y=y,ymin=ymin,ymax=ymax))+
    geom_errorbar(
      position="dodge"
    )+
    geom_vline(xintercept=true.params$psi,linetype=2)+
    geom_hline(
      yintercept=max(params$exact)-
        c(0,0.5*qchisq(p=0.95,df=1)),
      linetype=2
    )+
    scale_color_manual(values=pal)+
    labs(
      color="effort",
      y="log likelihood",
      x=expression(psi)
    )+
    theme(
      legend.position="inside",
      legend.position.inside=c(0.5,0.4),
      legend.background=element_rect(fill="white")
    ),
  plot_grid(
    stats |>
      ggplot(aes(x=Np,y=rmse))+
      geom_smooth(formula=y~x,method="lm")+
      geom_point()+
      scale_x_log10(labels=\(x)aakmisc::scinot(x,simplify=TRUE))+
      scale_y_log10()+
      coord_fixed(ratio=1)+
      labs(x="num. particles",y="RMSE"),
    stats |>
      ggplot(aes(x=Np,y=sd))+
      geom_smooth(formula=y~x,method="lm")+
      geom_point()+
      scale_x_log10(labels=\(x)aakmisc::scinot(x,simplify=TRUE))+
      scale_y_log10()+
      coord_fixed(ratio=1)+
      labs(x="num. particles",y="SD"),
    stats |>
      ggplot(aes(x=Np,y=bias))+
      geom_smooth(formula=y~x,method="lm")+
      geom_point()+
      scale_x_log10(labels=\(x)aakmisc::scinot(x,simplify=TRUE))+
      scale_y_log10()+
      coord_fixed(ratio=1)+
      labs(x="num. particles",y="bias"),
    ncol=1,
    align="h",
    axis="btlr",
    rel_heights=c(1,1,1)
    ),
  nrow=1,
  rel_widths=c(2,1)
)

ggsave(filename="lbdp_plot.png",device=png,width=14,height=9)

