library(tidyverse)
library(phylopomp)
library(pomp)
library(tictoc)
library(ape)
theme_set(theme_bw())

setwd("~/projects/yang-phylopomp")


mif2 <- readRDS("mif2_best_reduced_pop_results.rds")

mif2|>
  traces() |>
  melt() |>
  filter(name %in% c("Beta11","Beta21","Beta22","psi1", "psi2",
                     "omega1", "omega2",
                     "loglik")) |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")
