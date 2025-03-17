library(tidyverse)
library(pomp)
library(phylopomp)

time = 5

runTwoSpecies(
  time,
  t0 = 0,
  Beta11 = 4,
  Beta12 = 0,
  Beta21 = 0,
  Beta22 = 4,
  gamma1 = 1,
  gamma2 = 1,
  psi1 = 1,
  psi2 = 0,
  c1 = 1,
  c2 = 1,
  omega1 = 0,
  omega2 = 0,
  b1 = 0,
  b2 = 0,
  d1 = 0,
  d2 = 0,
  iota1 = 0,
  iota2 = 0,
  S1_0 = 100,
  S2_0 = 100,
  I1_0 = 0,
  I2_0 = 10,
  R1_0 = 0,
  R2_0 = 0
) -> s

s |>
  twospecies_pomp(
    Beta11=4,Beta12=1,
    Beta21=1,Beta22=4,
    gamma1=1,gamma2=1,
    c1=0.9,c2=0.9,
    psi1=1,psi2=1,
    S1_0=100,I1_0=2,R1_0=0,
    S2_0=100,I2_0=2,R2_0=0,
    b1=0.1,b2=0.1,
    d1=0.1,d2=0.1,
    omega1=0.5,omega2=0.5
  ) |>
  pfilter(Np=100) |>
  freeze(seed=468314464) |>
  logLik()

s |>
  twospecies_pomp(
    Beta11=4,Beta12=1,
    Beta21=1,Beta22=4,
    gamma1=1,gamma2=1,
    c1=0.9,c2=0.9,
    psi1=1,psi2=1,
    S1_0=100,I1_0=5,R1_0=0,
    S2_0=100,I2_0=5,R2_0=0,
    b1=0.1,b2=0.1,
    d1=0.1,d2=0.1,
    omega1=0.5,omega2=0.5
  ) -> p

p |>
  pfilter(Np=1000) |>
  replicate(n=5) |>
  concat() |>
  freeze(seed=468314464) -> pf

pf |>
  logLik() |>
  logmeanexp(se=TRUE,ess=TRUE)

try(
  s |>
    twospecies_pomp(
      S1_0=20,I1_0=-5,R1_0=0,
      S2_0=20,I2_0=5,R2_0=0
    )
)

stopifnot(
  s |>
    twospecies_pomp(
      Beta11=4,Beta12=1,
      Beta21=1,Beta22=4,
      gamma1=1,gamma2=1,
      c1=0.9,c2=0.9,
      psi1=1,psi2=1,
      S1_0=20,I1_0=5,R1_0=0,
      S2_0=20,I2_0=5,R2_0=0,
      b1=0.1,b2=0.1,
      d1=0.1,d2=0.1,
      omega1=0.5,omega2=0.5
    ) |>
    pfilter(Np=100) |>
    freeze(seed=407324464) |>
    logLik()==-Inf
)