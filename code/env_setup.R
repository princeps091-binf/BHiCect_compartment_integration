library(renv)

renv::init()

renv::install("tidyverse")
renv::install("vroom")
renv::install("Matrix")
renv::install("mgcv")

renv::install("bioc::GenomicRanges")
renv::snapshot()
