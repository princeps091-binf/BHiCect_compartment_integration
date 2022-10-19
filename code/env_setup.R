library(renv)

renv::init()

renv::install("tidyverse")
renv::install("vroom")
renv::install("Matrix")
renv::install("mgcv")
renv::install("furrr")
renv::install("data.tree")
renv::install("igraph")
renv::install("viridis")
renv::install("seriation")
renv::install("caret")

renv::install("bioc::GenomicRanges")
renv::snapshot()
