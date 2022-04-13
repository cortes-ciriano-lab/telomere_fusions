install.packages(pkgs = c("data.table","plyr","dplyr","stringr","tibble","argparse","ggplot2","stringdist"), repos = "http://cran.us.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install(c("Rsamtools","GenomicAlignments"))

