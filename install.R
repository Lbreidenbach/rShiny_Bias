update.packages(ask = FALSE, checkBuilt = TRUE)
install.packages(c(
    "BiocManager",
    "devtools",
    "DiagrammeR",
    "dplyr",
    "GGally",
    "ggplot2",
    "ggridges",
    "grid",
    "gridExtra",
    "gtable",
    "MatchIt",
    "plyr",
    "purrr",
    "reshape2",
    "rjags",
    "shiny"
    ))

BiocManager::install("graph")

devtools::install_github("nutterb/HydeNet@fcbb7d81f2359b98494f0712a5db15291193ae5f")
# remotes::install_url("https://cran.r-project.org/src/contrib/Archive/HydeNet/HydeNet_0.10.11.tar.gz")
