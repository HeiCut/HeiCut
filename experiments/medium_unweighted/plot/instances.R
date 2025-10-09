#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  if (!exists(".script_dir_initialized", envir=.GlobalEnv) || !.script_dir_initialized) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
    } else {
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- grep("^--file=", args, value = TRUE)
      if (length(file_arg) == 0) {
        wd <- getwd()  # fallback
      } else {
        wd <- normalizePath(dirname(sub("^--file=", "", file_arg[1])))
      }
    }
    setwd(wd)
    assign(".script_dir_initialized", TRUE, envir=.GlobalEnv)
  }
}))
source("common.R")

tex_mode <- exists("TEX") && TEX == TRUE

excluded_graphs <- c(
    "hcircuit",
    "coAuthorsDBLP",
    "coAuthorsCiteseer",
    "citationCiteseer",
    "wiki-Talk"
)

test_graphs <- c(
    "arabic-2005",
    "RHG-100000000-nodes-1000000000-edges",
    "RHG-100000000-nodes-2000000000-edges",
    "com-friendster",
    "it-2004-sorted",
    "nlpkkt240",
    "orkut",
    "rgg_n26",
    "sk-2005-sorted",
    "twitter-2010",
    "uk-2007-05",
    "webbase-2001"
)

`%notin%` <- Negate(`%in%`)

get_algo_name <- function(df) df$Algorithm[1]

print_algo <- function(df) paste0(
    " - ", df$Algorithm[1], ": ", nrow(df),
    "\n"
)

load_algorithm <- function(filename, print_name) {
    in_filename <- paste0("../generated/all_results/", filename, ".csv")

    cat(paste0("Loading algorithm ", print_name, ":\n")) 
    cat(paste0("  input: ", in_filename, "\n"))

    df <- read.csv(in_filename)
    df$Algorithm <- print_name
    df$RawAlgorithm <- df$algorithm_name
    df$Graph <- df$graph_name
    df$K <- df$k_value
    df$Time <- df$runtime
    df$Memory <- df$memory
    df$Cut <- df$mincut
    df$Failed <- df$Failed
    df$Timeout <- FALSE
    df$Feasible <- TRUE
    df$Infeasible <- FALSE

    #df <- df %>% dplyr::filter(Graph %notin% excluded_graphs)

    cat(paste0("  nrows: ", nrow(df), "\n"))

    return(df)
}

ilp_bname <- "Relaxed-BIP"
it0_bname <- "HeiCut"
it1_bname <- "HeiCut-LP"
tight_bname <- "Tight"
trimmer_bname <- "Trimmer"


ilp <- load_algorithm("ilp", ilp_bname)
it0 <- load_algorithm("kernelizer_IT0", it0_bname)
it1 <- load_algorithm("kernelizer_IT1", it1_bname)
tight <- load_algorithm("submodular_tight_single", tight_bname)
trimmer <- load_algorithm("trimmer", trimmer_bname)


unique(it1$Graph)

palette <- brewer.pal(n = 9, name = "Set1")
palette2 <- brewer.pal(n = 7, name = "Set2")

colors_perf <- c()
colors_perf[ilp_bname] = palette[[1]]
colors_perf[it0_bname] = palette[[2]]
colors_perf[it1_bname] = palette[[3]]
colors_perf[tight_bname] = palette[[4]]
colors_perf[trimmer_bname] = palette[[5]]
#colors_perf[rle_b1_k15_bname] = palette[[7]]
#colors_perf[rle_b1_k20_bname] = palette[[8]]
#colors_perf[rle_b1_k25_bname] = palette2[[6]]
#colors_perf[expq_bname] = palette2[[7]]
#colors_perf[vec_bname] = "black" # palette[[9]]
#colors_perf[c_bname] = "magenta" # palette[[9]]
#colors_perf[cl_bname] = "darkgreen" # palette[[9]]

