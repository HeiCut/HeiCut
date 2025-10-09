#!/usr/bin/env Rscript
.bootstrap_deps <- function() {
  # Required packages for plotting pipeline
  required_pkgs <- c(
    "ggplot2","plyr","dplyr","RColorBrewer","tikzDevice",
    "gridExtra","egg","ggpubr","stringr","stringi","ggrepel"
  )

  install_if_missing <- function(pkgs) {
    miss <- pkgs[!pkgs %in% rownames(installed.packages())]
    if (length(miss)) {
      cat("[deps] Installing missing R packages:", paste(miss, collapse=", "), "...\n")
      install.packages(miss, repos = "https://cloud.r-project.org")
    }
  }

  try_load_pkg <- function(pkgname) {
    tryCatch({
      suppressPackageStartupMessages(
        library(pkgname, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      )
      TRUE
    }, error = function(e) {
      assign(".last_pkg_error", e$message, .GlobalEnv)
      FALSE
    })
  }

  # 1) Ensure required packages are installed
  install_if_missing(required_pkgs)

  # 2) Special handling for stringi ICU mismatch
  if (!try_load_pkg("stringi")) {
    msg <- get(".last_pkg_error", envir = .GlobalEnv, inherits = FALSE)
    if (grepl("ICU|libicu|libicui18n", msg, ignore.case = TRUE)) {
      cat("[deps] Detected ICU mismatch for 'stringi'. Rebuilding from source (bundled ICU)…\n")
      install.packages("stringi", type = "source", repos = "https://cloud.r-project.org")
      if (!try_load_pkg("stringi")) {
        stop("[deps] Failed to (re)install 'stringi'.\n",
             "If possible, install ICU runtime (e.g. apt: libicu-dev) or ensure you have build tools,\n",
             "then rerun. Last error: ", msg, call. = FALSE)
      }
    } else {
      stop("[deps] Could not load 'stringi': ", msg, call. = FALSE)
    }
  }

  # 3) Ensure stringr loads (depends on stringi)
  if (!try_load_pkg("stringr")) {
    installIf <- tryCatch(install.packages("stringr", repos = "https://cloud.r-project.org"), error=function(e) e)
    stopifnot(try_load_pkg("stringr"))
  }
}

.bootstrap_deps()

options(show.error.locations = TRUE)
options(error = traceback)

library(ggplot2, warn.conflicts = FALSE)
library(plyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(RColorBrewer, warn.conflicts = FALSE)
library(tikzDevice, warn.conflicts = FALSE)
library(gridExtra, warn.conflicts = FALSE)
library(egg, warn.conflicts = FALSE)
library(ggpubr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(ggrepel, warn.conflicts = FALSE)

DATA_INPUT_DIR <- "../experiments"
CACHE_DIR <- "../cache/experiments"

TEX_INPUT <- "preamble.tex"

TEX_OUTPUT <- "../plot/"
DATA_OUTPUT <- "../data/"
PDF_OUTPUT <- "../plot/"

PDF_LABEL_TIMEOUT <- "T"
PDF_LABEL_INFEASIBLE <- "I"
PDF_LABEL_FAILED <- "F"

TEX_LABEL_TIMEOUT <- "T"
TEX_LABEL_INFEASIBLE <- "F"
TEX_LABEL_FAILED <- "I"

TIMELIMIT <- 8880000 
EPSILON <- 0.03

aggregator <- function(df)
    data.frame(
        AvgCut = ifelse(all(is.na(df$Cut)), NA, mean(df$Cut, na.rm = TRUE)),
        MinBalance = ifelse(all(is.na(df$Balance)), NA, min(df$Balance, na.rm = TRUE)),
        AvgTime = ifelse(all(is.na(df$Time)), NA, mean(df$Time, na.rm = TRUE)),
        Timeout = all(df$Timeout),
        Failed = all(df$Failed)
    )

normalize_pareto <- function(df) df %>% 
    dplyr::mutate(Graph = sub("\\.metis|\\.mtx|\\.mtx.hgr|\\.hgr|\\.graph|\\.scotch", "", Graph)) %>%
    dplyr::mutate(Cut = ifelse(Balance > EPSILON + .Machine$double.eps, NA, Cut)) %>%
    dplyr::mutate(Timeout = (as.numeric(Time) >= TIMELIMIT)) %>%
    dplyr::mutate(Cut = ifelse(Timeout, NA, Cut)) %>%
    dplyr::mutate(Balance = ifelse(Timeout, 0.0, Balance)) %>%
    dplyr::mutate(Time = ifelse(Timeout, TIMELIMIT, Time)) %>%
    dplyr::mutate(Failed = (Failed == "1")) %>%
    dplyr::mutate(Cut = ifelse(Failed, NA, Cut)) %>%
    dplyr::mutate(Time = ifelse(Failed, TIMELIMIT, Time)) %>%
    ddply(c("Algorithm", "Graph", "K", "Parallelism", "Epsilon", "Threads"), aggregator) %>%
    dplyr::mutate(AvgCut = ifelse(is.na(AvgCut), Inf, AvgCut)) %>% # *all* repetitions failed / timeouted / are infeasible: set cut to Inf
    dplyr::mutate(Infeasible = !Failed & (MinBalance != 1.0 | AvgTime < TIMELIMIT) & MinBalance > EPSILON + .Machine$double.eps) %>% 
    dplyr::mutate(Invalid = Failed | Timeout | Infeasible)

normalize_perf <- function(df) df %>% 
    dplyr::mutate(Graph = sub("\\.metis|\\.mtx|\\.mtx.hgr|\\.hgr|\\.graph|\\.scotch", "", Graph)) %>%
    dplyr::mutate(Failed = ifelse(Failed == "1", TRUE, FALSE)) %>%
    dplyr::mutate(Timeout = ifelse(Timeout == "1", TRUE, FALSE)) %>%
    dplyr::mutate(Balance = ifelse(Failed | Timeout, NA, Balance)) %>%  
    dplyr::mutate(Time = ifelse(Timeout, TIMELIMIT, ifelse(Failed, NA, Time))) %>% 
    dplyr::mutate(Cut = ifelse(Failed | Timeout | Balance > EPSILON + .Machine$double.eps, NA, Cut)) %>%
    plyr::ddply(c("Algorithm", "Graph", "K", "Parallelism", "Epsilon", "Threads"), aggregator) %>%
    dplyr::mutate(AvgCut = ifelse(is.na(AvgCut), Inf, AvgCut)) %>%       
    dplyr::mutate(Infeasible = (AvgCut == Inf)) %>% 
    dplyr::mutate(Feasible = !Failed & !Timeout & !Infeasible) %>%
    dplyr::mutate(Invalid = Failed | Timeout | Infeasible)

load_dataset <- function(filename, print_name, parallelism, aggregator = normalize_perf) {
    in_filename <- paste0(DATA_INPUT_DIR, "/", filename)
    cache_filename <- paste0(CACHE_DIR, "/", gsub("/", "_", filename))

    cat(paste0("Loading algorithm ", print_name, ":\n")) 
    cat(paste0("  input: ", in_filename, "\n"))
    cat(paste0("  cache: ", cache_filename, "\n"))

    if (file.exists(cache_filename)) {
        return(read.csv(cache_filename))
    }

    df <- read.csv(in_filename)

    if (!("Threads" %in% colnames(df))) {
        df$Threads <- 1
    }
    if (!("Failed" %in% colnames(df))) {
        df$Failed <- 0
    }
    if (!("Timeout" %in% colnames(df))) {
        df$Timeout <- 0
    }
    if (!("Epsilon" %in% colnames(df))) {
        df$Epsilon <- 0.03
    }

    df$Algorithm <- print_name
    df$Parallelism <- parallelism
    cat(paste0("  num threads: ", unique(df$Threads), "\n"))
    cat(paste0("  rows before aggregation: ", nrow(df), "\n"))
    ndf <- aggregator(df)
    cat(paste0("  rows after aggregation: ", nrow(ndf), "\n"))

    write.csv(ndf, cache_filename, row.names = FALSE)

    return(ndf)
}

DEFAULT_ASPECT_RATIO <- 2 / (1 + sqrt(5))

create_theme <- function(aspect_ratio = DEFAULT_ASPECT_RATIO)
    theme(aspect.ratio = aspect_ratio,
          legend.background = element_blank(),
          legend.title = element_blank(),
          #legend.margin = margin(-5, 0, 0, 0),
          #legend.spacing.x = unit(0.01, "cm"),
          #legend.spacing.y = unit(0.01, "cm"),
          legend.box.spacing = unit(0.1, "cm"),
          legend.title.align = 0.5,
          legend.text = element_text(size = 8, color = "black"), 
          plot.title = element_text(size = 10, hjust = 0.5, color = "black"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_line(linetype = "11", linewidth = 0.3, color = "grey"),
          panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 0.2, color = "black"),
          axis.title.y = element_text(size = 8, vjust = 1.5, color = "black"),
          axis.title.x = element_text(size = 8, vjust = 1.5, color = "black"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"))

if (!exists("tikzDeviceLoaded")) {
  options(tikzLatexPackages = c(getOption('tikzLatexPackages'), 
                                "\\usepackage[utf8]{inputenc}",  # Add this line
                                paste0("\\input{", getwd(), "/", TEX_INPUT, "}")))
  tikzDeviceLoaded <- T
  options(tikzMetricsDictionary = ".cache.db")
}


current_device_file <- ""
current_device_file_is_tikz <- FALSE

open_pdf <- function(file, width = 7) {
    current_device_file <<- paste0(PDF_OUTPUT, "/", file, ".pdf") 
    current_device_file_is_tikz <<- FALSE
    pdf(current_device_file, width = width)
}

open_tikz <- function(file, width = 7) {
    current_device_file <<- paste0(TEX_OUTPUT, "/", file, ".tex")
    current_device_file_is_tikz <<- TRUE
    if (width == 7) {
        tikz(current_device_file, pointsize = 12)
    } else {
        tikz(current_device_file, width = width, pointsize = 12)
    }
}

open_dev <- function(file, width = 7, tex = FALSE) {
    if (tex) {
        open_tikz(file, width)
    } else {
        open_pdf(file, width)
    }
}

open_sink <- function(file) {
    sink(paste0(DATA_OUTPUT, "/", file, "tex"))
}

dev_off <- function() {
    dev.off()
    if (current_device_file_is_tikz) {
        lines <- readLines(con = current_device_file)
        lines <- lines[-which(grepl("\\path\\[clip\\]*", lines, perl = F))]
        lines <- lines[-which(grepl("\\path\\[use as bounding box*", lines, perl = F))]
        writeLines(lines, con = current_device_file)
    }
}

sink_off <- function() {
    sink()
}

make_var <- function (name, value)  
    cat(paste0("\\newcommand{\\", name, "}{", value, "}\n"))
make_varf2 <- function (name, value) 
    cat(paste0("\\newcommand{\\", name, "}{", sprintf("%.2f", value), "}\n"))
make_varf1 <- function (name, value) 
    cat(paste0("\\newcommand{\\", name, "}{", sprintf("%.1f", value), "}\n"))
make_varf0 <- function (name, value) 
    cat(paste0("\\newcommand{\\", name, "}{", sprintf("%.0f", value), "}\n"))

gm_mean <- function(x, na.rm = TRUE, zero.propagate = FALSE){
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }

  if (zero.propagate) {
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    return(exp(mean(log(x[x != Inf]), na.rm = na.rm)))
  } else {
    return(exp(sum(log(x[x > 0 & x != Inf]), na.rm=na.rm) / length(x)))
  }
}

hm_mean <- function(x) {
  return(length(x) / sum( 1.0 / x[x > 0] ))
}

axis_text_size <- function(latex_export = F, small_size = F) {
  if (small_size) {
    return(6)
  } else if (latex_export) {
    return(8)
  } else {
    return(10)
  }
}

legend_text_size <- function(latex_export = F, small_size = F) {
  if (small_size) {
    return(6)
  } else if (latex_export) {
    return(8)
  } else {
    return(10)
  }
}

axis_title_size <- function(latex_export = F, small_size = F) {
  if (small_size) {
    return(6)
  } else if (latex_export) {
    return(8)
  } else {
    return(10)
  }
}

plot_text_size <- function(latex_export = F) {
  if (latex_export) {
    return(3)
  } else {
    return(2)
  }
}

plot_title_size <- function(latex_export = F) {
  if (latex_export) {
    return(10)
  } else {
    return(12)
  }
}

plot_point_size <- function(latex_export = F) {
  if (latex_export) {
    return(0.5)
  } else {
    return(1)
  }
}

plot_line_size <- function(latex_export = F) {
  if (latex_export) {
    return(0.5)
  } else {
    return(0.75)
  }
}

to_latex_math_mode <- function(x, latex_export = F) {
  if (latex_export) {
    return(paste("$",x,"$",sep=""))
  } else {
    return(x)
  }
}

pow_text <- function(base, exp, latex_export = F) {
  if ( latex_export ) {
    x <- paste(base, "^{", exp, "}", sep="")
  } else {
    x <- paste(base, "^", exp, sep="")
  }
  return(to_latex_math_mode(x, latex_export))
}

add_infeasible_break <- function(breaks, infeasible_value, show_infeasible_tick = F, latex_export = F) {
  if ( show_infeasible_tick ) {
    breaks <- c(breaks, infeasible_value)
  }
  return(breaks)
}

add_infeasible_label <- function(labels, infeasible_value, show_infeasible_tick = F, latex_export = F) {
  if ( show_infeasible_tick ) {
    if ( latex_export ) {
      labels <- c(labels, "\\ding{55}")
    } else {
      labels <- c(labels, "infeasible")
    }
  }
  return(labels)
}

add_timeout_break <- function(breaks, timeout_value, show_timeout_tick = F, latex_export = F) {
  if ( show_timeout_tick ) {
    breaks <- c(breaks, timeout_value)
  }
  return(breaks)
}

add_timeout_label <- function(labels, timeout_value, show_timeout_tick = F, latex_export = F) {
  if ( show_timeout_tick ) {
    if ( latex_export ) {
      labels <- c(labels, "\\ClockLogo")
    } else {
      labels <- c(labels, "timeout")
    }
  }
  return(labels)
}
