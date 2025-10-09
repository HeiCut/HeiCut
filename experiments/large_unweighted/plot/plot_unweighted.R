#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  if (!exists(".script_dir_initialized", envir=.GlobalEnv) || !.script_dir_initialized) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
    } else {
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- grep("^--file=", args, value = TRUE)
      wd <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
    }
    setwd(wd)
    assign(".script_dir_initialized", TRUE, envir=.GlobalEnv)
  }
}))
TEX = TRUE
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("common.R")
source("instances.R")
source("performance_profile.R")

ilp <- ilp
it0 <- it0
it1 <- it1
tight <- tight
trimmer <- trimmer

perf_comp_factors <- c(ilp_bname, it0_bname, it1_bname, tight_bname, trimmer_bname)

segments_default <- list(
  list(
    trans = "identity", 
    to = 1.1, 
    width = 2,
    breaks = seq(1.0, 1.1, by = 0.01),
    labels = c("1.0", "", "", "", "", "1.05", "", "", "", "", "1.1")),
  list(
    trans = "identity", 
    to = 2, 
    width = 2,
    breaks = c(1.25, 1.5, 1.75, 2),
    labels = c("", "1.5", "", "2")),
  list(
    trans = "log10", 
    to = 100,
    width = 2,
    breaks = c(10, 100),
    labels = c("10^1", "10^2"))
)

segments_cut <- list(
  list(
    trans = "identity", 
    to = 1.1, 
    width = 2,
    breaks = seq(1.0, 1.1, by = 0.01),
    labels = c("1.0", "", "", "", "", "1.05", "", "", "", "", "1.1")),
  list(
    trans = "identity", 
    to = 2, 
    width = 2,
    breaks = c(1.25, 1.5, 1.75, 2),
    labels = c("", "1.5", "", "2")),
  list(
    trans = "log10", 
    to = 100,
    width = 2,
    breaks = c(10, 100),
    labels = c("10^1", "10^2"))
)

mem_segments_adil <- list(
  list(
    trans = "identity",
    to = 7,
    width = 6,
    breaks = seq(1.0, 7.0, by = 1.0),
    labels = c("1", "2", "3", "4", "5", "6", "7")
  )
)
time_segments_adil <- list(
  list(
    trans = "identity",
    to = 90,
    width = 6,
    breaks = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90),
    labels = c("", "10", "20", "30", "40", "50", "60", "70", "80", "90")
  )
)
cut_segments_adil <- list(
  list(
    trans = "log2",
    to = 1024,
    width = 6,
    breaks = c(1, 4, 16, 64, 256, 1024),
    labels = c("1", "4", "16", "64", "256", "1024")
  )
)

segments_my_mem <- list(
  list(
    trans = "identity", 
    to = 2.0, 
    width = 2,
    breaks = seq(1.0, 2.0, by = 0.25),
    labels = c("1", "", "1.5", "", "2")
  ),
  list(
    trans = "identity", 
    to = 10, 
    width = 2,
    breaks = seq(2, 10, by = 2),
    labels = c("", "", "6", "", "")
  ),
  list(
    trans = "log10", 
    to = 100000,
    width = 2,
    breaks = c(10, 100, 1000, 10000,100000),
    labels = c("10^1", "", "10^3", "", "10^5")
  )
)

segments_my_time <- list(
  list(
    trans = "identity", 
    to = 2.0, 
    width = 2,
    breaks = seq(1.0, 2.0, by = 0.25),
    labels = c("1", "", "1.5", "", "2")
  ),
  list(
    trans = "identity", 
    to = 10, 
    width = 2,
    breaks = seq(2, 10, by = 2),
    labels = c("", "", "6", "", "")
  ),
  list(
    trans = "log10", 
    to = 1000000000,
    width = 2,
    breaks = c(10, 1000, 100000, 10000000, 1000000000),
    labels = c("10^1", "", "10^5", "", "10^9")
  )
)

segments_my_c <- list(
  list(
    trans = "identity", 
    to = 2.0, 
    width = 2,
    breaks = seq(1.0, 2.0, by = 0.25),
    labels = c("1", "", "1.5", "", "2")
  ),
  list(
    trans = "identity", 
    to = 10, 
    width = 2,
    breaks = seq(2, 10, by = 2),
    labels = c("", "", "6", "", "")
  ),
  list(
    trans = "log10", 
    to = 100000,
    width = 2,
    breaks = c(10, 100, 1000, 10000,100000),
    labels = c("10^1", "", "10^3", "", "10^5")
  )
)

mem_segments <- segments_my_mem #segments_default # mem_segments_adil
time_segments <- segments_my_time #segments_default # time_segments_adil
cut_segments <- segments_my_c #segments_cut # cut_segments_adil

perf_comp_all_mem <- create_performance_profile(
  ilp, 
  it0, 
  it1, 
  tight,
  trimmer,
  column.objective = "Memory",
  tex = TEX,
  tiny = TRUE,
  colors = colors_perf,
  factor.levels = perf_comp_factors,
  segments = mem_segments
) +
  theme_bw(base_size = 10) +
  create_theme() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    legend.spacing.x = unit(5, "mm"),
    legend.spacing.y = unit(2, "mm"),
    legend.key.spacing.x = unit(5, "mm"),
    legend.key.spacing.y = unit(2, "mm"),
    legend.box.spacing = unit(2, "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    #legend.text = element_text(margin = margin(l = 5, r = 5)),
    plot.margin = margin(t = 0, r = 2.5, b = 0, l = 0, unit = "mm")
  ) + 
  guides(color = guide_legend(nrow = 1)) +
  ggtitle("a) Memory")

perf_comp_all_time <- create_performance_profile(
  ilp, 
  it0, 
  it1, 
  tight,
  trimmer,
  column.objective = "Time",
  tex = TEX,
  tiny = TRUE,
  colors = colors_perf,
  factor.levels = perf_comp_factors,
  segments = time_segments
) +
  theme_bw(base_size = 10) +
  create_theme() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    legend.spacing.x = unit(-2, "mm"),
    legend.spacing.y = unit(-2, "mm"),
    legend.key.spacing.x = unit(-3, "mm"),
    legend.key.spacing.y = unit(-2, "mm"),
    legend.box.spacing = unit(-2, "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    plot.margin = margin(t = 0, r = 2.5 / 2, b = 0, l = 2.5 / 2, unit = "mm")
  ) + 
  guides(color = guide_legend(nrow = 1)) +
  ggtitle("b) Running Time")

perf_comp_all_cut <- create_performance_profile(
  ilp, 
  it0, 
  it1, 
  tight,
  trimmer,
  column.objective = "Cut",
  tex = TEX,
  tiny = TRUE,
  colors = colors_perf,
  factor.levels = perf_comp_factors,
  segments = cut_segments
) +
  theme_bw(base_size = 10) +
  create_theme() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    legend.spacing.x = unit(-2, "mm"),
    legend.spacing.y = unit(-2, "mm"),
    legend.key.spacing.x = unit(-3, "mm"),
    legend.key.spacing.y = unit(-2, "mm"),
    legend.box.spacing = unit(-2, "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 2.5, unit = "mm")
  ) + 
  guides(color = guide_legend(nrow = 1)) + 
  ggtitle("c) Min Cut")

perf_comp_all_legend <- ggpubr::get_legend(perf_comp_all_mem, position = "bottom")
perf_comp_all_mem <- perf_comp_all_mem + theme(legend.position = "none")

open_dev("perf_large_unw", width = 7 / 1.05, tex = TRUE)
egg::ggarrange(
  perf_comp_all_mem, perf_comp_all_time, perf_comp_all_cut,
  #heights = 1,
  nrow = 1
  #padding = unit(5.25, "line")
)
dev_off()
open_dev("perf_large_unw_legend", width = 7 / 1.05, tex = TRUE)
egg::ggarrange(
  ggpubr::as_ggplot(perf_comp_all_legend), 
  heights = 1.5,
  nrow = 1,
  padding = unit(0.25, "line")
)
dev_off()

