#!/usr/bin/env Rscript
create_performance_profile <- function(...,
                                       column.objective = "AvgCut", 
                                       column.algorithm = "Algorithm", 
                                       column.timeout = "Timeout",
                                       column.infeasible = "Infeasible",
                                       column.failed = "Failed",
                                       primary_key = c("Graph", "K"),
                                       segments = list(
                                         list(trans = "identity", 
                                              to = 1.1, 
                                              width = 3,
                                              breaks = seq(1.0, 1.1, by = 0.01),
                                              labels = c("1.0", "", "", "", "", "1.05", "", "", "", "", "1.1")),
                                         list(trans = "identity", 
                                              to = 2, 
                                              width = 2,
                                              breaks = c(1.25, 1.5, 1.75, 2),
                                              labels = c("", "1.5", "", "2")),
                                         list(trans = "log10", 
                                              to = 100,
                                              width = 1,
                                              breaks = c(10, 100),
                                              labels = c("10^1", "10^2"))
                                       ),
                                       segment.errors.width = 1,
                                       tick.timeout = "auto", 
                                       tick.infeasible = "auto",
                                       tick.failed = "auto",
                                       label.pdf.timeout = PDF_LABEL_TIMEOUT,
                                       label.pdf.infeasible = TEX_LABEL_INFEASIBLE,
                                       label.pdf.failed = PDF_LABEL_FAILED,
                                       label.tex.timeout = TEX_LABEL_TIMEOUT,
                                       label.tex.infeasible = TEX_LABEL_INFEASIBLE,
                                       label.tex.failed = TEX_LABEL_FAILED,
                                       colors = c(),
                                       tex = FALSE,
                                       tiny = FALSE,
                                       factor.levels = c()) {
  
  all_datasets <- list(...)
  stopifnot(length(all_datasets) > 0)
  
  # Sort by primary key & mark failed as Inf for objective
  for (i in seq_along(all_datasets)) {
    all_datasets[[i]] <- all_datasets[[i]] %>% dplyr::arrange_at(primary_key)
    all_datasets[[i]][[column.objective]] <-
      ifelse(all_datasets[[i]][[column.failed]], Inf, all_datasets[[i]][[column.objective]])
  }
  
  # Consistency checks
  first_dataset <- all_datasets[[1]]
  for (dataset in all_datasets) {
    stopifnot(column.objective %in% colnames(dataset))
    stopifnot(column.algorithm %in% colnames(dataset))
    stopifnot(column.timeout %in% colnames(dataset))
    stopifnot(column.infeasible %in% colnames(dataset))
    stopifnot(!(-Inf %in% dataset[[column.objective]]))
    stopifnot(nrow(dataset) == nrow(first_dataset))
    stopifnot(dataset[, primary_key] == first_dataset[, primary_key])
  }
  
  # Ratios vs best
  best <- do.call(pmin, lapply(all_datasets, \(df) df[[column.objective]]))
  all_ratios <- lapply(all_datasets, \(df) df[[column.objective]] / best)  # list per algorithm
  
  # Build aggregated profile per algorithm
  PSEUDO_RATIO_TIMEOUT    <- 30000000000
  PSEUDO_RATIO_INFEASIBLE <- 20000000000
  PSEUDO_RATIO_FAILED     <- 10000000000
  
  num_rows <- nrow(first_dataset)
  pp_data <- data.frame()
  for (i in seq_along(all_datasets)) {
    dataset <- all_datasets[[i]]
    ratios <- all_ratios[[i]]
    tmp <- data.frame(Ratio = ratios) %>%
      dplyr::mutate(Ratio = ifelse(dataset[[column.timeout]],    PSEUDO_RATIO_TIMEOUT,    Ratio)) %>%
      dplyr::mutate(Ratio = ifelse(dataset[[column.infeasible]], PSEUDO_RATIO_INFEASIBLE, Ratio)) %>%
      dplyr::mutate(Ratio = ifelse(dataset[[column.failed]] & !dataset[[column.timeout]],
                                   PSEUDO_RATIO_FAILED, Ratio)) %>%
      dplyr::group_by(Ratio) %>%
      dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(Ratio) %>%
      dplyr::mutate(Fraction = cumsum(N) / num_rows,
                    Algorithm = dataset[[column.algorithm]][1],
                    Transformed = 0)
    pp_data <- rbind(tmp, pp_data)
  }
  
  # Piecewise transformation across numeric segments
  offset <- 0
  from <- 1
  for (segment in segments) {
    min_value <- do.call(segment$trans, list(from))
    max_value <- do.call(segment$trans, list(segment$to))
    span <- max_value - min_value
    
    pp_data <- pp_data %>%
      dplyr::mutate(Transformed = ifelse(Transformed == 0 & Ratio >= from & Ratio < segment$to, 1, Transformed)) %>%
      dplyr::mutate(Ratio = ifelse(Transformed == 1,
                                   offset + segment$width * (do.call(segment$trans, list(Ratio)) - min_value) / span,
                                   Ratio)) %>%
      dplyr::mutate(Transformed = ifelse(Transformed == 1, 2, Transformed))
    
    offset <- offset + segment$width
    from <- segment$to
  }
  
  # Error segment mapping
  map_errors <- function(vals)
    sapply(vals, \(val)
           if (val == PSEUDO_RATIO_TIMEOUT) 2
           else if (val == PSEUDO_RATIO_INFEASIBLE) 1
           else if (val == PSEUDO_RATIO_FAILED) 1
           else if (val == Inf) 2
           else 0)
  
  pp_data <- pp_data %>%
    dplyr::mutate(Transformed = ifelse(Transformed == 0 & Ratio >= from, 1, Transformed)) %>%
    dplyr::mutate(Ratio = ifelse(Transformed == 1,
                                 offset + segment.errors.width * map_errors(Ratio) / 2,
                                 Ratio)) %>%
    dplyr::mutate(Transformed = ifelse(Transformed == 1, 2, Transformed))
  
  stopifnot(pp_data %>% dplyr::filter(Transformed != 2) %>% nrow() == 0)
  
  # Ensure vertical jump at x=1.0 where needed
  transformed_one <- NA_real_
  offset_tmp <- 0
  from_tmp <- 1
  for (segment in segments) {
    if (is.na(transformed_one) && 1 >= from_tmp && 1 < segment$to) {
      min_val <- do.call(segment$trans, list(from_tmp))
      span <- do.call(segment$trans, list(segment$to)) - min_val
      transformed_one <- offset_tmp + segment$width * (do.call(segment$trans, list(1)) - min_val) / span
    }
    offset_tmp <- offset_tmp + segment$width
    from_tmp <- segment$to
  }
  if (is.na(transformed_one)) transformed_one <- 0
  
  algo_names <- sapply(all_datasets, \(d) unique(d[[column.algorithm]]))
  need_jump <- logical(length(algo_names)); names(need_jump) <- algo_names
  for (i in seq_along(all_datasets)) {
    raw_ratios <- all_ratios[[i]]
    if (any(raw_ratios == 1, na.rm = TRUE)) {
      need_jump[ all_datasets[[i]][[column.algorithm]][1] ] <- TRUE
    }
  }
  if (any(need_jump)) {
    add_rows <- lapply(names(need_jump)[need_jump], function(algo) {
      alg_rows <- pp_data[pp_data$Algorithm == algo, ]
      tol <- 1e-12
      at_one <- alg_rows[abs(alg_rows$Ratio - transformed_one) < tol, ]
      if (nrow(at_one) > 0 && all(at_one$Fraction > 0) ) {
        data.frame(Algorithm = algo,
                   Ratio = transformed_one,
                   Fraction = 0,
                   Transformed = 2,
                   stringsAsFactors = FALSE)
      } else NULL
    })
    add_rows <- do.call(rbind, add_rows)
    if (!is.null(add_rows) && nrow(add_rows) > 0) {
      pp_data <- dplyr::bind_rows(pp_data, add_rows) %>%
        dplyr::arrange(Algorithm, Ratio, Fraction)
    }
  }
  
  # ### ADD: all-failed vertical line fix ###
  # Identify algorithms that failed/timeout/infeasible on *every* instance.
  all_failed_algos <- character(0)
  for (i in seq_along(all_datasets)) {
    ds <- all_datasets[[i]]
    invalid_flags <- (ds[[column.failed]] | ds[[column.timeout]] | ds[[column.infeasible]])
    if (all(invalid_flags)) {
      all_failed_algos <- c(all_failed_algos, ds[[column.algorithm]][1])
    }
  }
  all_failed_algos <- unique(all_failed_algos)
  if (length(all_failed_algos) > 0) {
    add_rows_failed <- lapply(all_failed_algos, function(algo) {
      alg_rows <- pp_data[pp_data$Algorithm == algo, ]
      if (nrow(alg_rows) == 0) return(NULL)
      # If only one unique Ratio and Fraction == 1, add pre (Fraction=0) point
      if (dplyr::n_distinct(alg_rows$Ratio) == 1 && max(alg_rows$Fraction) == 1 && min(alg_rows$Fraction) == 1) {
        data.frame(Algorithm = algo,
                   Ratio = alg_rows$Ratio[1],
                   Fraction = 0,
                   Transformed = 2,
                   stringsAsFactors = FALSE)
      } else {
        NULL
      }
    })
    add_rows_failed <- do.call(rbind, add_rows_failed)
    if (!is.null(add_rows_failed) && nrow(add_rows_failed) > 0) {
      pp_data <- dplyr::bind_rows(pp_data, add_rows_failed) %>%
        dplyr::arrange(Algorithm, Ratio, Fraction)
    }
  }
  # ### END ADD ###
  
  # Axis breaks & labels
  x_breaks <- c()
  x_labels <- c()
  offset <- 0
  from <- 1.0
  for (segment in segments) {
    stopifnot(length(segment$labels) == length(segment$breaks))
    min_value <- do.call(segment$trans, list(from))
    max_value <- do.call(segment$trans, list(segment$to))
    span <- max_value - min_value
    x_breaks <- c(x_breaks,
                  sapply(segment$breaks,
                         \(v) offset + segment$width *
                           (do.call(segment$trans, list(v)) - min_value) / span))
    x_labels <- c(x_labels,
                  if (tex) paste0("$", segment$labels, "$") else segment$labels)
    offset <- offset + segment$width
    from <- segment$to
  }
  x_breaks <- c(x_breaks, c(offset + 1/2, offset + 1))
  if (tex) {
    x_labels <- c(x_labels, c(label.tex.infeasible, label.tex.timeout))
  } else {
    x_labels <- c(x_labels, c(label.pdf.infeasible, label.pdf.timeout))
  }
  
  # Y labels
  y_labels <- if (tiny)
    c("0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0")
  else
    seq(0.0, 1.0, by = 0.1)
  
  if (length(factor.levels) > 0) {
    pp_data$Algorithm <- factor(pp_data$Algorithm, levels = factor.levels)
  }
  
  # Plot
  p <- ggplot(pp_data, aes(x = Ratio, y = Fraction, color = Algorithm, linetype = Algorithm)) +
    scale_x_continuous(expand = c(0, 0.01), breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0.01),
                       breaks = seq(0.0, 1.0, by = 0.1),
                       labels = y_labels) +
    geom_step(linewidth = 1.5)
  
  # Segment dividers
  x <- 0
  for (segment in segments) {
    x <- x + segment$width
    p <- p + geom_vline(xintercept = x)
  }
  
  # Colors
  if (length(colors) > 0) {
    p <- p + scale_color_manual(name = "Algorithm", values = colors)
  }
  
  # Distinct linetypes
  lt_vals <- c("solid", "longdash", "dotted", "dotdash", "twodash",
               "longdashdot", "123456", "42", "F1", "4313")
  p <- p + scale_linetype_manual(values = lt_vals[seq_along(levels(pp_data$Algorithm))],
                                 drop = TRUE)
  
  # Shapes removed/commented out as per previous decision
  
  return(p)
}
