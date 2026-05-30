#' @title Plot enrichment map
#' @description
#' Plot enrichment map through a vector or matrix of scores and a self-defined
#' set that summarizes a few groups of the names of the vector or matrix.
#'
#' @name enrichment
#'
#' @import ggplot2
#' @importFrom stats p.adjust
#'
#' @param x a vector or matrix of scores to be enriched
#' @param custom.set a self-defined set. The first column contains feature names,
#'   and the second column, preferably named "group", contains group names.
#' @param alpha exponent weight of the score of ordered features.
#' @param normalize logical value to determine if normalizing enrichment scores.
#' @param permute.n number of custom-set permutations for significance testing.
#' @param padj.method correction method passed to stats::p.adjust.
#' @param pvalue.cutoff cutoff for both unadjusted and adjusted p-value.
#' @param angle angle of rotating x-axis labels.
#' @param match.feature one of "rownames" or "colnames". Default is "rownames",
#'   which keeps the original EnrichIntersect behavior.
#' @param ... other arguments
#'
#' @return A list including S, pvalue, and g.
#'
#' @export
enrichment <- function(
    x,
    custom.set,
    alpha = 0,
    normalize = TRUE,
    permute.n = 100,
    padj.method = "none",
    pvalue.cutoff = 0.05,
    angle = 45,
    match.feature = c("rownames", "colnames"),
    ...
) {
  match.feature <- match.arg(match.feature)
  
  ## --------------------------------------------------------------------------
  ## Original-style input handling
  ## --------------------------------------------------------------------------
  
  if (is.matrix(x) || is.data.frame(x)) {
    if (any(colSums(is.na(x)) == nrow(x)) && ncol(x) > 1) {
      stop("The argument 'x' matrix has some columns with all missing values!")
    }
  } else {
    x <- as.matrix(x)
  }
  
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  if (is.null(colnames(x))) {
    colnames(x) <- "X"
  }
  
  if (is.null(rownames(x))) {
    stop("The argument 'x' must have rownames matching the feature names in 'custom.set'!")
  }
  
  if (is.matrix(custom.set) || is.data.frame(custom.set)) {
    if (dim(custom.set)[2] != 2) {
      stop("The argument 'custom.set' has to have two columns!")
    }
  } else {
    stop("The argument 'custom.set' has to be a matrix or dataframe!")
  }
  
  custom.set <- as.data.frame(custom.set, stringsAsFactors = FALSE)
  
  ## Original code uses custom.set$group.
  ## If the second column is not named "group", name it "group" for compatibility.
  if (!"group" %in% colnames(custom.set)) {
    colnames(custom.set)[2] <- "group"
  }
  
  custom.set[[1]] <- trimws(as.character(custom.set[[1]]))
  custom.set$group <- trimws(as.character(custom.set$group))
  
  ## --------------------------------------------------------------------------
  ## Support original orientation and optional transposed orientation
  ##
  ## Original behavior:
  ##   features are rownames(x)
  ##
  ## Optional:
  ##   features are colnames(x), so use t(x) internally
  ## --------------------------------------------------------------------------
  
  if (match.feature == "colnames") {
    if (is.null(colnames(x))) {
      stop("When match.feature = 'colnames', x must have colnames.")
    }
    
    x_work <- t(x)
  } else {
    x_work <- x
  }
  
  if (anyDuplicated(rownames(x_work))) {
    stop("Feature names must be unique.")
  }
  
  features <- intersect(rownames(x_work), custom.set[[1]])
  
  if (length(features) == 0L) {
    stop("None of the features in 'custom.set' matched the feature names of 'x'.")
  }
  
  groups <- unique(custom.set$group[custom.set[[1]] %in% features])
  n_groups <- length(groups)
  
  if (n_groups == 0L) {
    stop("No valid groups were found in 'custom.set'.")
  }
  
  ## Restrict x to matched features, preserving intersect() order.
  x_work <- x_work[features, , drop = FALSE]
  
  ## Build group-wise feature index list for the C++ core.
  set_indices <- lapply(groups, function(g) {
    idx <- match(
      custom.set[[1]][custom.set$group == g],
      features,
      nomatch = 0L
    )
    
    as.integer(unique(idx[idx > 0L]))
  })
  
  names(set_indices) <- groups
  
  ## --------------------------------------------------------------------------
  ## RcppArmadillo replacement for original heavy for-loop
  ##
  ## C++ returns matrices with:
  ##   rows    = columns/profiles of x_work
  ##   columns = custom-set groups
  ##
  ## This matches the original S and pvalue orientation before t().
  ## --------------------------------------------------------------------------
  
  core <- enrichment_core_original(
    x = x_work,
    set_indices = set_indices,
    alpha = alpha,
    normalize = normalize,
    permute_n = as.integer(permute.n)
  )
  
  S <- core$S
  pvalue <- core$pvalue
  
  rownames(S) <- colnames(x_work)
  colnames(S) <- groups
  
  rownames(pvalue) <- colnames(x_work)
  colnames(pvalue) <- groups
  
  ## --------------------------------------------------------------------------
  ## Original post-processing
  ## --------------------------------------------------------------------------
  
  pvalue[is.na(S)] <- NA
  pvalue <- t(pvalue)
  
  if (!padj.method == "none") {
    pvalue <- p.adjust(pvalue, method = padj.method)
  }
  
  S <- t(S)
  
  ## Define plotting dataframe
  dat <- data.frame(
    x = factor(rep(colnames(S), each = nrow(S))),
    y = rep(rownames(S), ncol(S)),
    ks = as.vector(S),
    pvalue = as.vector(pvalue)
  )
  
  dat$y <- factor(
    dat$y,
    levels = levels(factor(dat$y))[c(length(unique(dat$y)):1)]
  )
  
  ks.min <- min(dat$ks, na.rm = TRUE)
  
  dat$ks[dat$ks < 0 & !is.na(dat$ks)] <- ks.min - 0.001
  dat$ks <- dat$ks - (ks.min - 0.001)
  
  if (sum(!is.na(dat$ks)) >= 2L) {
    dat$ks[which.max(dat$ks)] <- round(sort(dat$ks, decreasing = TRUE)[2] + 0.5)
  }
  
  dat$border <- rep("red", nrow(dat))
  dat$border[dat$pvalue >= pvalue.cutoff] <- "gray"
  
  y <- border <- ks <- NULL
  
  if (normalize) {
    ES.name <- "Normalized\nEnrichment\nScore"
  } else {
    ES.name <- "Enrichment\nScore"
  }
  
  ## --------------------------------------------------------------------------
  ## Original plot format, with legend order:
  ##
  ##   1. significance legend, e.g. p<0.05
  ##   2. p-value colorbar
  ##   3. enrichment score size legend
  ## --------------------------------------------------------------------------
  
  if (any(dat$pvalue < pvalue.cutoff, na.rm = TRUE)) {
    g <- ggplot(data = dat) +
      geom_point(
        aes(x = x, y = y, color = border, fill = pvalue, size = ks),
        shape = 21
      ) +
      scale_fill_gradientn(
        name = "p-value",
        na.value = "black",
        colours = c("blue", "white"),
        limits = c(0, 1),
        guide = guide_colorbar(
          barheight = 3,
          barwidth = 1,
          order = 2
        )
      ) +
      geom_point(
        aes(x = x, y = y, size = ks),
        color = dat$border,
        shape = 21
      ) +
      guides(
        colour = guide_legend(
          override.aes = list(size = 5),
          order = 1
        )
      ) +
      scale_size(
        name = ES.name,
        range = c(1, 5),
        breaks = c(0, 1, 2, 3),
        guide = guide_legend(
          keyheight = .8,
          order = 3
        )
      ) +
      theme(
        axis.text.x = element_text(
          size = 8,
          angle = angle,
          vjust = 1,
          hjust = 1
        ),
        legend.margin = margin(-0.1, 0, 0, 0, unit = "cm"),
        legend.box = "vertical"
      ) +
      xlab("") +
      ylab("")
    
    if (pvalue.cutoff == 0.05) {
      g <- g +
        scale_color_manual(
          name = NULL,
          values = c(`red` = "red"),
          labels = "p<0.05",
          guide = guide_legend(order = 1)
        )
    } else {
      if (pvalue.cutoff == 0.1) {
        g <- g +
          scale_color_manual(
            name = NULL,
            values = c(`red` = "red"),
            labels = "p<0.1",
            guide = guide_legend(order = 1)
          )
      } else {
        g <- g +
          scale_color_manual(
            name = NULL,
            values = c("gray", "red"),
            labels = paste0("p", c(">=", "<"), pvalue.cutoff),
            guide = guide_legend(order = 1)
          )
      }
    }
  } else {
    g <- ggplot(data = dat) +
      geom_point(
        aes(x = x, y = y, fill = pvalue, size = ks),
        shape = 21
      ) +
      scale_fill_gradientn(
        name = "p-value",
        na.value = "black",
        colours = c("blue", "white"),
        limits = c(round(min(dat$pvalue, na.rm = TRUE), 2), 1),
        guide = guide_colorbar(
          barheight = 3,
          barwidth = 1,
          order = 1
        )
      ) +
      scale_size(
        name = ES.name,
        range = c(1, 5),
        breaks = c(0, 1, 2, 3),
        guide = guide_legend(
          keyheight = .8,
          order = 2
        )
      ) +
      theme(
        axis.text.x = element_text(
          size = 8,
          angle = 45,
          vjust = 1,
          hjust = 1
        ),
        legend.box = "vertical"
      ) +
      xlab("") +
      ylab("")
  }
  
  print(g)
  
  return(
    list(
      S = S,
      pvalue = pvalue,
      g = g
    )
  )
}