# ChIP-seq Large Margin Supervised Penalty Learning Visualization
# Ported to animint2 for the gallery
# Source: https://github.com/tdhock/PeakSegFPOP-paper/blob/master/figure-large-margin.R

# Install required packages if needed:
# install.packages("animint2", dep=TRUE)
# install.packages("data.table")
# install.packages("quadmod", repos="http://R-Forge.R-project.org")
# remotes::install_github("tdhock/PeakError")
# remotes::install_github("tdhock/penaltyLearning")

library(data.table)
library(ggplot2)
library(animint2)
library(PeakError)
library(penaltyLearning)
library(quadmod)

### Download data files if they don't exist
set.name <- "H3K4me3_TDH_immune"
chunk.id <- 5

data.dir <- file.path("data", set.name, chunk.id)
dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)

data.files <- c("regions.RData", "counts.RData", "dp.model.RData")
base.url <- "https://rcdata.nau.edu/genomic-ml/chip-seq-chunk-db"

for(file.name in data.files) {
  local.file <- file.path(data.dir, file.name)
  if(!file.exists(local.file)) {
    remote.url <- paste0(base.url, "/", set.name, "/", chunk.id, "/", file.name)
    cat("Downloading", remote.url, "\n")
    download.file(remote.url, local.file)
  }
}

### Load data
chunk.name <- paste0(set.name, "/", chunk.id)
regions.file <- file.path(data.dir, "regions.RData")
load(regions.file)
counts.file <- file.path(data.dir, "counts.RData")
load(counts.file)
model.file <- file.path(data.dir, "dp.model.RData")
load(model.file)

region.list <- split(regions, regions$sample.id)
sample.ids <- c("McGill0002", "McGill0091", "McGill0322", "McGill0004")
chunk.counts <- data.table(counts)[sample.id %in% sample.ids,]
chunk.counts[, count := as.integer(coverage)]
chunk.counts[, sample.id := factor(sample.id, sample.ids)]
signal.list <- split(chunk.counts, chunk.counts$sample.id, drop=TRUE)

### Process model data
chunk.peak.list <- list()
chunk.region.list <- list()
chunk.loss.list <- list()

for(sample.id in sample.ids) {
  model.info <- dp.model[[sample.id]]
  chunk.loss.list[[sample.id]] <- data.table(
    sample.id=factor(sample.id, sample.ids), model.info$error)
  signal <- signal.list[[sample.id]]
  sample.regions <- region.list[[sample.id]]
  for(n.peaks in names(model.info$peaks)) {
    peak.df <- model.info$peaks[[n.peaks]]
    region.df <- PeakErrorChrom(peak.df, sample.regions)
    if(nrow(peak.df)) {
      chunk.peak.list[[paste(sample.id, n.peaks)]] <- data.table(
        sample.id=factor(sample.id, sample.ids), peaks=as.integer(n.peaks), peak.df)
    }
    chunk.region.list[[paste(sample.id, n.peaks)]] <- data.table(
      sample.id=factor(sample.id, sample.ids), peaks=as.integer(n.peaks), region.df)
  }
}

chunk.region <- do.call(rbind, chunk.region.list)
chunk.peak <- do.call(rbind, chunk.peak.list)
chunk.loss <- do.call(rbind, chunk.loss.list)

### Calculate error counts
error.counts <- chunk.region[, list(
  incorrect.labels=sum(fp+fn)
), by=.(sample.id, peaks)]

loss.list <- split(chunk.loss, chunk.loss$sample.id, drop=TRUE)
err.list <- split(error.counts, error.counts$sample.id, drop=TRUE)

exact.dfs.list <- list()
intervals.list <- list()

for(sample.id in names(loss.list)) {
  one <- loss.list[[sample.id]]
  # Use penaltyLearning::modelSelectionC for exact model selection
  exact.df <- penaltyLearning::modelSelectionC(one$error, one$peaks, one$peaks)
  err <- err.list[[sample.id]]
  setkey(err, peaks)
  signal <- signal.list[[sample.id]]
  # Match by peaks integer value
  exact.df$errors <- err[J(exact.df$model.id), ]$incorrect.labels
  # Find largest continuous minimum
  indices <- penaltyLearning::largestContinuousMinimumC(
    exact.df$errors, exact.df$max.log.lambda - exact.df$min.log.lambda)
  meta <- data.frame(sample.id, log.max.count=log(max(signal$coverage)))
  # Rename model.id to peaks for consistency
  exact.df$peaks <- exact.df$model.id
  exact.dfs.list[[sample.id]] <- data.table(meta, exact.df)
  intervals.list[[sample.id]] <- data.table(
    meta,
    min.log.lambda=exact.df$min.log.lambda[indices[1]],
    max.log.lambda=exact.df$max.log.lambda[indices[2]])
}

exact.dfs <- do.call(rbind, exact.dfs.list)
intervals <- do.call(rbind, intervals.list)

### Max margin function
max.margin <- function(features, limits, verbose=0, ...) {
  stopifnot(nrow(features)==nrow(limits))
  if(ncol(limits)!=2) stop("limits should be a 2-column matrix")
  stopifnot(is.matrix(features))
  
  has.limits <- apply(is.finite(limits), 1, any)
  some.limits <- limits[has.limits,]
  some.features <- features[has.limits,,drop=FALSE]
  
  scaled <- scale(some.features)
  mu <- attr(scaled,"scaled:center")
  sigma <- attr(scaled,"scaled:scale")
  
  n <- nrow(scaled)
  p <- ncol(scaled)
  vars <- make.ids(margin=1, intercept=1, weights=p)
  constraints <- list(vars$margin*1 >= 0)
  
  for(i in 1:n) {
    left <- some.limits[i,1]
    if(is.finite(left)) {
      ivars <- with(vars, {
        intercept * 1 + sum(weights)*scaled[i,] + margin*-1
      })
      constraints <- c(constraints, list(ivars >= left))
    }
    right <- some.limits[i,2]
    if(is.finite(right)) {
      ivars <- with(vars, {
        intercept * -1 + sum(weights)*scaled[i,]*-1 + margin*-1
      })
      constraints <- c(constraints, list(ivars >= -right))
    }
  }
  
  const.info <- standard.form.constraints(constraints, vars)
  n.vars <- length(unlist(vars))
  Dvec <- rep(1e-10, n.vars)
  D <- diag(Dvec)
  d <- rep(0, n.vars)
  d[vars$margin] <- 1
  sol <- quadprog::solve.QP(D, d, const.info$A, const.info$b0)
  
  sol$mu <- mu
  sol$sigma <- sigma
  sol$scaled <- scaled
  sol$log.limits <- some.limits
  sol$features <- some.features
  sol$weights <- sol$solution[vars$weights]
  sol$intercept <- sol$solution[vars$intercept]
  sol$margin <- sol$solution[vars$margin]
  sol$normalize <- function(X) {
    mu.mat <- matrix(mu, nrow(X), ncol(X), byrow=TRUE)
    s.mat <- matrix(sigma, nrow(X), ncol(X), byrow=TRUE)
    (X-mu.mat)/s.mat
  }
  sol$f <- function(x) sum(x*sol$weights) + sol$intercept
  sol$predict <- function(X) {
    stopifnot(is.matrix(X))
    X.norm <- sol$normalize(X)
    weights.mat <- matrix(sol$weights, nrow(X), ncol(X), byrow=TRUE)
    L.hat <- rowSums(X.norm * weights.mat) + sol$intercept
    L.hat
  }
  sol$L.pred <- apply(scaled, 1, sol$f)
  sol$lambda.pred <- sol$predict(features)
  sol
}

### Fit the model
fit <- with(intervals, {
  max.margin(cbind(log.max.count), cbind(min.log.lambda, max.log.lambda))
})

zero.error <- subset(exact.dfs, errors==0)
intervals$predicted <- fit$L.pred

### Prepare data for visualization
what.peaks <- "peaks"
what.error <- "errors"
exact.peaks <- data.frame(exact.dfs, what=what.peaks)
exact.error <- data.frame(exact.dfs, what=what.error)

zero.peaks <- subset(exact.peaks, errors==0)
z.error <- subset(exact.error, errors==0)

### Build duration list for per-sample peaks selection (1 selector per sample)
duration.list <- list()
duration.list[paste0(sample.ids, "peaks")] <- 2000

### Add per-sample peaks variable column (variable name per row)
exact.peaks$sample.peaks  <- paste0(exact.peaks$sample.id,  "peaks")
exact.error$sample.peaks  <- paste0(exact.error$sample.id,  "peaks")
exact.dfs$sample.peaks    <- paste0(exact.dfs$sample.id,    "peaks")
zero.peaks$sample.peaks   <- paste0(zero.peaks$sample.id,   "peaks")
z.error$sample.peaks      <- paste0(z.error$sample.id,      "peaks")
chunk.peak$sample.peaks   <- paste0(chunk.peak$sample.id,   "peaks")
chunk.region$sample.peaks <- paste0(chunk.region$sample.id, "peaks")

ann.colors <- c(noPeaks="#f6f4bf",
                peakStart="#ffafaf",
                peakEnd="#ff4c4c",
                peaks="#a445ee")

### Create display counts for smoother rendering
disp.counts <- chunk.counts[, {
  base <- seq(min(chromStart), max(chromEnd), l=500)
  data.table(base, count=approx(chromStart, coverage, base)$y)
}, by=sample.id]

### Prepare model comparison data
min.log.count <- min(intervals$log.max.count)
max.log.count <- max(intervals$log.max.count)
seg.dt <- data.table(
  min.log.lambda=9, min.log.count,
  max.log.lambda=c(11, 11.5), max.log.count)
seg.dt[, slope := (min.log.lambda-max.log.lambda)/(min.log.count-max.log.count)]
seg.dt[, intercept := min.log.lambda-slope*min.log.count]
reg.dt <- rbind(
  seg.dt[, .(slope, intercept)],
  data.frame(slope=0, intercept=8.5))
count.grid <- c(3, 7)
penalty.grid.list <- list()
for(reg.i in 1:nrow(reg.dt)) {
  r <- reg.dt[reg.i, ]
  penalty.grid.list[[reg.i]] <- data.table(
    count.grid, log.lambda=r$slope * count.grid + r$intercept, reg.i)
}
penalty.grid <- do.call(rbind, penalty.grid.list)

### Prepare interval data with predictions
setkey(intervals, sample.id)
exact.ids <- paste(exact.dfs$sample.id)
exact.dfs[, predicted := intervals[exact.ids,]$predicted]

### Helper function for plotting
notInf <- function(x) {
  min.log.pen <- min(exact.peaks$max.log.lambda)
  max.log.pen <- max(exact.peaks$min.log.lambda)
  ifelse(x==Inf, max.log.pen+1, ifelse(x==-Inf, min.log.pen-1, x))
}

### Create text labels
tsize <- 2.5
zero.error[, `:=`(
  label=sprintf("%d peak%s", peaks,
                ifelse(peaks==1, "", "s")),
  label.penalty=(min.log.lambda + max.log.lambda)/2)]
text.df <- data.frame(zero.error)
rownames(text.df) <- with(text.df, paste(sample.id, peaks))
text.df["McGill0004 2", "label.penalty"] <- 8.75
text.df["McGill0091 1", "label.penalty"] <- 10
text.df["McGill0002 2", "label.penalty"] <- 12

### Create the interactive visualization
viz <- animint(
  title="Max-margin supervised penalty learning for peak detection in ChIP-seq data",
  source="https://github.com/ANAMASGARD/chip-seq-large-margin-code/blob/main/chip-seq-large-margin.R",
  
  coverage=ggplot()+
    ggtitle("ChIP-seq data and peaks")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    theme_animint(width=800, height=500)+
    facet_grid(sample.id ~ ., scales="free")+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits) floor(limits[2]))+
    xlab("position on chr11 (kilo base pairs)")+
    animint2::geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                                fill=annotation),
                            data=chunk.region[peaks==0,],
                            color="grey",
                            alpha=0.5)+
    scale_linetype_manual("error type",
                          values=c(correct=0, "false negative"=3, "false positive"=1))+
    animint2::geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                                key=paste(sample.id, chromStart),
                                linetype=status),
                            showSelected=c("sample.peaks"="peaks"),
                            data=chunk.region,
                            color="black", fill=NA, size=1.5)+
    geom_line(aes(base/1e3, count),
              data=disp.counts, color="grey50")+
    geom_segment(aes(chromStart/1e3, 0, xend=chromEnd/1e3, yend=0,
                     key=paste(sample.id, chromStart, chromEnd)),
                 showSelected=c("sample.peaks"="peaks"),
                 data=chunk.peak, color="deepskyblue", size=4)+
    geom_point(aes(chromStart/1e3, 0,
                   key=paste(sample.id, chromStart)),
               showSelected=c("sample.peaks"="peaks"),
               data=chunk.peak, color="black", fill="deepskyblue")+
    scale_fill_manual("label", values=ann.colors, breaks=names(ann.colors)),
  
  penalty=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800, height=700)+
    facet_grid(sample.id ~ what, scales="free")+
    ggtitle("Select sample and number of peaks")+
    geom_label_aligned(aes(feature, log.penalty, label=label, vjust=vjust, hjust=hjust),
              data=data.table(
                what="regression",
                feature=5.8,
                log.penalty=c(12.5, 11, 9),
                hjust=c(0.5, 0, 0),
                label=c("0 errors\nlarge margin", "0 errors\nsmall margin", "1 error\nconstant"),
                vjust=c(0, 1, 0.5)),
              color="blue", size=3, fill=NA, label.size=0)+
    # Target intervals (clickable) - thick segment for easier clicking
    geom_segment(aes(log.max.count, min.log.lambda,
                     yend=max.log.lambda, xend=log.max.count),
                 clickSelects="sample.id",
                 data=data.table(intervals, what="regression"),
                 size=12, alpha=0.6, color="green")+
    # Zero error points
    geom_point(aes(log.max.count, min.log.lambda),
               clickSelects="sample.id",
               data=data.table(zero.error, what="regression")[is.finite(min.log.lambda),],
               size=tsize, color="grey")+
    geom_point(aes(log.max.count, min.log.lambda),
               clickSelects="sample.id",
               data=data.table(intervals, what="regression"),
               size=tsize, fill="white")+
    geom_point(aes(log.max.count, max.log.lambda),
               clickSelects="sample.id",
               data=data.table(intervals, what="regression"),
               size=tsize, fill="black")+
    # Sample labels (increased size for readability)
    geom_text(aes(log.max.count, max.log.lambda, label=sample.id,
                  hjust=ifelse(log.max.count==min(log.max.count), 0,
                               ifelse(log.max.count==max(log.max.count), 1, 0.5))),
              clickSelects="sample.id",
              data=data.table(intervals, what="regression"), vjust=-0.5, size=5)+
    # Margin line
    geom_segment(aes(log.max.count, min.log.lambda, yend=predicted, xend=log.max.count),
                 data=data.table(intervals["McGill0002",], what="regression"),
                 color="red", size=3)+
    # Model selection segments
    geom_segment(aes(log.max.count, notInf(min.log.lambda),
                     key=sample.id,
                     yend=notInf(max.log.lambda), xend=log.max.count),
                 showSelected=c("sample.peaks"="peaks"),
                 data=data.table(exact.dfs, what="regression"), size=1)+
    # Penalty lines
    geom_line(aes(count.grid, log.lambda, group=reg.i),
              data=data.table(penalty.grid, what="regression"),
              size=1, color="blue")+
    geom_text(aes(x, y, label=label),
              data=data.table(what="regression", x=5, y=6, label="log(max(coverage))"))+
    # Zero error segments (visible green bars)
    geom_segment(aes(peaks, min.log.lambda, yend=max.log.lambda, xend=peaks,
                     key=paste(sample.id, peaks)),
                 data=zero.peaks, size=4, color="green")+
    geom_segment(aes(errors, min.log.lambda, yend=max.log.lambda, xend=errors,
                     key=paste(sample.id, peaks)),
                 data=z.error, size=4, color="green")+
    # Model complexity segments
    geom_segment(aes(peaks, notInf(min.log.lambda), yend=notInf(max.log.lambda), xend=peaks,
                     key=peaks),
                 data=exact.peaks, size=2)+
    geom_segment(aes(errors, notInf(min.log.lambda), yend=notInf(max.log.lambda), xend=errors,
                     key=peaks),
                 data=exact.error, size=2)+
    # Clickable widerect for easier peak selection (one per sample row via faceting)
    animint2::geom_widerect(aes(ymin=notInf(min.log.lambda), ymax=notInf(max.log.lambda),
                                key=peaks),
                            clickSelects=c("sample.peaks"="peaks"),
                            alpha=0.2, data=exact.peaks)+
    animint2::geom_widerect(aes(ymin=notInf(min.log.lambda), ymax=notInf(max.log.lambda),
                                key=peaks),
                            clickSelects=c("sample.peaks"="peaks"),
                            alpha=0.2, data=exact.error)+
    scale_y_continuous("log(penalty)")+
    scale_x_continuous("", breaks=0:9),
  
  duration=duration.list,
  first=list(
    sample.id="McGill0002",
    McGill0002peaks=2,
    McGill0091peaks=1,
    McGill0322peaks=1,
    McGill0004peaks=2
  )
)

# Deploy to GitHub Pages
Sys.setenv(CHROMOTE_CHROME="/usr/bin/google-chrome")
animint2pages(viz, "chip-seq-large-margin", chromote_sleep_seconds=5)
