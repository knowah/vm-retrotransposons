library(tidyr)
library(dplyr)
library(msir)
library(FourCSeq)
library(cowplot)

existing.viewpoint.frags <- function(proj.path=".")
{
  frag.file <- file.path(proj.path, "restriction_frag_info", "primerFragments.txt")
  if (!file.exists(frag.file)) return(character(0))
  return(read.table(frag.file, sep="\t", header=TRUE, as.is=TRUE)$viewpoint)
}

add.fragments <- function(object, minSize=20, filter=TRUE, force=FALSE)
{
  stopifnot(class(object) == "FourC")
  
  gr.file = file.path(metadata(object)$fragmentDir,
    sprintf("valid_fragments.%s_%s.minSize_%d.%sfilter.RDS",
            metadata(object)$reSequence1,
            metadata(object)$reSequence2,
            minSize,
            ifelse(filter, "", "no_")
    )
  )
  if (force || !file.exists(gr.file))
  {
    getFragments <- getFromNamespace("getFragments", "FourCSeq")
    getSites <- getFromNamespace("getSites", "FourCSeq")
    ref = getReferenceSeq(object)
    frag = getFragments(metadata(object)$reSequence1, ref)
    site = getSites(metadata(object)$reSequence2, ref)
    fragWithSite = subsetByOverlaps(frag, site, minoverlap = nchar(metadata(object)$reSequence2))
    firstSites = site[findOverlaps(fragWithSite, site, select = "first")]
    lastSites = site[findOverlaps(fragWithSite, site, select = "last")]
    mcols(fragWithSite)$leftSize = start(firstSites) - start(fragWithSite)
    mcols(fragWithSite)$rightSize = end(fragWithSite) - end(lastSites)
    mcols(fragWithSite)$leftSize[start(firstSites) < start(fragWithSite)] = -1
    mcols(fragWithSite)$rightSize[end(fragWithSite) < end(lastSites)] = -1
    mcols(fragWithSite)$leftValid = mcols(fragWithSite)$leftSize >= 
      minSize
    mcols(fragWithSite)$rightValid = mcols(fragWithSite)$rightSize >= 
      minSize
    fragNoSite = setdiff(frag, fragWithSite)
    mcols(fragNoSite)$leftSize = -1
    mcols(fragNoSite)$rightSize = -1
    mcols(fragNoSite)$leftValid = FALSE
    mcols(fragNoSite)$rightValid = FALSE
    frag = sort(c(fragNoSite, fragWithSite))
    if (filter) 
      frag = frag[frag$leftValid | frag$rightValid]
    mcolsRows <- DataFrame(type = rep("referenceFragments", ncol(mcols(frag))), 
                           description = rep("", ncol(mcols(frag))))
    mcols(mcols(frag)) <- mcolsRows
    if (!dir.exists(metadata(object)$fragmentDir)) {
    	dir.create(metadata(object)$fragmentDir)
    }
    saveRDS(frag, gr.file)
  }
  else
  {
    message("  [loading from ", gr.file, "]")
    frag <- readRDS(gr.file)
  }

  rowRanges(object) <- frag
  object
}

get.meta.and.col.data <- function(
	vp.name, ref.genome, cutters, proj.path="data/4C-seq",
	primer.fa=file.path(proj.path, "primers.fa"), bam.path=file.path(proj.path, "bams"), seq.pr="first"
) {
  metadata <- list(
  	projectPath = proj.path,
    fragmentDir = file.path(proj.path, "restriction_frag_info"),
    referenceGenomeFile = ref.genome,
    reSequence1 = cutters[1],
    reSequence2 = cutters[2],
    primerFile = primer.fa,
    bamFilePath = bam.path
  )
  
  bams <- basename(Sys.glob(file.path(metadata$bamFilePath, paste0("4C_*.",vp.name,".bam"))))
  if (length(bams) == 0) {
  	stop("No matching bam files in ", metadata$bamFilePath, ".")	
  }
  
  indivs <- gsub("^4C_", "", unlist(purrr::map(strsplit(bams, "\\."),1)))
  
  colData <- DataFrame(
  	viewpoint = vp.name,
    condition = factor(indivs),
    replicate = 1,
    bamFile = bams,
    sequencingPrimer=seq.pr
  )
  
  return(
    list(
      "metadata" = metadata,
      "colData" = colData
    )
  )
}

smoothed.counts.fc <- function(col.d, meta.d, trim=0, minMapq=30, force.find=FALSE) 
{
  message("Getting smoothed counts from FourCSeq")
  vp.name <- col.d$viewpoint[1]
  fc <- FourC(col.d, meta.d)
  message("  Finding genomic fragments")
  fc <- add.fragments(fc)
  
  if (force.find || !vp.name %in% existing.viewpoint.frags(meta.d$projectPath))
  {
    message("  Generating viewpoint fragments for ", vp.name, "...")
    findViewpointFragments(fc)
  }
  else
  {
    message("  Using existing fragment dataset for ", vp.name)
  }
  message("  Adding viewpoint frags")
  fc <- addViewpointFrags(fc)
  
  message("  Counting frag overlaps")
  fc <- countFragmentOverlaps(fc, trim=trim, minMapq=minMapq)
  message("  Combining frag ends & smoothing counts")
  fc <- combineFragEnds(fc)
  fc <- smoothCounts(fc)
  
  fc
}

subset.fc.by.vp.chrom <- function(fc)
{
  vp.chrom <- colData(fc)$chr[1]
  fc[seqnames(fc)==vp.chrom,]
}

trim.frag.data <- function(fc, minCount=40) 
{
  fc <- subset.fc.by.vp.chrom(fc)
  fragData = getDistAroundVp(colData(fc)$viewpoint[1], colData(fc), rowRanges(fc))
  medianCounts <- apply(counts(fc), 1, median)
  toLeft <- fragData$dist > -20000 & fragData$dist < 0 & 
    !is.na(fragData$dist)
  afterMin <- 1:sum(toLeft) > tail(which(sign(diff(medianCounts[toLeft])) < 
                                           0), 1) + 1
  toExclude <- which(toLeft)[afterMin]
  toRight <- fragData$dist < 20000 & fragData$dist > 0 & 
    !is.na(fragData$dist)
  beforeMin <- 1:sum(toRight) < which(sign(diff(medianCounts[toRight])) > 
                                        0)[1]
  toExclude <- c(toExclude, which(toRight)[beforeMin])
  toExclude = union(toExclude, which(abs(fragData$dist) < 
                                       1000))
  tooClose = rep(FALSE, length(fragData$dist))
  tooClose[toExclude] = TRUE
  
  fragData$posLeft[tooClose] = FALSE
  fragData$posRight[tooClose] = FALSE
  fragData$tooClose = tooClose
  lowCounts = medianCounts < minCount
  fragData$posLeft[lowCounts] = FALSE
  fragData$posRight[lowCounts] = FALSE
  fragData$lowCounts = lowCounts
  fragData$selectedForFit <- (fragData$posLeft | fragData$posRight)
  
  fragData
}

trim.frag.data.simple <- function(fc, minCount=1, ignoreInner=10000, maxDist=2000000)
{
  fc <- subset.fc.by.vp.chrom(fc)
  fragData <- getDistAroundVp(colData(fc)$viewpoint[1], colData(fc), rowRanges(fc))
  medianCounts <- apply(counts(fc), 1, median)
  lowCounts <- medianCounts < minCount
  not.selected <- medianCounts < minCount | abs(fragData$dist) < ignoreInner | abs(fragData$dist) >= maxDist
  fragData$selectedForFit <- TRUE
  fragData[not.selected,]$selectedForFit <- FALSE
  fragData
}

log.norm.counts <- function(fc, frag.trim, norm.type="max") 
{
  vals <- counts(fc[frag.trim$selectedForFit,])
  if ("sum" %in% norm.type) 
    vals <- apply(vals, 2, function(x) x / sum(x))
  if ("max" %in% norm.type)
    vals <- apply(vals, 2, function(x) {log((x+1) / max(x+1))})
  if ("median" %in% norm.type)
    vals <- apply(cals, 2, function(x) {log((x+1) / median(x+1))})
  return(vals)
}

get.plot.df <- function(vp.name, count.data, frag.trim)
{
  data.cols <- which(startsWith(colnames(count.data), paste0(vp.name, "_")))
  trim.df.long <- pivot_longer(
    bind_cols(
      as.data.frame(count.data),
      as.data.frame(frag.trim[frag.trim$selectedForFit])[,c("start", "end", "mid", "dist")]
    ),
    cols=data.cols,
    names_to="sample",
    values_to="counts"
  )
  trim.df.long$sample <- factor(trim.df.long$sample, levels=colnames(count.data)[data.cols])
  trim.df.long
}

get.loess.and.sd <- function(df, var.dep, var.indep, alpha=0.25, sigma=1.96)
{
  l.sd <- loess.sd(df[[var.dep]] ~ df[[var.indep]], span=alpha, nsigma=sigma)
  L_SD <- data.frame(x=as.numeric(l.sd$x),
                     y=as.numeric(l.sd$y),
                     ylow=as.numeric(l.sd$lower),
                     yupp=as.numeric(l.sd$upper)
  ) %>% distinct()
  L_SD
}

get.loess.df <- function(plot.df, ...)
{
  trim.loess.df <- bind_rows(
    sapply(
      levels(plot.df$sample),
      function(samp) get.loess.and.sd(filter(plot.df, sample==samp), "counts", "mid", ...),
      simplify=FALSE,
      USE.NAMES=TRUE
    ),
    .id="sample"
  )
  trim.loess.df$sample <- factor(trim.loess.df$sample, levels=levels(plot.df$sample))
  
  trim.loess.df
}

plot.l_sd <- function(L_SD) 
{
  c(
    geom_ribbon(aes(x=x, ymin=ylow, ymax=yupp, group=sample), data=L_SD, alpha=0.2, fill="salmon"),
    #geom_line(aes(x=x, y=ylow, group=sample), data=L_SD, alpha=0.3, color="red"),
    #geom_line(aes(x=x, y=yupp, group=sample), data=L_SD, alpha=0.3, color="red"),
    geom_line(aes(x=x, y=y, group=sample), data=L_SD, size=1, color="blue")
  )
}

plot.loess_color <- function(L_SD) 
{
  geom_line(aes(x=x, y=y, group=sample, color=sample), data=L_SD, size=1)
}

make.indiv.plot <- function(plot.df, loess.df, vp.coords, x_span=500000, y_max="counts", clip=FALSE) 
{
  if (y_max=="counts")
  {
    ylims <- c(min(plot.df$counts)-0.1, max(plot.df$counts)+0.1)
  }
  else if (y_max=="loess")
  {
    ylims <- c(min(loess.df$y) * 0.98, max(loess.df$y)*2)
  }
  xlims <- mean(vp.coords)+c(-1,1)*x_span
  if (clip)
  {
    plot.df <- plot.df[plot.df$mid > xlims[1] & plot.df$mid < xlims[2],]
    loess.df <- loess.df[loess.df$x > xlims[1] & loess.df$x < xlims[2],]
  }
  ggplot(plot.df) +
    geom_point(aes(x=mid, y=counts), alpha=0.3) +
    facet_wrap(vars(sample), ncol=1) +
    theme_minimal() +
    geom_vline(xintercept=c(vp.coords[1], vp.coords[2]), linetype="dashed") +
    coord_cartesian(xlim=xlims, ylim=ylims, expand=FALSE) +
    plot.l_sd(loess.df)
}

make.combined.plot <- function(plot.df, loess.df, vp.coords, x_span=500000, y_max="counts", clip=FALSE) 
{
  if (y_max=="counts")
  {
    ylims <- c(min(plot.df$counts)-0.1, max(plot.df$counts)+0.1)
  }
  else if (y_max=="loess")
  {
    ylims <- c(min(loess.df$y) * 0.98, max(loess.df$y)*2)
  }
  xlims <- mean(vp.coords)+c(-1,1)*x_span
  if (clip)
  {
    plot.df <- plot.df[plot.df$mid > xlims[1] & plot.df$mid < xlims[2],]
    loess.df <- loess.df[loess.df$x > xlims[1] & loess.df$x < xlims[2],]
  }
  ggplot(plot.df) +
    geom_point(aes(x=mid, y=counts, color=sample), alpha=0.15) +
    theme_minimal() +
    geom_vline(xintercept=c(vp.coords[1], vp.coords[2]), linetype="dashed") +
    coord_cartesian(xlim=xlims, ylim=ylims, expand=FALSE) +
    plot.loess_color(loess.df)
}

get.vp.coords <- function(fc) 
{
  c(colData(fc)$start[1], colData(fc)$end[1])
}

run.analysis <- function(vp.name, data.path, genome.fa, cutters)
{
  fc.meta <- get.meta.and.col.data(vp.name, data.path, genome.fa, cutters)
  fc <- smoothed.counts.fc(fc.meta$colData, fc.meta$metadata)
  message("Subsetting data using FourCSeq method")
  fc <- subset.fc.by.vp.chrom(fc)
  frag.trim <- trim.frag.data(fc)
  message("Calculating log-norm counts")
  lognorm <- log.norm.counts(fc, frag.trim)
  
  lognorm.plt.df <- get.plot.df(vp.name, lognorm, frag.trim)
  lognorm.loess.df <- get.loess.df(lognorm.plt.df)
  
  vp.coords <- get.vp.coords(fc)
  
  message("Generating plots")
  vp.chrom <- colData(fc)$chr[1]
  indiv.plt <- make.indiv.plot(lognorm.plt.df, lognorm.loess.df, vp.coords) 
  + xlab(vp.chrom) + ylab("Normalized fragment counts")
  combo.plt <- make.combined.plot(lognorm.plt.df, lognorm.loess.df, vp.coords) + xlab(vp.chrom) + ylab("Normalized fragment counts")
  message("Done.")
  
  return(
    list(
      fc = fc,
      frag.trim = frag.trim,
      lognorm.counts = lognorm,
      indiv.plt = indiv.plt,
      combo.plt = combo.plt
    )
  )
}

arrange.combo.and.gene.plt <- function(fourc.obj_list) {
	gene.y.ranges <- layer_scales(fourc.obj_list$gene.plt)$y$range$range
	gene.rows <- if (!is.null(gene.y.ranges)) floor(gene.y.ranges[2]) else 0
	plt.heights <- c(0.4, 0.3, 0.3)
	if(gene.rows < 4) {
		dummy.height <- plt.heights[3]/2
		dummy.plt <- ggplot() + theme_void()
		ap <- align_plots(fourc.obj_list$combo.plt + ylab(NULL) + theme(legend.position="none", axis.text.y=element_blank()),
											fourc.obj_list$hmm.plt, 
											fourc.obj_list$gene.plt,
											dummy.plt,
											align="v", axis="lr")
		plt.heights <- c(plt.heights[1:2], dummy.height, dummy.height)
	}
	else {
		ap <- align_plots(fourc.obj_list$combo.plt + ylab(NULL) + theme(legend.position="none", axis.text.y=element_blank()),
											fourc.obj_list$hmm.plt, 
											fourc.obj_list$gene.plt,
											align="v", axis="lr")
	}
	gridExtra::grid.arrange(grobs=ap, ncol=1, heights=plt.heights)
}

arrange.indiv.and.gene.plt <- function(fourc.obj_list) {
	gene.y.ranges <- layer_scales(fourc.obj_list$gene.plt)$y$range$range
	gene.rows <- if (!is.null(gene.y.ranges)) floor(gene.y.ranges[2]) else 0
	plt.heights <- c(0.6, 0.3, 0.3)
	if(gene.rows < 4) {
		dummy.height <- plt.heights[3]/2
		dummy.plt <- ggplot() + theme_void()
		ap <- align_plots(fourc.obj_list$indiv.plt, #+ ylab(NULL) + theme(legend.position="none", axis.text.y=element_blank()),
											fourc.obj_list$hmm.plt, 
											fourc.obj_list$gene.plt,
											dummy.plt,
											align="v", axis="lr")
		plt.heights <- c(plt.heights[1:2], dummy.height, dummy.height)
	}
	else {
		ap <- align_plots(fourc.obj_list$indiv.plt, #+ ylab(NULL) + theme(legend.position="none", axis.text.y=element_blank()),
											fourc.obj_list$hmm.plt, 
											fourc.obj_list$gene.plt,
											align="v", axis="lr")
	}
	gridExtra::grid.arrange(grobs=ap, ncol=1, heights=plt.heights)
	#ap <- align_plots(fourc.obj_list$indiv.plt + ylab(NULL) + theme(legend.position="none", axis.text.y=element_blank()),
	#                  fourc.obj_list$hmm.plt, 
	#                  fourc.obj_list$gene.plt,
	#                  align="v", axis="lr")
	#gridExtra::grid.arrange(grobs=ap, ncol=1, heights=c(0.6, 0.3, 0.3))
}

refactor_col <- function(old.values, names.df) {
	new.order <- names.df[sort(match(unique(old.values), names.df$old.name)),]$new.name
	new.values <- names.df[match(old.values, names.df$old.name),]$new.name
	return(factor(new.values, levels=new.order))
}

refactored_df <- function(df, col.to.refactor, names.df) {
	df[,col.to.refactor] <- refactor_col(df[[col.to.refactor]], names.df)
	return(df)
}

process_4C_data <- function(vp.name, loess.alpha, ensdb, chromhmm, ignoreInner=1000, window.size=400000) {
	message("Generating 4C data for viewpoint: ", vp.name)
	# get fragment counts
	fc.meta <- get.meta.and.col.data(vp.name, GENOME.4C, CUTTERS)
	fc <- smoothed.counts.fc(fc.meta$colData, fc.meta$metadata)
	
	# subset and trim data
	message(vp.name, ": Subsetting data using FourCSeq method")
	fc <- subset.fc.by.vp.chrom(fc)
	#frag.trim <- trim.frag.data(fc)
	frag.trim <- trim.frag.data.simple(fc, ignoreInner = ignoreInner)
	
	# calculate lognorm counts and get data for plots
	message(vp.name, ": log-normalising count data")
	lognorm.sum <- log.norm.counts(fc, frag.trim, "sum")
	lognorm.plt.df <- get.plot.df(vp.name, lognorm.sum, frag.trim)
	lognorm.loess.df <- get.loess.df(lognorm.plt.df, alpha=loess.alpha)
	
	vp.coords <- get.vp.coords(fc)
	vp.chrom <- colData(fc)$chr[1]
	Xlims <- mean(vp.coords) + c(-1, 1) * (window.size/2)
	region.gr <- GRanges(sprintf("%s:%d-%d", vp.chrom, Xlims[1], Xlims[2]))
	
	# generate plots
	message(vp.name, ": generating plots")
	indiv.plt <- make.indiv.plot(refactored_df(lognorm.plt.df,"sample",fourc.methy), refactored_df(lognorm.loess.df,"sample",fourc.methy), vp.coords, x_span=200000, y_max="loess", clip=TRUE) + xlab(vp.chrom) + ylab("Normalized fragment counts")
	combo.plt <- make.combined.plot(refactored_df(lognorm.plt.df,"sample",fourc.methy), refactored_df(lognorm.loess.df,"sample",fourc.methy), vp.coords, x_span=200000, y_max="loess", clip=TRUE) + xlab(vp.chrom) + ylab("Normalized fragment counts")
	gene.plt <- get.gene.plt(ensdb, conv.chr.u2e(region.gr), x.buffer=0.1, stagger=FALSE)
	hmm.plt <- get.gr.tile.plot(chromhmm, "state", region.gr, chromHMM.colors, legend.label="ChromHMM states") + facet_wrap(~ category, ncol=1, drop=FALSE)
	
	return(
		list(
			fc = fc,
			frag.trim = frag.trim,
			plt.df = lognorm.plt.df,
			loess.df = lognorm.loess.df,
			lognorm.counts = lognorm.sum,
			indiv.plt = indiv.plt,
			combo.plt = combo.plt,
			gene.plt = gene.plt,
			hmm.plt = hmm.plt
		)
	)
}


