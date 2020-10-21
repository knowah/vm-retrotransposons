df.as.gr <- function(df, ignore.NA=FALSE, use.names=FALSE, remove.invalid=TRUE) {
	# include a strand column if missing to simplify downstream steps
  if (!"strand" %in% names(df))
    df$strand <- factor("*", levels=c("+", "-", "*"))
  
  # identify seqnames column if missing
  if (!"seqnames" %in% names(df)) {
    if ("chrom" %in% names(df))
      names(df)[which(names(df)=="chrom")] <- "seqnames"
    else if ("chr" %in% names(df))
      names(df)[which(names(df)=="chr")] <- "seqnames"
    else
      stop("data frame contains no seqnames/chrom/chr column")
  }
  
  # check for start/end/pos information
  if (!"start" %in% names(df) || !"end" %in% names(df)) {
  	# if start and end missing, check for pos column
  	if (!"start" %in% names(df) && !"end" %in% names(df) && "pos" %in% names(df)) {
  		names(df)[which(names(df)=="pos")] <- "start"
  		df$end <- df$start
  	}
    stop("data frame missing start and/or end column")
  }
  
  if (ignore.NA) {
    bad.rows <- which(with(df, nchar(as.character(seqnames))==0 | 
                             is.na(seqnames) | 
                             is.na(start) |
                             is.na(end)))
    if (length(bad.rows) > 0) {
      df <- df[-bad.rows,]
    }
  }
  
  gr <- with(df, GRanges(as.character(seqnames), IRanges(start, end), strand=strand))
  
  if (ncol(df) > 4) {
    ignore.cols <- which(names(df) %in% c("seqnames", "start", "end", "strand"))
    if (remove.invalid)
    	ignore.cols <- c(ignore.cols, which(names(df) %in% c("ranges", "seqlevels", "seqlengths", "isCircular", "width", "element")))
    	
    mcols(gr) <- df[,-ignore.cols]
  }
  
  if (use.names)
  	names(gr) <- rownames(df)
  
  gr
}

ucsc2gr <- function(u) {
  col.split <- strsplit(u, ":")[[1]]
  coords <- as.integer(strsplit(gsub(",", "", col.split[2]),"-")[[1]])
  GRanges(col.split[1], IRanges(coords[1], coords[2]))
}

get.mm10.CpGs <- function(expected.file="data/R_objects/mm10.CpG.gr.RDS", write.if.missing=TRUE) {
	if (!file.exists(expected.file)) {
		if (!require("BSgenome.Mmusculus.UCSC.mm10")) {
			stop("Could not load mm10 CpGs; Bioconductor package 'BSgenome.Mmusculus.UCSC.mm10' not installed")
		}
		
		# identify CpG locations
		Mmusculus <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
		CGs <- Biostrings::vmatchPattern("CG", Mmusculus)
		
		# filter redundant (- strand) entries
		CGs <- CGs[strand(CGs)=="+"]
		strand(CGs) <- "*"
		
		if (write.if.missing) {
			saveRDS(CGs, expected.file)
		}
	} else {
		CGs <- readRDS(expected.file)
	}
	
	return(CGs)
}


extent.gr <- function(gr, by.col=NULL) {
	# assumption - all elements with the same grouping are on same chrom / strand
	if (is.null(by.col)) {
		return(GRanges(seqnames(gr)[1],
									 IRanges(min(start(gr)), max(end(gr))),
									 strand=strand(gr)[1]
		))
	}
	extent.df <- gr %>%
		as.data.frame %>%
		group_by_at(vars(by.col, seqnames, strand)) %>%
		summarize(st = min(start), en=max(end))
	new.gr <- with(extent.df, GRanges(seqnames, IRanges(st, en), strand=strand))
	mcols(new.gr)[,by.col] <- extent.df[,by.col]
	new.gr
}

gr.inside.flank.tiles <- function(gr, inside.tile.n, outside.flank.dist, outside.tile.width, keep.cols=NULL) {
	N.flank.tiles <- ceiling(outside.flank.dist / outside.tile.width)
	
	# clear existing seqlengths to simplify arithmetic (much faster implementation)
	# seqlengths will be restored at end of function
	old.seqlens <- seqlengths(gr)
	seqlengths(gr) <- NA
	
	# get tiles on the inside
	inside.tiles <- unlist(tile(gr, n=inside.tile.n))
	inside.tiles$pos <- unlist(lapply(
		as.character(strand(gr))=="-",
		function(minus) if (minus) (inside.tile.n-1):0 else 0:(inside.tile.n-1)
	))
	for (kc in keep.cols) {
		if (kc %in% colnames(mcols(gr))) {
			mcols(inside.tiles)[,kc] <- rep(mcols(gr)[,kc], each=inside.tile.n)
		}
	}
	
	five_prime.tiles <- unlist(
		tile(
			flank(gr, width=outside.flank.dist, start=TRUE, both=FALSE, ignore.strand=FALSE),
			width=outside.tile.width
		)
	)
	five_prime.tiles$pos <- unlist(lapply(
		as.character(strand(gr))=="-",
		function(minus) if (minus) -1:-N.flank.tiles else -N.flank.tiles:-1
	))
	for (kc in keep.cols) {
		if (kc %in% colnames(mcols(gr))) {
			mcols(five_prime.tiles)[,kc] <- rep(mcols(gr)[,kc], each=N.flank.tiles)
		}
	}
	
	three_prime.tiles <- unlist(
		tile(
			flank(gr, width=outside.flank.dist, start=FALSE, both=FALSE, ignore.strand=FALSE),
			width=outside.tile.width
		)
	)
	three_prime.tiles$pos <- unlist(lapply(
		as.character(strand(gr))=="-",
		function(minus)
			if (minus)
				(N.flank.tiles+inside.tile.n-1):inside.tile.n
			else 
				inside.tile.n:(N.flank.tiles+inside.tile.n-1)
	))
	for (kc in keep.cols) {
		if (kc %in% colnames(mcols(gr))) {
			mcols(three_prime.tiles)[,kc] <- rep(mcols(gr)[,kc], each=N.flank.tiles)
		}
	}
	
	combined.gr <- c(inside.tiles, five_prime.tiles, three_prime.tiles)
	
	# reinstate seqlengths and prune zero-width / out of bounds tiles
	suppressWarnings(seqlengths(combined.gr) <- old.seqlens)
	zero.width.tiles <- width(combined.gr) == 0 | end(combined.gr) < 1
	if (any(zero.width.tiles)) {
		combined.gr <- combined.gr[-which(zero.width.tiles)]
	}
	
	return(combined.gr)
}

score.features <- function(score.gr, feature.gr, score.col, score.fn="weighted.mean", new.score.col=score.col) {
	# get overlapping entries
	score.df <- as.data.frame(findOverlaps(score.gr, feature.gr))
	score.df$score <- mcols(score.gr[score.df$queryHits])[,score.col]
	
	if (is.function(score.fn)) {
		# apply scoring function to overlapping entries
		score.df <- score.df %>%
			group_by(subjectHits) %>%
			summarize(feat.score=score.fn(.data[[score.col]]))
	} else if (score.fn=="weighted.mean") {
		# get size of overlap between feature and score bin
		score.df$overlap.width <- width(
			pintersect(score.gr[score.df$queryHits], feature.gr[score.df$subjectHits])
		)
		
		# calculate weighted average of scores
		score.df <- score.df %>%
			mutate(wt.score=score*overlap.width) %>%
			group_by(subjectHits) %>%
			summarize(feat.score=sum(wt.score)/sum(overlap.width))
	} else {
		stop("score.fn must either be a function or 'weighted.mean'")
	}
	
	
	# add score data to original GR
	mcols(feature.gr)[,new.score.col] <- NA
	mcols(feature.gr)[score.df$subjectHits, new.score.col] <- score.df$feat.score
	
	return(feature.gr)
}


get.gene.plt <- function(db, gr, xlims=NULL, x.buffer=0.05, stagger=TRUE) {
	# obtain gene annotation data from database and convert to df
	gene.df <- biovizBase::crunch(db, gr)
	gene.df$component_id <- names(gene.df)
	names(gene.df) <- NULL
	gene.df <- filter(as.data.frame(gene.df), type %in% c("utr", "cds", "exon"))
	if (nrow(gene.df) == 0) return(ggplot() + theme_void())
	
	# some genes (e.g. miRNA) have exons but no CDS
	# all other genes will have 'exon' entries ignore
	# since CDS and UTR together make up the exons
	non.cds.genes <- gene.df %>%
		group_by(gene_name, type) %>%
		summarize(N=n()) %>%
		group_by(gene_name) %>%
		filter(n()==1 & type[1]=="exon") %>%
		pull(gene_name)
	
	# get df of CDS/UTR/exons to plot
	gene.plt.df <- gene.df %>%
		filter(gene_name %in% non.cds.genes | type != "exon") %>%
		dplyr::select(gene_name, type, seqnames, start, end, strand) %>%
		distinct()
	
	# make separate df with one entry per gene
	gene.plt.lines <- gene.plt.df %>%
		group_by(gene_name, seqnames, strand) %>%
		summarize(min.start=min(start),
							max.end=max(end)) %>%
		ungroup()
	
	# if limiting the x range, remove unnecessary exons / genes
	if (!is.null(xlims)) {
		gene.plt.df <- filter(gene.plt.df, end >= xlims[1] & start <= xlims[2])
		if (nrow(gene.df) == 0) return(ggplot() + theme_void())
		gene.plt.lines <- filter(gene.plt.lines, gene_name %in% unique(gene.plt.df$gene_name))
	} else {
		xlims <- c(start(gr), end(gr))
	}
	plt.width <- diff(xlims)
	
	# index the genes
	gene.plt.df$gene_idx <- as.numeric(factor(gene.plt.df$gene_name))
	gene.plt.lines <- left_join(gene.plt.lines, distinct(gene.plt.df[,c("gene_name", "gene_idx")]), by="gene_name")
	
	# define height of boxes (UTRs are shorter)
	gene.plt.df$h <- 0.3
	if ("utr" %in% gene.plt.df$type) {
		gene.plt.df[gene.plt.df$type=="utr",]$h <- 0.2
	}
	
	# set a display name for the genes including strand info
	gene.plt.lines$disp_name <- with(gene.plt.lines, 
																	 ifelse(
																	 	strand=="+", paste0(gene_name, " >"),
																	 	ifelse(
																	 		strand=="-", paste0("< ", gene_name),
																	 		# else (missing strand)
																	 		gene_name
																	 	)))
	
	# place genes on as many rows as necessary to avoid conflicts with nearby genes
	row.block.span <- round(plt.width * x.buffer) # buffer distance between adjacent genes
	gene.plt.gr <- with(gene.plt.lines, GRanges(seqnames, IRanges(min.start, max.end), gene_idx=gene_idx, gene_name=gene_name))
	gene.ovlps <- filter(as.data.frame(findOverlaps(gene.plt.gr, gene.plt.gr+row.block.span)), queryHits != subjectHits)
	g.rows <- c()
	for (g.idx in gene.plt.gr$gene_idx) {
		# identify already-placed genes in conflict with this one
		conflict <- with(gene.ovlps, subjectHits[queryHits==g.idx & subjectHits < g.idx])
		if (length(conflict) == 0) {
			g.rows <- c(g.rows, 1)
			next
		}
		blocked <- unique(g.rows[conflict]) # rows not to be used
		chosen <- min(setdiff(1:max(blocked+1), blocked))
		g.rows <- c(g.rows, chosen)
	}
	gene.plt.lines$gene_y <- g.rows
	gene.plt.df <- left_join(gene.plt.df, gene.plt.lines[,c("gene_idx", "gene_y")], by="gene_idx")
	
	# determine where to plot gene label text (x-axis)
	gene.plt.lines$text_x <- with(gene.plt.lines, (min.start+max.end)/2)
	xlims.text <- c(xlims[1] + plt.width * 0.05, xlims[2] - plt.width * 0.05)
	
	## identify off-screen labels
	midpt.off.left <- which(gene.plt.lines$text_x < xlims.text[1])
	midpt.off.right <- which(gene.plt.lines$text_x > xlims.text[2])
	midpt.off.both <- intersect(midpt.off.left, midpt.off.right)
	
	## fix off-screen labels
	if (length(midpt.off.both) > 0) {
		gene.plt.lines[midpt.off.both,]$text_x <- sum(xlims)/2
		midpt.off.left <- setdiff(midpt.off.left, midpt.off.both)
		midpt.off.right <- setdiff(midpt.off.right, midpt.off.both)
	}
	if (length(midpt.off.left) > 0) {
		gene.plt.lines[midpt.off.left,]$text_x <- gene.plt.lines[midpt.off.left,]$max.end - plt.width * 0.05
	}
	if (length(midpt.off.right) > 0) {
		gene.plt.lines[midpt.off.right,]$text_x <- gene.plt.lines[midpt.off.right,]$min.start + plt.width * 0.05
	}
	gene.plt.lines$text_x <- sapply(gene.plt.lines$text_x, function(x) max(x, xlims.text[1]))
	gene.plt.lines$text_x <- sapply(gene.plt.lines$text_x, function(x) min(x, xlims.text[2]))
	
	# stagger vertical position of gene label text (helps avoid collisions)
	if (stagger) {
		gene.plt.lines <- gene.plt.lines %>%
			arrange(gene_y, text_x) %>%
			group_by(gene_y) %>%
			mutate(text_y=0.125*(row_number() %% 2) + 0.375 + gene_y)
	} else {
		gene.plt.lines <- gene.plt.lines %>%
			mutate(text_y=gene_y+0.5)
	}
	
	
	# make plot
	ggplot() +
		# include horizontal line running entire length of gene
		geom_segment(data=gene.plt.lines, aes(x=min.start, y=gene_y, xend=max.end, yend=gene_y)) +
		# plot exon/cds/utr boxes
		geom_rect(data=gene.plt.df, aes(xmin=start, xmax=end, ymin=gene_y-h, ymax=gene_y+h), fill="black", color="black") +
		# plot label text
		geom_text(data=gene.plt.lines, aes(x=text_x, y=text_y, label=disp_name), hjust="center") +
		# limit x-axis to given limits...
		coord_cartesian(xlim=xlims) +
		# ... and prevent it being expanded
		scale_x_continuous(expand=c(0,0)) +
		# remove decoration
		theme_void() +
		theme(panel.grid.major.x = element_line(color="gray80"))
}

get.gr.tile.plot <- function(gr, color.col=NULL, limiting.gr=NULL, fill.colors=NULL, tile.border=0.1, legend.label=NULL) {
	if (length(gr) == 0) return(ggplot() + theme_void())
	
	if (!is.null(limiting.gr)) {
		stopifnot(length(limiting.gr) == 1)
		gr <- subsetByOverlaps(gr, limiting.gr)
		if (length(gr) == 0) return(ggplot() + theme_void())
		xlims <-c(start(limiting.gr)[1], end(limiting.gr)[1])
	} else {
		stopifnot(length(unique(seqnames(gr))) == 1)
		xlims <-c(min(start(gr)), max(end(limiting.gr)))
	}
	
	if (!is.null(color.col)) color.col <- sym(color.col)
	plt <- ggplot(as.data.frame(gr)) +
		geom_rect(
			aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=!!color.col, color=!!color.col),
			size=tile.border
		) +
		coord_cartesian(xlim=xlims, expand=FALSE) +
		theme_void() +
		theme(
			panel.grid.major.x = element_line(color="gray80"),
			legend.position="bottom",
			strip.text=element_text(vjust=1)
		)
	
	if (!is.null(fill.colors)) {
		plt <- plt +
			scale_fill_manual(name=legend.label, values=fill.colors) +
			scale_color_manual(name=legend.label, values=fill.colors)
	}
	
	plt
}


conv.chr.e2u <- function(gr) {
	seqlevelsStyle(gr) <- "UCSC"
	gr
}

conv.chr.u2e <- function(gr) {
	seqlevelsStyle(gr) <- "NCBI"
	gr
}
