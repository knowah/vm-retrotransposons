identify.LTR.components <- function(ltr.seq, ccatcc.tata.dist=40, tata.r_g.dist.seed=30, tata.r_g.dist.max=10, tata.polyA.dist.range=c(40, 100), polyA.min.start.dist=100, polyA.CA.dist.range=c(10, 70)) {
	
	# identify possible polyA tails and downstream CA site (marks end of R region)
	poss.polyA.sites <- matchPattern("AATAAA", ltr.seq, max.mismatch=0)
	if (length(poss.polyA.sites) == 0) {
		poss.polyA.sites <- matchPattern("AATAAA", ltr.seq, max.mismatch=1)
	}
	poss.polyA.sites <- poss.polyA.sites[start(poss.polyA.sites) >= polyA.min.start.dist]
	CA.sites <- matchPattern("CA", ltr.seq)
	poss.polyA_CA.pairs <- list()
	for (p.pA in start(poss.polyA.sites)) {
		any.ca <- FALSE
		for (ca in start(CA.sites)) {
			if (ca >= p.pA + polyA.CA.dist.range[1] && ca <= p.pA + polyA.CA.dist.range[2]) {
				any.ca <- TRUE
				poss.polyA_CA.pairs <- c(poss.polyA_CA.pairs, list(c(p.pA, ca)))
				break # only take the first acceptable CA hit
			}
		}
		if (!any.ca) {
			poss.polyA_CA.pairs <- c(poss.polyA_CA.pairs, list(c(p.pA, NA)))
		}
	}
	poss.polyA_CA.pairs.pA_st <- purrr::map(poss.polyA_CA.pairs, 1)
	
	# identify possible CCAATC sites
	poss.CCAATC.sites <- matchPattern("CCAATC", ltr.seq, max.mismatch=1, with.indels=TRUE)
	
	valid.total.matches <- list()
	if (length(poss.CCAATC.sites) > 0) {
		# find possible TATA boxes downstreaam of CCAATC site
		poss.CCAATC.3p_flank <- DNAStringSet(flank(poss.CCAATC.sites, width=ccatcc.tata.dist, start=FALSE, both=FALSE))
		TATA.in.CCAATC <- vmatchPattern("ATATAA", poss.CCAATC.3p_flank)
		tata.hits <- which(unlist(lapply(TATA.in.CCAATC, length)) > 0)
		
		# if no canonical seq, try AGGATAA
		if (length(tata.hits) == 0) {
			TATA.in.CCAATC <- vmatchPattern("AGGATAA", poss.CCAATC.3p_flank)
			tata.hits <- which(unlist(lapply(TATA.in.CCAATC, length)) > 0)
		}
		
		# try AGGATAA
		if (length(tata.hits) == 0) {
			TATA.in.CCAATC <- vmatchPattern("AGGAGAA", poss.CCAATC.3p_flank)
			tata.hits <- which(unlist(lapply(TATA.in.CCAATC, length)) > 0)
		}
		
		# next try TATAA (missing initial A)
		if (length(tata.hits) == 0) {
			TATA.in.CCAATC <- vmatchPattern("TATAA", poss.CCAATC.3p_flank)
			tata.hits <- which(unlist(lapply(TATA.in.CCAATC, length)) > 0)
		}
		
		if (length(tata.hits) > 0 && length(poss.polyA_CA.pairs) > 0) {
			for (t.hit in tata.hits) {
				# get absolute position of TATA box in LTR
				tata.ltr.starts <- end(poss.CCAATC.sites)[t.hit] + start(TATA.in.CCAATC[[t.hit]])
				tata.ltr.ends <- end(poss.CCAATC.sites)[t.hit] + end(TATA.in.CCAATC[[t.hit]])
				
				for (t.l.idx in length(tata.ltr.starts)) {
					t.l.st <- tata.ltr.starts[t.l.idx]
					# find R-initiating 'G' site downstream of TATA box
					r_g.pos <- NA
					g.guess.loc <- t.l.st + tata.r_g.dist.seed
					if (g.guess.loc >= length(ltr.seq)) {
						next
					}
					g.search.seq <- ltr.seq[(g.guess.loc-tata.r_g.dist.max):min(nchar(ltr.seq), g.guess.loc+tata.r_g.dist.max)]
					# also verify r_g dist max is actually respected! 
					tata.r_init_g.dists <- start(matchPattern("G", g.search.seq))-(tata.r_g.dist.max+1)
					if (length(tata.r_init_g.dists) > 0) {
						best.r_g <- tata.r_init_g.dists[which(abs(tata.r_init_g.dists)==min(abs(tata.r_init_g.dists)))[1]]
						r_g.pos <- g.guess.loc + best.r_g
					}
					
					valid.pA_CA.idx <- integer(0)
					# find polyA/CA sites within reasonable distance of TATA box
					if (length(poss.polyA_CA.pairs) > 1) {
						valid.pA_CA.idx <- which(poss.polyA_CA.pairs.pA_st >= t.l.st + tata.polyA.dist.range[1] & poss.polyA_CA.pairs.pA_st <= t.l.st + tata.polyA.dist.range[2])
					} else if (length(poss.polyA_CA.pairs) == 1 && poss.polyA_CA.pairs.pA_st >= t.l.st + tata.polyA.dist.range[1]) {
						valid.pA_CA.idx <- 1
					}
					for (v.i in valid.pA_CA.idx) {
						valid.total.matches <- c(
							valid.total.matches, 
							list(list(
								"CCAATC" = c(start(poss.CCAATC.sites[t.hit]), end(poss.CCAATC.sites[t.hit])),
								"TATA"   = c(t.l.st, tata.ltr.ends[t.l.idx]),
								"R.initial.G" = r_g.pos,
								"polyA"  = c(poss.polyA_CA.pairs[[v.i]][1], poss.polyA_CA.pairs[[v.i]][1]+6),
								"R.terminal.CA" = c(poss.polyA_CA.pairs[[v.i]][2], poss.polyA_CA.pairs[[v.i]][2]+1)
							))
						)
					}
				}
			}
		}
	}
	
	return(valid.total.matches)
}

##choose.top.component.breakdown <- function(comps, ltr.size, ideal.r.size.frac=c(0.1, 0.5)) {
#	if (length(comps) == 0) {
#		# return blank
#	}
#	
#	unique.boundaries <- unique(lapply(comps, function(x) list(x$R.initial.G, x$R.terminal.CA[2])))
#	r.sizes <- lapply(unique.boundaries, function(x) (x[2]-x[1]+1) / ltr.size)
#	ideal.r <- which(r.sizes >= ideal.r.size.frac[1] & r.sizes <= ideal.r.size.frac[2])
#	
#	if (length(ideal.r) == 1) {
#		# only one ideal R
#		
#	} else if (length(ideal.r) == 0) {
#		# no ideal R's
#		if (length(comps) == 1) {
#			# use this
#		} else {
#			biggest.small.r <- 
#		}
#	} else {
#		# ideal.r > 1
#	}
#}
#
