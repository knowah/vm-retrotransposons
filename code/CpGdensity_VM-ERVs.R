library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(purrr)
library(effsize)
#upload ERV annotations (only ERV families that are represented in VM-ERV and false positive lists) and cat into single GRange
ERV1 = read_tsv("/data/genomes/mm10/annotation/rmskOutCurrent.Dfam_2_0.v4_0_7.ERV1.fixed_100kb_s_20_r_20.tsv",col_names = FALSE)
ERV1 = with(ERV1, GRanges(X1, IRanges(start=X2, end=X3), strand=X4, type=X5, family=X7, ID=X10))
ERVK = read_tsv("/data/genomes/mm10/annotation/rmskOutCurrent.Dfam_2_0.v4_0_7.ERVK.fixed_100kb_s_20_r_20.tsv",col_names = FALSE)
ERVK = with(ERVK, GRanges(X1, IRanges(start=X2, end=X3), strand=X4, type=X5, family=X7, ID=X10))
ERVL = read_tsv("/data/genomes/mm10/annotation/rmskOutCurrent.Dfam_2_0.v4_0_7.ERVL.fixed_100kb_s_20_r_20.tsv",col_names = FALSE)
ERVL = with(ERVL, GRanges(X1, IRanges(start=X2, end=X3), strand=X4, type=X5, family=X7, ID=X10))
ERVL_MaLR = read_tsv("/data/genomes/mm10/annotation/rmskOutCurrent.Dfam_2_0.v4_0_7.ERVL-MaLR.fixed_100kb_s_20_r_20.tsv",col_names = FALSE)
ERVL_MaLR = with(ERVL_MaLR, GRanges(X1, IRanges(start=X2, end=X3), strand=X4, type=X5, family=X7, ID=X10))
ERV_combine <- c(ERV1, ERVK, ERVL, ERVL_MaLR)

#upload VM-ERV and False positive beds
VM_ERV <- read_tsv("/data/reference/VM-ERVs/VM-ERVs_coord.bed", col_names = FALSE)
VM_ERV = with(VM_ERV, GRanges(X1, IRanges(start=X2, end=X3)))
#FP_ERV <- read_tsv("/data/reference/VM-ERVs/falsePos_coord.bed", col_names = FALSE)
#FP_ERV = with(FP_ERV, GRanges(X1, IRanges(start=X2, end=X3)))
cVM_IAP <- read_tsv("/data/reference/ME/IAP_validation.July2019.True_ME.meta_sub.5pLTR.bed", col_names = FALSE) 
cVM_IAP = with(cVM_IAP, GRanges(X1, IRanges(start=X2, end=X3)))
tsVM_IAP = read_tsv("/data/reference/ME/IAP_validation.July2019.Tissue_spec.meta_sub.5pLTR.bed", col_names = FALSE)
tsVM_IAP = with(tsVM_IAP, GRanges(X1, IRanges(start=X2, end=X3)))

#add VM-ERV and false positive information to the ERV GRange
VM_ERV_overlapCombine <- findOverlaps(VM_ERV, ERV_combine)
#FP_ERV_overlapCombine <- findOverlaps(FP_ERV, ERV_combine)
cVM_IAP_overlapCombine <- findOverlaps(cVM_IAP, ERV_combine)
tsVM_IAP_overlapCombine <- findOverlaps(tsVM_IAP, ERV_combine)
ERV_combine$validations <- NA
ERV_combine[subjectHits(VM_ERV_overlapCombine)]$validations <- "VM-ERVs"
#ERV_combine[subjectHits(FP_ERV_overlapCombine)]$validations <- "False Positives"
ERV_combine[subjectHits(cVM_IAP_overlapCombine)]$validations <- "cVM-IAPs"
ERV_combine[subjectHits(tsVM_IAP_overlapCombine)]$validations <- "tsVM-IAPs"


#upload mm10 cpgs
mm10CpG = readRDS("/data/genomes/mm10/annotation/mm10.CpG.gr.RDS")

#calculate and add cpg density information to ERV GRange
ERV_combine.cpgOverlaps <- as.data.frame(table(queryHits(findOverlaps(ERV_combine, mm10CpG))),stringsAsFactors=FALSE)
ERV_combine$cpgCounts <- 0
ERV_combine[as.numeric(ERV_combine.cpgOverlaps$Var1)]$cpgCounts <- ERV_combine.cpgOverlaps$Freq
ERV_combine$norm <- (ERV_combine$cpgCounts / width(ERV_combine))*100

#add alpha info for colouring plot
ERV_combine$alpha <- ifelse(is.na(ERV_combine$validations), 0, 1)
ERV_combine$validations[is.na(ERV_combine$validations)] <- " "

#filter ERV GRange to LTRs 
ERV_combine_fixed <- ERV_combine[which(!grepl(".*[-int|_I]$", ERV_combine$type)&width(ERV_combine)<700&width(ERV_combine)>200)]
list_of_relevant_types <- unique(ERV_combine_fixed[which(ERV_combine_fixed$validations=="VM-ERVs"|ERV_combine_fixed$validations=="cVM-IAPs"|ERV_combine_fixed$validations=="tsVM-IAPs")]$type)
list_of_relevant_families <- unique(ERV_combine_fixed[which(ERV_combine_fixed$validations=="VM-ERVs"|ERV_combine_fixed$validations=="cVM-IAPs"|ERV_combine_fixed$validations=="tsVM-IAPs")]$family)
#take out one LTR
ERV_combine_fixed_5pLTR<-as.data.frame(ERV_combine_fixed) %>% group_by(ID) %>% arrange(ifelse(strand=="+",start,-start)) %>% filter(row_number()==1) %>% ungroup() %>% arrange(ID)

#plot of types
ggplot(ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$type %in% list_of_relevant_types,], aes(x=type, y=norm)) + 
  geom_boxplot() +
  geom_jitter(aes(color=validations, alpha=alpha), shape=21, size=2, position=position_jitter(seed=1000)) +
  geom_jitter(aes(color=validations, alpha=alpha), shape=20, size=1.2, position=position_jitter(seed=1000)) +
  scale_alpha_identity() +
  scale_color_manual(values=c("cVM-IAPs"="#469C8D","tsVM-IAPs"="#BF7F41"," "=NA,"VM-ERVs"="blue")) +
  guides(color=guide_legend(title=NULL)) +
  xlab("LTR type") +
  ylab("CpG dinucleotides per 100bp") +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black"),
        axis.ticks.x=element_line(color="black"))

#plot of families
ggplot(ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$family %in% list_of_relevant_families,], aes(x=family, y=norm)) + 
  geom_boxplot() +
  geom_jitter(aes(color=validations, alpha=alpha), shape=21, size=2, position=position_jitter(seed=1000)) +
  geom_jitter(aes(color=validations, alpha=alpha), shape=20, size=1.2, position=position_jitter(seed=1000)) +
  scale_alpha_identity() +
  scale_color_manual(values=c("cVM-IAPs"="#469C8D","tsVM-IAPs"="#BF7F41"," "=NA,"VM-ERVs"="blue")) +
  guides(color=guide_legend(title=NULL)) +
  xlab("ERV family") +
  ylab("CpG dinucleotides per 100bp") +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black"),
        axis.ticks.x=element_line(color="black"))

#ttest for VM-ERVs,cVM-IAPs, and tsVM-IAPs compared to all ERVs
VM_ERV.ttest<-t.test(ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$validations=="VM-ERVs",]$norm,ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$validations==" ",]$norm)
cVM_IAP.ttest<-t.test(ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$validations=="cVM-IAPs",]$norm,ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$validations==" ",]$norm)
tsVM_IAP.ttest<-t.test(ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$validations=="tsVM-IAPs",]$norm,ERV_combine_fixed_5pLTR[ERV_combine_fixed_5pLTR$validations==" ",]$norm)


