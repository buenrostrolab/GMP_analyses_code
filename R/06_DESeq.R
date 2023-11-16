### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University

library(dplyr)
library(DESeq2)
library(BuenColors)
library(BuenRTools)



setwd("<data_analysis_folder>")


# Load filtered SE
countsSE <- readRDS("./processed_results/countFiles/scATAC_SE_filt.rds")
cellMeta <- readRDS("./data/annot/cellMeta.rds")
stopifnot(all.equal(rownames(cellMeta),colnames(countsSE)))

table(countsSE$Sort,countsSE$Condition)
table(countsSE$SampleID,countsSE$Sort)
table(countsSE$SampleID,countsSE$Condition)

table(countsSE$cellSort,cellMeta$Louvain.h)

library(ggrepel)
labels.d <- cellMeta %>% group_by(Louvain.h) %>% summarise(UMAP1.h=median(UMAP1.h),UMAP2.h=median(UMAP2.h),cellType_GS=unique(cellType_GS)) 

gLouv <- shuf(cellMeta) %>% filter(cellSort %in% "GMP") %>% ggplot() + 
  geom_point_rast(aes(UMAP1.h,UMAP2.h,color=Louvain.h),size=0.1,shape=16) + theme_classic() + 
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  geom_label_repel(data=labels.d,aes(UMAP1.h,UMAP2.h,label=cellType_GS,color=Louvain.h),alpha=0.8) + 
  theme(legend.position = "none")

gLouv


cellMeta$cellType_GS <- gsub(x = cellMeta$cellType_GS,pattern=" ",replacement = "_")
table(cellMeta$cellType_GS,cellMeta$Louvain.h)

  
# Only keep non differentiated cells (GMP sorted)
myCellTypes <- c("Neu-biased_progenitors","Mono-biased_progenitors","Early_progenitors","HSC_and_early_progenitors","Early_GMP_progenitors") 
stopifnot(all(myCellTypes %in% cellMeta$cellType_GS))

# Make pseudobulks, per replicate, per condition, per cell type
counts.L <- list()

for(i in myCellTypes){
  cat("\n\nCell type: ",i,"\n")
  
  # Subset to celltype
  myCells <- cellMeta$cellType_GS %in% i
  
  for(j in unique(cellMeta$SampleID)){
    cond <- as.character(unique(cellMeta$Condition[cellMeta$SampleID %in% j]))
    cat("Sample: ",j,"\n")
    cat("Condition: ",cond,"\n\n")
    pseudoCells <- cellMeta$SampleID %in% j & myCells
    
    if(sum(pseudoCells) < 10){
      cat("Too few cells. Skipping ..\n")
      next 
      }
    counts.L[[paste(i,j,cond,sep="__")]] <- Matrix::rowSums(assay(countsSE)[,pseudoCells])
  }
}

names(counts.L)



pseudo.mat <- do.call('cbind',counts.L)

sampleData <- data.frame("Condition"=splitAndFetch(colnames(pseudo.mat),"__",3),
                         "CellType"=splitAndFetch(colnames(pseudo.mat),"__",1),
                         row.names = colnames(pseudo.mat))


dds <- DESeqDataSetFromMatrix(countData = pseudo.mat,
                              colData = sampleData,
                              design= ~ Condition + CellType) # Adjusting for cell type differences
dds$Condition <- factor(dds$Condition, levels = c("WT","Tet2KO","Tet2KO_Sox4OE"))
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Condition_Tet2KO_vs_WT")
dim(res)



res$Peak <- 1:nrow(res)
rownames(res) <- 1:nrow(res)


sum(res$padj <= 0.01,na.rm = TRUE)

sigPeaks <- which(res$padj <= 0.01)

# Score peaks (taking mean normalized counts and then Z-scoring)

tet2up <- rownames(res)[intersect(sigPeaks,which(res$log2FoldChange > 0))]
tet2dn <- rownames(res)[intersect(sigPeaks,which(res$log2FoldChange < 0))]
intersect(tet2up,tet2dn) # 0

length(tet2up)
length(tet2dn)


MAcols <- c("Up"="firebrick3","Down"="steelblue4","NS"="gray")

MAplotTet2 <- res %>% as.data.frame() %>%
  mutate("Group"=ifelse(rownames(res) %in% tet2up,"Up",ifelse(rownames(res) %in% tet2dn,"Down","NS"))) %>%
  ggplot(aes(x=log2(baseMean+1),y=log2FoldChange,color=Group)) + 
  geom_point_rast(size=0.2,shape=16) + scale_color_manual(values = MAcols)+
  theme_classic() + ggtitle("Tet2 vs WT") + 
  labs(subtitle = paste0("Up: ",length(tet2up),"\nDown: ",length(tet2dn)))+
  theme(legend.position = "none",plot.title = element_text(hjust=0.5),
        axis.text = element_text(color="black",size=6),
        axis.title = element_text(color="black",size=7)) 
MAplotTet2


siglist <- list("Tet2KO_vs_WT.up"=as.numeric(tet2up),"Tet2KO_vs_WT.dn"=as.numeric(tet2dn))


# Do the same for Sox4OE + Tet2KO vs WT
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Condition_Tet2KO_Sox4OE_vs_WT")

res$Peak <- 1:nrow(res)
rownames(res) <- 1:nrow(res)

sum(res$padj <= 0.01,na.rm = TRUE)

sigPeaks <- which(res$padj <= 0.01)


tet2soxup <- rownames(res)[intersect(sigPeaks,which(res$log2FoldChange > 0))]
tet2soxdn <- rownames(res)[intersect(sigPeaks,which(res$log2FoldChange < 0))]
intersect(tet2soxup,tet2soxdn) # 0

length(tet2soxup)
length(tet2soxdn)

MAplotTet2Sox4 <- res %>% as.data.frame() %>%
  mutate("Group"=ifelse(rownames(res) %in% tet2soxup,"Up",ifelse(rownames(res) %in% tet2soxdn,"Down","NS"))) %>%
  ggplot(aes(x=log2(baseMean+1),y=log2FoldChange,color=Group)) + 
  geom_point_rast(size=0.2,shape=16) + scale_color_manual(values = MAcols)+
  theme_classic() + ggtitle("Tet2_Sox4OE vs WT") + 
  labs(subtitle = paste0("Up: ",length(tet2soxup),"\nDown: ",length(tet2soxdn)))+
  theme(legend.position = "none",plot.title = element_text(hjust=0.5),
        axis.text = element_text(color="black",size=6),
        axis.title = element_text(color="black",size=7))

MAplotTet2Sox4



siglist$Tet2KOSox4OE_vs_WT.up <- as.numeric(tet2soxup)
siglist$Tet2KOSox4OE_vs_WT.dn <- as.numeric(tet2soxdn)
  
saveRDS(siglist,"./processed_results/DE/DESeq_tet2KO_vs_WT_GMPlike_sig.rds")

# Make annotation matrix
# We use this with chromVAR to derive Tet2 implicated accessibility scores
tet2_annot_ix <- Matrix(0,nrow=nrow(res),ncol = length(siglist),sparse = TRUE,dimnames = list(rownames(res),names(siglist)))
for(i in 1:length(siglist)){
  tet2_annot_ix[siglist[[i]],i] <- TRUE
}

head(tet2_annot_ix)

saveRDS(tet2_annot_ix,"./processed_results/DE/DESeq_tet2KO_vs_WT_GMPlike_annot_ix.rds")