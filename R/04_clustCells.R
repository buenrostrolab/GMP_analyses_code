### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University

library(dplyr)
library(cisTopic)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(patchwork)


cacheUMAP <- function(cachePath = NULL, dataSet, runAnyway = FALSE, seed = 123, 
          useUWOT = TRUE, ...){
  if (!is.null(cachePath)) {
    cat("Using provided cache root path: ", cachePath, "\n")
    R.cache::setCacheRootPath(cachePath)
  }
  else {
    R.cache::setCacheRootPath()
    cat("Setting cache root path as: ", R.cache::getCacheRootPath(), 
        "\n")
  }
  umapparams <- list(...)
  key <- list(dataSet = dataSet, seed = seed, UWOT = useUWOT, 
              umapparams = umapparams)
  umap <- R.cache::loadCache(key)
  if (is.null(umap) | runAnyway) {
    if (is.null(umap)) 
      cat("No existing UMAP run found for the given dataset and settings.. \n")
    if (runAnyway) 
      warning("runAnyway set to TRUE. May overwrite existing UMAP key-value pairs.. \n")
    if (!useUWOT) {
      validParams <- names(umap::umap.defaults)
      if (!is.null(umapparams)) 
        if (any(!names(umapparams) %in% validParams)) 
          stop(paste("Invalid UMAP parameter provided.. See following options and defaults: ", 
                     umap::umap.defaults, sep = "\n"))
      cat("Running UMAP on data..\n")
      custom.config = umap::umap.defaults
      if (!is.null(umapparams)) {
        custom.config[names(umapparams)] <- unlist(umapparams)
        int_vecs <- c("n_neighbors", "n_components", 
                      "min_dist", "local_connectivity", "bandwidth", 
                      "alpha", "gamma", "negative_sample_rate", "spread", 
                      "knn_repeats")
        custom.config[int_vecs] <- as.numeric(custom.config[int_vecs])
      }
      custom.config$random_state = seed
      umap <- umap::umap(dataSet, config = custom.config)
    }
    else {
      cat("Running uwot's UMAP implementation ..\n")
      set.seed(seed)
      umap <- uwot::umap(X = dataSet, ...)
    }
    cat("Saving to cache ..\n")
    R.cache::saveCache(umap, key = key)
  }
  else {
    cat("Found existing UMAP run for the given dataset and settings..\n")
  }
  cat("Finished!\n")
  return(umap)
}

setwd("<data_analysis_folder>")

# Load filtered SE
countsSE <- readRDS("./processed_results/countFiles/scATAC_SE_filt.rds")

# Load cisTopics output
cisOut <- readRDS("./processed_results/DR/cisTopics/scATAC_SE_filt_cisTopics.rds")

# Select model
cisOut <- selectModel(cisOut, type='maximum')

# Run UMAP based on selected output
nNeighbors <- 50
cisOut <- cisTopic::runUmap(cisOut, 
                  target='cell', seed=123, 
                  method='Z-score',n_neighbors=nNeighbors,
                  metric='cosine')

# Update so we don't have to re-run UMAP (takes a while)
saveRDS(cisOut,"./processed_results/DR/cisTopics/scATAC_SE_filt_cisTopics.rds")

# Originally, this won't have anything clustering related
cellMeta <- as.data.frame(colData(countsSE))

umap.d <- as.data.frame(cisOut@dr$cell$Umap)
colnames(umap.d) <- paste0("UMAP",1:2)

# Visualize by cluster / annotation
stopifnot(all.equal(rownames(umap.d),rownames(cellMeta)))

cellMeta <- cbind(umap.d,cellMeta)

cellMeta$Condition <- factor(cellMeta$Condition,levels=c("WT","Tet2KO","Tet2KO_Sox4OE"))

cisZ <- modelMatSelection(cisOut, 'cell', 'Z-score')

saveRDS(cisZ,"./processed_results/DR/cisTopics/scATAC_filt_cisTopics_cellZ.rds")

# Visualize

# Plot UMAPs
library(ggrastr)
gBase <- cellMeta %>% dplyr::select(UMAP1,UMAP2) %>% ggplot() + 
  geom_point_rast(aes(UMAP1,UMAP2),size=1,color="gray90") + 
  theme_classic()
gBase

# By condition
conditionCols <- c("Tet2KO"="mediumpurple","Tet2KO_Sox4OE"="darkorchid4","WT"="gray","WT_Sox4OE"="darkslategray")

gCondition <- shuf(cellMeta) %>% ggplot(aes(UMAP1,UMAP2,color=Condition)) + 
  geom_point_rast(size=0.5,shape=16) + theme_classic() + scale_color_manual(values = conditionCols) + 
  guides(colour = guide_legend(override.aes = list(size=3)))

gCondition


gConditionSplit <- gBase + geom_point(data=cellMeta,aes(UMAP1,UMAP2,color=Condition),size=0.5,shape=16) + 
  facet_wrap(~Condition,nrow=1) + scale_color_manual(values=conditionCols) + theme_classic() + 
  theme(strip.background = element_blank(),legend.position = "none",strip.text = element_text(size=12))
gConditionSplit


sortCols <- c("GMP"="darkorange","Lin_neg_Cd11b"="#00AF99")

gSortSplit <- gBase + geom_point(data=cellMeta,aes(UMAP1,UMAP2,color=cellSort),size=0.5,shape=16) + 
  facet_wrap(~cellSort,nrow=1) + scale_color_manual(values=sortCols) + theme_classic() + 
  theme(strip.background = element_blank(),legend.position = "none",strip.text = element_text(size=12))
gSortSplit


table(cellMeta$Condition,cellMeta$cellSort)

gRep <- gBase + geom_point(data=cellMeta,aes(UMAP1,UMAP2,color=SampleID),size=0.5,shape=16) + 
  facet_wrap(~SampleID,nrow=4) + theme_classic() + 
  theme(strip.background = element_blank(),legend.position = "none",strip.text = element_text(size=12))
gRep




############################# Harmony batch correction
# Batch correction is needed
# A+B, then C, and E as separate batches
countsSE$Batch <- ifelse(grepl(paste(c("^A","^B"),collapse = "|"),countsSE$SampleID),"Batch1",
                               ifelse(grepl("^C",countsSE$SampleID),"Batch2","Batch3"))
                                      
                                      

table(countsSE$SampleID,countsSE$Batch)

# Run Harmony
library(harmony)
set.seed(123)
cisZ.harmonized <- HarmonyMatrix(data_mat = t(cisZ), # Cistopics Z-score matrix input to harmony
                                 meta_data = as.data.frame(colData(countsSE)),
                                 vars_use = "Batch",
                                 do_pca = FALSE)
dim(cisZ.harmonized)
saveRDS(cisZ.harmonized,"./processed_results/DR/cisTopics/scATAC_filt_cisTopics_cellZ_harmony.rds")

# Re-cluster
umap.harmony <- cacheUMAP(cachePath = "./processed_results/UMAP_cache/",
                           dataSet = cisZ.harmonized,
                           metric="cosine",
                           n_neighbors=50)

colnames(umap.harmony) <- c("UMAP1.h","UMAP2.h")
umap.harmony <- cbind(umap.harmony,cellMeta)

gBase2 <- umap.harmony %>% dplyr::select(UMAP1.h,UMAP2.h) %>% ggplot() + 
  geom_point_rast(aes(UMAP1.h,UMAP2.h),size=1,color="gray90") + theme_classic()
gBase2

gRep2 <- gBase2 + geom_point(data=umap.harmony,aes(UMAP1.h,UMAP2.h,color=SampleID),size=0.5,shape=16) + 
  facet_wrap(~SampleID,nrow=4) + theme_classic() + 
  theme(strip.background = element_blank(),legend.position = "none",strip.text = element_text(size=12))
gRep2


# Louvain clustering with harmonized PCs

# Get new cell KNN
dim(cisZ.harmonized)
cellKNN.h <- FNN::get.knn(data = cisZ.harmonized,k = nNeighbors)$nn.index
dim(cellKNN.h)
rownames(cellKNN.h) <- rownames(cellMeta)

saveRDS(cellKNN.h,"./processed_results/DR/cisTopics/scATAC_filt_cisTopics_harmony_cellKNN.rds")

set.seed(123)
igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(cellKNN.h)[
  ,c("Var1", "value")]), directed=FALSE)), mode = "undirected")

# Louvain
clusters <- igraph::cluster_louvain(igraphObj)
Kmemberships <- as.numeric(igraph::membership(clusters))
table(Kmemberships)

umap.harmony$Louvain.h <- factor(Kmemberships)

labels.d <- umap.harmony %>% group_by(Louvain.h) %>% summarise(UMAP1=median(UMAP1.h),UMAP2=median(UMAP2.h))

library(ggrepel)
gLouv <- shuf(umap.harmony) %>% ggplot() + 
  geom_point_rast(aes(UMAP1.h,UMAP2.h,color=Louvain.h),size=0.1,shape=16) + theme_classic() + 
  guides(colour = guide_legend(override.aes = list(size=1.2))) + 
  geom_label_repel(data=labels.d,aes(UMAP1,UMAP2,label=Louvain.h,color=Louvain.h),size=4) + 
  theme(legend.position = "none")

gLouv


# Update cell meta (use this henceforth since it has new UMAP and clustering)
all.equal(umap.harmony[,3:(ncol(umap.harmony)-1)],cellMeta)
saveRDS(umap.harmony,"./data/annot/cellMeta.rds")



############################## Visualize Motifs
# See 05_runchromVAR.R
motif_dev <- readRDS("./processed_results/chromVAR/scATAC_filt_motif_dev.rds")
motifZ <- chromVAR::deviationScores(motif_dev)

rownames(motifZ) <- extractTFNames(rownames(motifZ))

stopifnot(all.equal(colnames(motifZ),rownames(cellMeta)))
myMotifs <- c("Sox4","Spi1","Cebpe","Irf8","Gata1","Klf4","Ets2","Bcl11a","Runx1")

# Smooth
library(doParallel)
motifZ.s <- smoothGeneScoresNN(NNmat = cellKNN.h,TSSmat = motifZ,geneList = myMotifs,nCores = 4)

gMotif <- plotMarker2D(df = umap.harmony,
             markers = myMotifs,
             markerMat = motifZ.s,
             rasteRize = TRUE,
             pointSize=0.01,
             minCutoff = -3,maxCutoff = 3,
             combine = FALSE)

cowplot::plot_grid(plotlist = gMotif,align="hv")


for(i in 1:length(gMotif)){
  gMotif[[i]]
}


# Load cluster annotations (celltype labels)
labelAnnot <- read.csv("./data/annot/bioRad_in_vivo_annotation_Gschiroli.csv",header=FALSE,row.names = 1,stringsAsFactors = FALSE)

stopifnot(all(rownames(labelAnnot) %in% umap.harmony$Louvain.h))

umap.harmony$cellType_GS <- labelAnnot[umap.harmony$Louvain.h,"V2"] 

library(RColorBrewer)
cellTypeCols <- brewer.pal(name = "Set3",n=6)
cellTypeCols <- c(cellTypeCols,brewer.pal(name = "Paired",n=7))
names(cellTypeCols) <- levels(as.factor(umap.harmony$cellType_GS))

labels.d <- umap.harmony %>% group_by(cellType_GS) %>% summarise(UMAP1=median(UMAP1.h),UMAP2=median(UMAP2.h))

gAnnot <- shuf(umap.harmony) %>% ggplot() + 
  geom_point_rast(aes(UMAP1.h,UMAP2.h,color=cellType_GS),size=0.1,shape=16) + theme_classic() + 
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  geom_label_repel(data=labels.d,aes(UMAP1,UMAP2,label=cellType_GS,color=cellType_GS),
                   color="black",min.segment.length = 0.1,size=3) + 
  scale_color_manual(values=cellTypeCols)  +
  theme(legend.position = "none")

gAnnot

# Make version without labels, and with legend

gAnnot2 <- shuf(umap.harmony) %>% ggplot() + 
  geom_point_rast(aes(UMAP1.h,UMAP2.h,color=cellType_GS),size=0.1,shape=16) + theme_classic() + 
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  scale_color_manual(values=cellTypeCols) + labs(color="Annotation")

gAnnot2



# Re-save cellMeta data finally, to include cell type annotation
saveRDS(umap.harmony,"./data/annot/cellMeta.rds")

# Also export PDFs of cells highlighting condition / cell sort / sample

gCondition2 <- shuf(umap.harmony) %>% filter(Condition %in% c("WT","Tet2KO")) %>% 
  ggplot(aes(UMAP1.h,UMAP2.h,color=Condition)) + 
  geom_point_rast(size=0.1,shape=16) + theme_classic() + scale_color_manual(values = conditionCols) + 
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  labs(color="Genotype")

gCondition2



gSort2 <- shuf(umap.harmony) %>% filter(Condition %in% c("WT","Tet2KO")) %>%  
  ggplot(aes(UMAP1.h,UMAP2.h,color=Sort)) + 
  geom_point_rast(size=0.1,shape=16) +
  scale_color_manual(values=sortCols) + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  labs(color="Sort")

gSort2



#Compute % cell KNNs that are Tet2 KO vs WT


cellKNN <- FNN::get.knn(cisZ.harmonized[!umap.harmony$Condition %in% "Tet2KO_Sox4OE",],k = 50)$nn.index
rownames(cellKNN) <- rownames(umap.harmony)[!umap.harmony$Condition %in% "Tet2KO_Sox4OE"]


cellMeta.sub <- umap.harmony[!umap.harmony$Condition %in% "Tet2KO_Sox4OE",]
tet2.pt <- apply(cellKNN,1,function(x) {t(table(cellMeta.sub[x,"Condition"])/50)})
dim(tet2.pt)
tet2.pt[,1:5]

cellMeta.sub$tet2KOpt <- tet2.pt[2,]

gtet2per <- cellMeta.sub %>% ggplot(aes(UMAP1.h,UMAP2.h,color=tet2KOpt)) + 
  geom_point_rast(size=0.05,shape=16) + scale_color_gradientn(colours = jdb_palette("wolfgang_basic")) +
  theme_classic() + labs(color="% Tet2 KO\nneighbors")

gCondition2 + gtet2per

ggcustom(gtet2per,clean = TRUE,legend.key.size = unit(0.3,"cm"),splitLeg = TRUE)



# Same, but with only Tet2KO + Sox4OE and WT data (no TET2KO)
# Re-derive KNN without Tet2KO_Sox4OE data

cellKNN <- FNN::get.knn(cisZ.harmonized[!umap.harmony$Condition %in% "Tet2KO",],k = 50)$nn.index
rownames(cellKNN) <- rownames(umap.harmony)[!umap.harmony$Condition %in% "Tet2KO"]


cellMeta.sub <- umap.harmony[!umap.harmony$Condition %in% "Tet2KO",]
table(cellMeta.sub$Condition)
levels(cellMeta.sub$Condition)

tet2sox4.pt <- apply(cellKNN,1,function(x) {t(table(cellMeta.sub[x,"Condition"])/50)})
dim(tet2sox4.pt)

tet2sox4.pt[1:3,1:5]
levels(cellMeta.sub$Condition) # 3rd row is Tet2KO Sox4OE percentage

cellMeta.sub$tet2KOsox4pt <- tet2sox4.pt[3,]

gtet2sox4per <- cellMeta.sub %>% ggplot(aes(UMAP1.h,UMAP2.h,color=tet2KOsox4pt)) + 
  geom_point_rast(size=0.05,shape=16) + scale_color_gradientn(colours = jdb_palette("wolfgang_basic")) +
  theme_classic() + labs(color="% Tet2 KO + Sox4 OE\nneighbors")

gtet2sox4per

gtet2per + gtet2sox4per