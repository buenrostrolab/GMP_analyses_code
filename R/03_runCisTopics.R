### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University


library(cisTopic)
library(SummarizedExperiment)
library(Matrix)

cisTopicFromSE <- function(ATAC.se, assay = NULL, projectName = NULL, ...) {
  peakRanges <- rowRanges(ATAC.se)
  rownames(ATAC.se) <- paste(as.character(seqnames(peakRanges)), 
                             paste(start(peakRanges), end(peakRanges), sep = "-"), 
                             sep = ":")
  if (is.null(projectName)) {
    projectName <- "cisTopicProject"
  }
  if (is.null(assay)) {
    peakMat <- assay(ATAC.se)
  }
  else {
    peakMat <- assays(ATAC.se)[[assay]]
  }
  cisObj <- createcisTopicObject(count.matrix = peakMat, project.name = projectName,...)
  
  if (!is.null(colData(ATAC.se))) 
    cisObj <- addCellMetadata(cisObj, cell.data = as.data.frame(colData(ATAC.se)[cisObj@cell.names,]))
  if (!is.null(rowData(ATAC.se))) 
    cisObj <- addRegionMetadata(cisObj, region.data = as.data.frame(rowData(ATAC.se)[cisObj@region.names,]))
  cisObj
}

setwd("<data_analysis_folder>")

# Load filtered SE
countsSE <- readRDS("./processed_results/countFiles/scATAC_SE_filt.rds")

# Convert to cisTopics object
countsCis <- cisTopicFromSE(countsSE)

# Run cisTopics (http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/WarpLDA_10X_workflow.html)
cat("Running cisTopics ..\n")
countsCis <- runWarpLDAModels(countsCis, topic=c(50), seed=987, nCores=4, iterations = 500, addModels=FALSE)

cat("Saving output cisTopics object  ..\n")
saveRDS(countsCis,"./processed_results/DR/cisTopics/scATAC_SE_filt_cisTopics.rds")
cat("Finished! \n")

# For peak x topic annotation matrix
cisOut <- getRegionsScores(cisOut)

cisOut <- binarizecisTopics(
  cisOut,
  method = "GammaFit",
  thrP = 0.995,
  plot = FALSE,
  cutoffs = NULL
)

sapply(cisOut@binarized.cisTopics,nrow)

nTopics <- length(cisOut@binarized.cisTopics)

# Peak x Topic binary annotation matrix
topic_annot_ix <- Matrix(0,
                         ncol=nTopics,
                         nrow=nrow(countsSE))

rownames(topic_annot_ix) <- Signac::GRangesToString(granges(countsSE),sep = c(":","-"))

# Just making sure all the peak names are identical as in the cistopics object
length(intersect(rownames(cisOut@binarized.cisTopics$Topic1),rownames(topic_annot_ix)))
length(intersect(rownames(cisOut@binarized.cisTopics$Topic2),rownames(topic_annot_ix)))

for(i in 1:nTopics){
  if(nrow(cisOut@binarized.cisTopics[[i]])==0)
    next
  
  topic_annot_ix[rownames(cisOut@binarized.cisTopics[[i]]),i] <- 1
  
}

colnames(topic_annot_ix) <- paste0("Topic",1:nTopics)

saveRDS(topic_annot_ix,"./processed_results/DR/cisTopics/scATAC_filt_cisTopics_annot_ix.rds")
