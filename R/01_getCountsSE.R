### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University

library(SummarizedExperiment)
library(GenomicRanges)
library(dplyr)
library(BuenColors)



setwd("<data_analysis_folder>")
  
# Frag files downloaded under analysis folder data subfolder, 1 per sample
fragFiles <- list.files("./data/fragFiles",full.names = TRUE)

# SampleNames
sampleNames <- splitAndFetch(basename(fragFiles),".fragments.tsv.gz",1)

# PeakRanges based on bed file of peaks
peakRanges <- makeGRangesFromDataFrame(read.table("./data/peakFiles/GMP_TET2_sc_ImmGen_combined_peaks_final_300.bed",sep="\t",header=FALSE),seqnames.field = "V1",start.field = "V2",end.field = "V3")

counts.L <- list()

# Get counts for each
for(i in 1:length(sampleNames)){
  cat("Sample: ",sampleNames[i],"\n")
  cat("Frag file: ",fragFiles[i],"\n")
 
  SE <-   getCountsFromFrags(fragFile = fragFiles[i],peaks = peakRanges)
  SE <- SE[,sort(colnames(SE))]
  
  # Add metadata
  
  SE$SampleID <- sampleNames[i]
  # Mouse ID
  SE$Sample <- stringr::str_extract_all(sampleNames[i],pattern = "[ABCE][012]",simplify = TRUE)[,1]
  # Sort
  SE$Sort <- ifelse(grepl("gmp",sampleNames[i],ignore.case = TRUE),"GMP","Lin_neg_Cd11b")
  
  counts.L[[i]] <-  SE
  gc()
 
  cat("\n\n\n")
}

# Merge
countsSE <- do.call('cbind',counts.L)
countsSE

# Add genotype info (WT/Tet2 KO) based on ID

countsSE$Condition <- ifelse(countsSE$Sample %in% c("E0","C0"),"WT",
                             ifelse(countsSE$Sample %in% c("E2","C1","A0","A2"),"Tet2KO","Tet2KO_Sox4OE"))

# Check
table(countsSE$Condition,countsSE$Sample)

# Save unfiltered SE object
saveRDS(countsSE,"./processed_results/countFiles/scATAC_SE.rds")