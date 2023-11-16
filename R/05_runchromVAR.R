### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University


library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(motifmatchr)


setwd("<data_analysis_folder>")


cat("Loading SE ..\n")
countsSE <- readRDS("./processed_results/countFiles/scATAC_SE_filt.rds")

# Non-zero peaks
cat("Non-zero peaks \n")
table(rowSums(assay(countsSE))!=0)

peaksToKeep <- rowSums(assay(countsSE))!=0

cat("Keeping non-zero peaks .. \n")
countsSE <- countsSE[peaksToKeep,]

cat("Adding GC bias ..\n")
countsSE <- addGCBias(countsSE,genome=BSgenome.Mmusculus.UCSC.mm10)

# Motifs
cat("Getting motif matches ..\n")
motif_ix <- matchMotifs(pwms = mouse_pfms_v4,
                        subject = countsSE,
                        genome=BSgenome.Mmusculus.UCSC.mm10)

cat("Getting background peaks ..\n")
set.seed(123)
bg_peaks <- getBackgroundPeaks(countsSE,niterations=100)

# Run chromVAR
cat("Running chromVAR ..\n")
BiocParallel::register(BiocParallel::MulticoreParam(4, progressbar = TRUE))
motif_dev <- computeDeviations(object = countsSE,
                               background_peaks=bg_peaks,
                               annotations = motif_ix)

cat("Saving output ..\n")
saveRDS(motif_dev,"./processed_results/chromVAR/scATAC_filt_motif_dev.rds")
cat("Finished!\n")
