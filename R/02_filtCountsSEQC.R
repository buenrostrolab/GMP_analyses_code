### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University


library(SummarizedExperiment)
library(ggplot2)
library(BuenColors)


setwd("<data_analysis_folder>")


# Load unfiltered SE (see 01_getCountsSE.R)

countsSE <- readRDS("./processed_results/countFiles/scATAC_SE.rds")
peakRanges <- granges(countsSE)

# Plot QC
cellMeta <- as.data.frame(colData(countsSE))
cellMeta$density <- get_density(x = log10(cellMeta$uniqueNuclearFrags),y=cellMeta$FRIP,n = 500)

gQC <- cellMeta %>% ggplot(aes(x=log10(uniqueNuclearFrags),y=FRIP,color=density)) + 
  geom_point(size=0.1) + theme_bw() + scale_color_gradientn(colours = jdb_palette("flame_light")) + labs(x="log10 UNF")

# Apply QC filters
FRIP.cut <- 0.4
UNF.cut <- 1500
dupProp.cut <- 0.15

cellsToKeep <- cellMeta$FRIP >= FRIP.cut & cellMeta$uniqueNuclearFrags >= UNF.cut & cellMeta$duplicateProportion >= dupProp.cut
table(cellsToKeep)

percentPass <- table(cellsToKeep) / length(cellsToKeep)
names(percentPass) <- c("Fail","Pass")

gQC <- gQC + geom_hline(yintercept = FRIP.cut,color="red",linetype="dashed",size=0.8) +
  geom_vline(xintercept = log10(UNF.cut),color="red",linetype="dashed",size=0.8) + 
  ggtitle(label = paste0("No. cells: ",ncol(countsSE),"\nPass %: ",round(percentPass[2]*100,3)))

table(cellMeta$SampleID,cellsToKeep) / as.numeric(table(cellMeta$SampleID) )

countsSE.filt <- countsSE[,cellsToKeep]

# Save filtered
saveRDS(countsSE.filt,"./processed_results/countFiles/scATAC_SE_filt.rds")