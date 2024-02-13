# GMP_analyses_code

Analysis code related to work presented in [Schiroli et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.09.04.556230v1)

Code repository:


01_getCountsSE.R - Code for deriving single cell reads-in-peaks (RIPs) counts matrix from individual samples' fragments files for single cell ATAC-seq (scATAC-seq) data

02_filtCountsSEQC.R - Code for QC-filtering of single cell data   

03_runCisTopics.R  - Code for running cisTopics on scATAC-seq data

04_clustCells.R  - Code for dimensionality reduction and clustering application to the scATAC-seq data

05_runchromVAR.R  - Code for running chromVAR TF motif accessibility analysis using scTAC-seq data

06_DESeq.R - Code for running differential peak accessibility analysis between Tet2KO and WT GMP cells using DESeq2
