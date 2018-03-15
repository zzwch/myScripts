## ----set workspace -----------
setwd(dir = "/data1/users/lizc07/HumanData/bzj180314/scenic/")
GENIE3_nCores = 20
SCENIC_nCores = 20
SCENIC_steps = c("1.2", "2", "3.1", "3.2", "4")
## ----setup, echo=FALSE, message=FALSE, warning=FALSE---------------------
# Suppress loading messages when building the HTML
suppressPackageStartupMessages({
    library(Biobase)
    library(data.table)
    library(reshape2)
})

# Do not convert strings to factors (IMPORTANT! Specially if reading-in GENIE3 text output)
# It is also important NOT to convert strings to factors (e.g. when reading from text files)
options(stringsAsFactors=FALSE)

## ----loadExprMat---------------------------------------------------------
load("humanFL.scenic.Rdata") # exprMat and annot
## ----setwd, results='hide', warning=FALSE, eval=FALSE--------------------
dir.create("data")
dir.create("int")
dir.create("output")

## ----cellInfo, fig.height=4, fig.width=4---------------------------------
# Color to assign to the variables (same format as for NMF::aheatmap)
pal_dolphin <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", "#8B0101", "#007502", "#FE0000", "#FFFF01", "#FF99CB", "#4A95FB", "#61FE69", "#9A7A01", "#017F8B", "#05FDFF", "grey")


colVars <- lapply(cellInfo, function(x) {tmp = sort(unique(x)); setNames(pal_dolphin[1:length(tmp)], tmp)})
## ----chooseOrg-----------------------------------------------------------
org <- "hg19"

## ----LoadDbs-------------------------------------------------------------
if(org=="hg19")
{
  library(RcisTarget.hg19.motifDatabases.20k)
  
  # Get genes in databases:
  data(hg19_500bpUpstream_motifRanking) # or 10kbp, they should have the same genes
  genesInDatabase <- hg19_500bpUpstream_motifRanking@rankings$rn
  
  # Get TFS in databases:
  data(hg19_direct_motifAnnotation)
  allTFs <- hg19_direct_motifAnnotation$allTFs
}

if(org=="mm9")
{
  library(RcisTarget.mm9.motifDatabases.20k)
  
  # Get genes in databases:
  data(mm9_500bpUpstream_motifRanking) # or 10kbp, they should have the same genes
  genesInDatabase <- mm9_500bpUpstream_motifRanking@rankings$rn
  
  # Get TFS in databases:
  data(mm9_direct_motifAnnotation)
  allTFs <- mm9_direct_motifAnnotation$allTFs
}

## ----filterGenesInDb-----------------------------------------------------
genesLeft_minCells <- rownames(exprMat)[rowSums(exprMat > 0) > 3]
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)

## ----saveFilteredExprMat-------------------------------------------------
exprMatrix_filtered <- exprMat[genesLeft_minCells_inDatabases, ]
save(exprMatrix_filtered, file="int/1.1_exprMatrix_filtered.RData")

## ----rmExprMat-----------------------------------------------------------
rm(exprMat)

## ----TFlist--------------------------------------------------------------
inputTFs <- allTFs[allTFs%in% rownames(exprMatrix_filtered)]
save(inputTFs, file="int/1.2_inputTFs.RData")

c(allTFs=length(allTFs), inputTFs=length(inputTFs))

## ----genie3, eval=FALSE--------------------------------------------------
# setwd("SCENIC_HumanEC")
load("int/1.1_exprMatrix_filtered.RData")
load("int/1.2_inputTFs.RData")
library(GENIE3)
library(doParallel)
library(doRNG)
# Run GENIE3
weightMatrix <- GENIE3(as.matrix(exprMatrix_filtered), regulators=inputTFs,nCores=GENIE3_nCores, verbose = T)
save(weightMatrix, file="int/1.3_GENIE3_weightMatrix.RData")

## ----correlationMat, eval=FALSE------------------------------------------
load("int/1.1_exprMatrix_filtered.RData")
corrMat <- cor(t(exprMatrix_filtered), method="spearman")
save(corrMat, file="int/1.4_corrMat.RData")

# To save storage space, you may save only the rows for TFs:
# corrMat <- corrMat[which(rownames(corrMat) %in% inputTFs),]

## ----runScenicWrapper, eval=FALSE----------------------------------------
library(SCENIC)
runSCENIC(exprMat=exprMatrix_filtered, org=org, cellInfo=cellInfo, colVars=colVars, nCores=SCENIC_nCores,  stepsToRun= SCENIC_steps)

## ----sessionInfo---------------------------------------------------------
date()
sessionInfo()

