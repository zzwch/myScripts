options(stringsAsFactors = F)
rm(list = ls())
setwd(dir = "/data1/users/lizongcheng/Tanglab_data/HSC_DNA_Methylation_Data/methylKit_pooledSamples/pipeline20161231/analysis/")
require(GenomicRanges)
regulatoryRegions <- readRDS(file = "ensemblRegulatoryRegions.rds")
tfbsRegions <- readRDS(file = "ensemblTFBSRegions.rds")

my.options <- list(
  dir.analysis = "/data1/users/lizongcheng/Tanglab_data/HSC_DNA_Methylation_Data/methylKit_pooledSamples/pipeline20161231/analysis/",
  dir.report = "/data1/users/lizongcheng/Tanglab_data/HSC_DNA_Methylation_Data/methylKit_pooledSamples/pipeline20161231/analysis/report_test/",
  dir.methylDB = "/data1/users/lizongcheng/Tanglab_data/HSC_DNA_Methylation_Data/methylKit_pooledSamples/pipeline20161231/pooled_methylDB/",
  sample.annotation = "sample_annotation.pooled.csv",
  sample.table.separator = ",",
  sample.ids = "Sample_ID",
  sample.files = "bismarkCovFile",
  sample.types = "Cell_type",
  sample.colors = c("DarkOrange","red","chartreuse4","magenta4","#66ccfd"),
  sample.levels = c("EC","T1","T2","E14","BM"),
  sample.compareList = list(ECvsT1 = c("EC","T1"),	T1vsT2 = c("T1","T2"),	T2vsE14=c("T2","E14"), E14vsBM =c("E14","BM"), ECvsT2 = c("EC","T2")),
  getBaseDB = T,
  getBaseDB.test = F, # read 1000 line for pipeline test
  readBaseDB = F,
  baseDB.options = list(read.dir = "/data1/users/lizongcheng/Tanglab_data/HSC_DNA_Methylation_Data/methylKit_pooledSamples/pipeline20161231/pooled_methylDB/methylDB 2016-12-31 qeE/",read.suffix = ".txt.bgz",
                        dbtype = "tabix",pipeline = "bismarkCoverage",assembly = "mm10",header = F,
						            skip = 0,sep = "\t",context = "CpG", resolution = "base", mincov = 0),
  
  getTileDB = T,
  getTileDB.options = list(win.size = 5000, step.size = 5000, cov.bases = 5, mc.cores = 1, save.db = T, suffix = "pool.baseCov0.tile5kCov5"),
  readTileDB = F,
  readTileDB.options = list(read.dir = "/data1/users/lizongcheng/Tanglab_data/HSC_DNA_Methylation_Data/methylKit_pooledSamples/pipeline20161231/pooled_methylDB/methylDB 2016-12-31 qeE/",read.suffix = "_pool.baseCov0.tile5kCov5.txt.bgz",
                            dbtype = "tabix",pipeline = "bismarkCoverage",assembly = "mm10",header = F,
						                skip = 0,sep = "\t",context = "CpG", resolution = "region", mincov = 0),
  
  getRegionDB = T,
  regions.list = list(regu = regulatoryRegions),
  regions.select = list(regu = names(regulatoryRegions)),
  regionDB.prefix = "pool.baseCov0.regionCov5.",
  getRegionDB.options = list(cov.bases = 5, strand.aware = F, chunk.size = 1e6, save.db = T),
  readRegionDB = F,
  readRegionDB.options = list(read.dir = "/data1/users/lizongcheng/Tanglab_data/HSC_DNA_Methylation_Data/methylKit_pooledSamples/pipeline20161231/pooled_methylDB/",read.suffix = ".txt.bgz",
                              dbtype = "tabix",pipeline = "bismarkCoverage",assembly = "mm10",header = F,
						                  skip = 0,sep = "\t",context = "CpG", resolution = "region", mincov = 0),
  dmrTile = T,
  dmrTile.rds = "pool.baseCov0.tile5kCov5.dmr.rds",
  
  dmrRegion = T,
  dmrRegion.rds = "pool.baseCov0.regionCov5.dmr.rds",
  
  dmr.cores = 6
)
source(file = "pipeline/functions.R")
source(file = "pipeline/run.analysis.R")
instance.options <- my.options
run.analysis(instance.options)
saveRDS(my.options, file = paste0(my.options$dir.report, "/my.options.rds"))

# initial <- list(getBaseDB = T, getBaseDB.test = T, readBaseDB = F,
#                      getTileDB = T, readTileDB = F,
#                      getRegionDB = T, readRegionDB = F)
# rerun <- list(getBaseDB = F, getBaseDB.test = F, readBaseDB = T,
#               getTileDB = F, readTileDB = T,
#               getRegionDB = F, readRegionDB = T)
# 
# instance.options <- my.options
# instance.options[names(initial)] <- initial
# run.analysis(instance.options)
# 
# instance.options <- my.options
# instance.options[names(rerun)] <- rerun
# run.analysis(instance.options)
