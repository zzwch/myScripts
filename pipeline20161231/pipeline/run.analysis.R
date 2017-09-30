# this script is designed for semi-automatic methylKit analysis
# Required R packages: methylKit, 
#
# NOTICE: currentyly, only support for MANUAL MODE 
#         that means you should specify everything in the setting section by yourself
run.analysis <- function(my.options){
  options(stringsAsFactors = F)
  require(methylKit)
  attach(my.options)
  # set analysis path as workspace
  if(!dir.exists(dir.analysis)){
    stop(paste0("dir.analysis ",dir.analysis," did not exist. please create it!"))
  }
  setwd(dir.analysis)
  print(paste0(toupper("current workdir: "),getwd()))
  # creating report directory
  if(dir.exists(dir.report)){
    detach(my.options)
    stop("dir.report directory existed!")
  }else{
    dir.create(dir.report)
    print("INFO: dir.report directory created!")
  }
  
  sampleInfo <- read.table(file = paste0(dir.analysis, sample.annotation),header = T, sep = sample.table.separator)
  sampleIDs <- as.list(sampleInfo[[sample.ids]])
  treatment <- as.numeric(factor(sampleInfo[[sample.types]])) - 1
  
  # first run
  if(getBaseDB){
    if(dir.exists(dir.methylDB)){
      print("INFO: methylDB directory exists, and continue to the next step.")
    }else{
      dir.create(dir.methylDB)
      print("INFO: methylDB directory created!")
    }
    locations <- as.list(sampleInfo[[sample.files]])
    if(getBaseDB.test){
      locations <- lapply(locations, function(x){
                                      gz1 <- gzfile(x)
                                      gz2 <- gzfile(paste0(x,".test.gz"), "w")
                                      writeLines(text = readLines(gz1, n = 1000), gz2)
                                      close(gz1);close(gz2)
                                      return(paste0(x,".test.gz"))
                                      })
    }
    myBaseData <- f.getBaseDB(locations, sampleIDs, treatment, baseDB.options, dir.methylDB)
    
  }else if(getTileDB || getRegionDB || readBaseDB){
    locations <- as.list(paste0(baseDB.options$read.dir, sampleInfo[[sample.ids]], baseDB.options$read.suffix))
    myBaseData <- f.readBaseDB(locations, sampleIDs, treatment, baseDB.options)
  }
  
  # get tile data
  if(getTileDB){
    myTileData <- f.getTileDB(myBaseData, getTileDB.options, dir.methylDB)
  }else if(readTileDB){
    locations <- as.list(paste0(readTileDB.options$read.dir, sampleInfo[[sample.ids]], readTileDB.options$read.suffix))
    myTileData <- f.readTileDB(locations, sampleIDs, treatment, readTileDB.options)
  }
  # get region data
  myRegions <- list()
  for(i in names(regions.list)){
    tmp.region <- regions.list[[i]][regions.select[[i]]]
    myRegions[[i]] <- tmp.region
  }
  if(getRegionDB){
    myRegionData <- f.getRegionDB(myBaseData, myRegions, regionDB.prefix, getRegionDB.options, dir.methylDB)
  }else if(readRegionDB){
    myRegionData <- list()
    for(i in names(myRegions)){
      for(j in names(myRegions[[i]])){
        locations <- as.list(paste0(readRegionDB.options$read.dir, regionDB.prefix, j, "/", sampleInfo[[sample.ids]], readRegionDB.options$read.suffix))
        myRegionData[[j]] <- f.readRegionDB(locations, sampleIDs, treatment, readRegionDB.options)
      }
    }
  }
  
  # dmr tile
  myTileUnite <- unite(object = myTileData, destrand = F, save.db = F, min.per.group = 0L)
  if(nrow(na.omit(myTileUnite)) > 0){
    saveRDS(f.getDMR(sample.compareList, myTileUnite, dmr.cores),file = paste0(dir.report, "tile.dmr.rds"), compress = F)
  }
  # dmr region 
  for(i in names(myRegionData)){
    print(i)
    myRegionUnite <- unite(object = myRegionData[[i]], destrand = F, save.db = F, min.per.group = 0L)
    if(nrow(na.omit(myRegionUnite)) > 0){
      print("dd")
      saveRDS(f.getDMR(sample.compareList, myRegionUnite, dmr.cores),file = paste0(dir.report, i,".dmr.rds"), compress = F)
    }
  }
  detach(my.options)
}

