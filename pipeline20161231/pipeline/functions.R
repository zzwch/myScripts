# useful functions
f.getBaseDB <- function(locations, sampleIDs, treatment, baseDB.options, dir.methylDB){
  
  setwd(dir.methylDB)# let the methylDB folders under dir.methylDB
  myBaseData <- with(baseDB.options, 
                     methRead(location = locations, sample.id = sampleIDs, assembly = assembly, dbtype = dbtype,
                         pipeline = pipeline,header = header, skip = skip, sep = sep, context = context, resolution = resolution,
                         treatment = treatment, dbdir = getwd(), mincov = mincov)
  )
  print("methylBaseDB files created!")
  return(myBaseData)
}
f.readBaseDB <- function(locations, sampleIDs, treatment, baseDB.options){
  myBaseData <- with(baseDB.options,
                     methRead(location = locations, sample.id = sampleIDs, assembly = assembly, dbtype = dbtype,
                         pipeline = pipeline,header = header, skip = skip, sep = sep, context = context, resolution = resolution,
                         treatment = treatment, dbdir = getwd(), mincov = mincov)
  )
  print("read methylBaseDB successfully!")
  return(myBaseData)
}



f.getTileDB <- function(myBaseData, getTileDB.options, dir.methylDB){
  require(methylKit)
  setwd(dir.methylDB)
  myTileData <- with(getTileDB.options,tileMethylCounts(myBaseData, win.size, step.size, cov.bases, mc.cores, save.db,suffix = suffix))
  print("tileDB files created!")
  return(myTileData)
}
f.readTileDB <- function(locations, sampleIDs, treatment,readTileDB.options){
  require(methylKit)
  myTileData <- with(readTileDB.options,
                     methRead(location = locations, sample.id = sampleIDs, assembly = assembly, dbtype = dbtype,
                         pipeline = pipeline,header = header, skip = skip, sep = sep, context = context, resolution = resolution,
                         treatment = treatment, dbdir = getwd(), mincov = mincov)
  )
  print("read methylTileDB successfully!")
  return(myTileData)
}

f.getRegionDB <- function(myBaseData, myRegions, regionDB.prefix, getRegionDB.options, dir.methylDB){
  require(methylKit)
  myRegionData <- list()
  for(i in names(myRegions)){
    for(j in names(myRegions[[i]])){
      tmpRegionData <- with(getRegionDB.options,regionCounts(myBaseData, regions = unique(myRegions[[i]][[j]]), cov.bases = cov.bases, strand.aware = strand.aware, chunk.size = chunk.size,save.db = F))
      if(getRegionDB.options$save.db){
        setwd(dir.methylDB)
        lapply(tmpRegionData, makeMethylDB, paste0(regionDB.prefix,j))
      }
      myRegionData[[j]] <- tmpRegionData
    }
  }
  print("regionDB files created!")
  return(myRegionData)
}
f.readRegionDB <- function(locations, sampleIDs, treatment, readRegionDB.options){
  require(methylKit)
  myRegionData <- with(readRegionDB.options,
                       methRead(location = locations, sample.id = sampleIDs, assembly = assembly, dbtype = dbtype,
                           pipeline = pipeline,header = header, skip = skip, sep = sep, context = context, resolution = resolution,
                           treatment = treatment, dbdir = getwd(), mincov = mincov)
  )
  print("read methylRegionDB successfully!")
  return(myRegionData)
}

f.getDMR <- function(sample.compareList, regionData, dmr.cores){
  myDiffData <- list()
  myDiffRes <- list()
  for(i in names(sample.compareList)){
    myDiffData[[i]] <- reorganize(regionData, sample.ids = sample.compareList[[i]], treatment = c(0,1),save.db = F)
    myDiffRes[[i]] <- calculateDiffMeth(na.omit(myDiffData[[i]]), mc.cores = dmr.cores, slim = F, save.db = F)
  }
  myDiffRes25p <- list()
  for(i in names(sample.compareList)){
    myDiffRes25p[[i]] <- getMethylDiff(myDiffRes[[i]], difference = 25, qvalue = 0.05, save.db = F)
  }
  print("dmr calculated!")
  return(list(diffData = myDiffData, diffRes = myDiffRes, diffRes25p = myDiffRes25p))
}