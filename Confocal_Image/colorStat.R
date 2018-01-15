# Functional code, which had been adapted to shiny server code.
# Save it here as a reminder.
options(stringsAsFactors = F)
rm(list = ls())
library(tiff)
tiff1 <- readTIFF(source = "20-1(20m)-l-3---Cy5(对比度后).tif", native = F)
tiff2 <- readTIFF(source = "20-1(20m)-l-3---FITC（对比度后）.tif", native = F)
tiff3 <- readTIFF(source = "20-1(20m)-l-3---TRITC（对比度后）.tif", native = F)
## save tiff per channel
cutoff <- 0.5
tiff <- tiffBinary("test", tiff1, tiff2, tiff3, cutoff, tiff.save = T)
dims <- dim(tiff1)

## select points
points.red <- triplet(tiff[,,1], binary.th = cutoff)
points.green <- triplet(tiff[,,2], binary.th = cutoff)
points.blue <- triplet(tiff[,,3], binary.th = cutoff)

tiff.resize <- EBImage::resize(EBImage::as.Image(tiff), dims[1]/3, dims[2]/3)@.Data #

op <- par(bg = "black")
plot(c(0, dims[2]), c(0, -dims[1]), type = "n", xlab = "", ylab = "")
rasterImage(tiff.resize, 0, -dims[1], dims[2], 0, interpolate = F)
par(op)

polygon <- locator(n = 512, type = 'l')
points(polygon, col = "red", type = "l", new = F)
points.all <- SDMTools::pnt.in.poly(
  pnts = points.red[, c(1, 2)],
  poly.pnts = polygon
)

###functions
triplet <- function(mat, binary.th = 0.5){
  ind <- which(mat > binary.th, arr.ind = T)
  return(data.frame(x = ind[,2], y = -ind[,1]))
}
tiffBinary <- function(file.prefix, red, green, blue, cutoff = 0.5, tiff.save = T){
  require(tiff)
  red[red < cutoff] <- 0
  green[green < cutoff] <- 0
  blue[blue < cutoff] <- 0
  rgb <- red
  rgb[,,2] <- green[,,2]
  rgb[,,3] <- blue[,,3]
  if(tiff.save){
    writeTIFF(what = red, where = paste0(file.prefix, ".red.tiff"), compression = "LZW")
    writeTIFF(what = green, where = paste0(file.prefix, ".green.tiff"), compression = "LZW")
    writeTIFF(what = blue, where = paste0(file.prefix, ".blue.tiff"), compression = "LZW")
    writeTIFF(what = rgb, where = paste0(file.prefix, ".rgb.tiff"), compression = "LZW")
  }
  return(rgb)
}
