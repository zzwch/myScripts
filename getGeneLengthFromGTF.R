## Transcript Length
library(GenomicRanges)
library(rtracklayer)

## reading your mm10 gtf file and read only exons
GTFfile <- gzfile(description = "E:/gencode.vM13.annotation.gtf.gz", open = "r") ## your mm10 gtf file
GTF <- import.gff(GTFfile, format="gtf", genome="mm10", feature.type="exon") 
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

## calculation of length
calc_length <- function(x) {
  sum(elementMetadata(x)$widths)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length))
output <- data.frame(t(output))
output$id = rownames(output)
output$id = gsub(pattern = "(.*)[.](.*)", replacement = "\\1", x = output$id)
output = output[,c(2,1)]
colnames(output)[2] = "Transcript.Length"
head(output)

write.csv(output, file = "gencode.vM13.geneLenght.csv")
