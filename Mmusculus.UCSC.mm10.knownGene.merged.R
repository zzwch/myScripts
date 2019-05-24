## merge trascript bed to be gene-centric bed. see https://www.biostars.org/p/70927/
## you may also try this: GTFtools: a Python package for analyzing various modes of gene models 
## use -m Option to Merge all exons of all isoforms of the same gene.
## https://www.biorxiv.org/content/biorxiv/early/2018/02/11/263517.full.pdf 
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(org.Mm.eg.db)
library(clusterProfiler)
exonsByGene <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,'gene')
mergedGene <- reduce(exonsByGene) 
mergedGene <- 
  ## Remove a few hundred 'tricky cases' like:
  ##  * trans-spliced genes (?? are there any in hg19.knownGene ??)
  ##  * genes on multiple chromosomes (eg: alt. haplotypes "chr_ctg9_hap1" )
  ##    Arguably(?) these should not be considered the 'same gene'.
  mergedGene[1==elementNROWS(runValue(strand(mergedGene))) &
               1==elementNROWS(runValue(seqnames(mergedGene)))]  

idMap <- bitr(names(mergedGene), 'ENTREZID', "SYMBOL", "org.Mm.eg.db")
geneid <- names(mergedGene)
genename <- idMap$SYMBOL[match(names(mergedGene), idMap$ENTREZID)]
anyDuplicated(genename)
#!(duplicated(idMap$ENTREZID) | duplicated(idMap$SYMBOL))
anyNA(genename)
which(is.na(genename))
genename[which(is.na(genename))] <- geneid[which(is.na(genename))]

names(mergedGene) <- genename


export(mergedGene,'Hsapiens.UCSC.hg19.knownGene.merged.bed')
