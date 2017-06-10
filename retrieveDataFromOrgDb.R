# using select() to extract data from AnnotationDb object.
# Here we take org.Hs.eg.db for example.
# we will extract genes included in Cell cycle GO term(GO:0007049) and its child terms. 
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
--------------
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
[11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
[21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"     
--------------

keytypes(org.Hs.eg.db)
--------------
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
[11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
[21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"    
--------------

gogenes <- select(org.Hs.eg.db, keys = c("GO:0007049"), columns = c("SYMBOL"), keytype = "GOALL") 
sort(unique(gogenes$SYMBOL))
