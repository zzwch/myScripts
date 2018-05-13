## An in-house Command Line Interface to process tag-based scRNA-Seq data.

### Installation
`pip install --editable .`
### Prerequisite
1. some common genome analysis tools, such as hisat2, samtools, htseq-count, bamtools, bam2fastx, R and so on.
   for more information see [tools] section in configs/mm10.config file.
2. genome and trancscriptome index files build by hisat2-build
### Usage
`smartliu -c mm10`
