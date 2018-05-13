## An in-house Command Line Interface to process tag-based scRNA-Seq data.

### Installation
`pip install --editable .`
### Prerequisite
1. some common genome analysis tools, see [tools] section in configs/mm10.config file to find more information.
   Including but no limited to: `hisat2`, `samtools`, `htseq-count`, `bamtools`, `bam2fastx`, `R`.
   You may use `conda install your-tool-name` to install them on your Linux server.
2. genome and trancscriptome index files build by hisat2-build
### Usage
`smartliu -c mm10`
