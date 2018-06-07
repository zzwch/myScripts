## An in-house Command Line Interface to process tag-based scRNA-Seq data.

### Installation
`cd path-to-smartliu`   
`pip install --editable .`
### Prerequisite
1. some common genome analysis tools, Including but no limited to: `hisat2`, `samtools`, `htseq-count`, `bamtools`, `bam2fastx`, `R`,`multiqc`.   
see [tools] section in configs/mm10.config file to find more information.   
You may use `conda install your-tool-name` to install them on your Linux server.
2. genome and trancscriptome index files build by hisat2-build
### Usage    
use `smartliu --help` to see how to start
eg.   
`smartliu -c mm10 -i raw_data -o smart_mm10 -p n_sample`
