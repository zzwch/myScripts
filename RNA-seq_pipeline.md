# RNA-seq data analysis
## course
https://bioinformatics-core-shared-training.github.io/RNAseq_September_2019/

## RNA-seq preprocess shell script
```
#!/user/bin bash
# hisat2 https://daehwankimlab.github.io/hisat2/manual/
# rnaseq https://nf-co.re/rnaseq
# featureCounts http://bioinf.wehi.edu.au/featureCounts/
# sam2bam https://www.biostars.org/p/377961/

samples=('GFP_HETER1' 'GFP_HETER2' 'GFP_HETER3' 'GFP_HETER4' 'GFP_U1' 'GFP_U2' 'GFP_U3' 'GFP_U4' 'KG_HETER1' 'KG_HETER2' 'KG_HETER3' 'KG_U1' 'KG_U2' 'KG_U3')
fastq_dir='/data1/users/lizc07/WuhanMouse/fastq_data/longRNA'

reference='/data1/Database/genome/fromTangLab/mm10/mm10_ERCC92_RGC.refGene'
#samples=('test1' 'test2')
#fastq_dir='/data1/users/lizc07/WuhanMouse/fastq_data/testRNA'
for i in "${samples[@]}"; 
do 
  #r1=`ls $fastq_dir/${i}*_R1_*.fastq.gz | xargs echo | sed 's/ /,/g'`
  #r2=`ls $fastq_dir/${i}*_R2_*.fastq.gz | xargs echo | sed 's/ /,/g'`
  
  echo "RUN $i"
  echo "step0: cat multiple fastq.gz into one"
  mkdir -p rawdata
  cat $fastq_dir/${i}*_R1*.fastq.gz > rawdata/${i}_R1.fastq.gz
  cat $fastq_dir/${i}*_R2*.fastq.gz > rawdata/${i}_R2.fastq.gz
  
  r1="rawdata/${i}_R1.fastq.gz"
  r2="rawdata/${i}_R2.fastq.gz"
  
  echo "step1: trim_galore"
  #mkdir -p fastqc
  #fastqc -t 2 -o ./fastqc rawdata/${i}_R1.fastq.gz rawdata/${i}_R2.fastq.gz
  trim_galore --illumina --paired --fastqc -o trim_galore/ $r1 $r2
  
  r1="trim_galore/${i}_R1_val_1.fq.gz"
  r2="trim_galore/${i}_R2_val_2.fq.gz"
  echo "step2: hisat2"
  hisat2 -x $reference -1 $r1 -2 $r2 | \
    tee >(samtools flagstat - > ${i}.flagstat) | \
    samtools sort -O BAM | \
    tee ${i}.bam | \
    samtools index - ${i}.bam.bai
done
```

## use `multiqc .` to summary those reports.

## use featureCounts program included in the SourceForge Subread package
http://bioinf.wehi.edu.au/featureCounts/   

Summarize multiple paired-end datasets:   
`featureCounts -p -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam`

exon-level counting https://groups.google.com/g/subread/c/aPM6_4B6bCA   
`featureCounts -p -t exon -s 0 -T 12 -f -O -a genes.gtf -o featurecounts.txt sample.starAligned.sortedByCoord.out.bam`
