# Using samtools and awk to convert bam to fastq
## Final version
`samtools view T0-ZCC-2_TKR180800591_HT3GYCCXY_L4.mapped.sorted.bam | awk 'BEGIN {FS="\t"} !(and($2,256) || and($2,2048)) { r1=substr($1, 1, 16); r1q=r1; gsub("A","F",r1q);gsub("T","F",r1q);gsub("G","F",r1q);gsub("C","F",r1q);gsub("N","#",r1q); sn=substr($1,18); r2=$10; r2q=$11; print sn,"\n" r1 "\n+\n" r1q >> "T0-ZCC-2_R1.fq"; print sn "\n" r2 "\n+\n" r2q >> "T0-ZCC-2_R2.fq"; }'`

## for test (you may start from this first)
`samtools view T0-ZCC-2_TKR180800591_HT3GYCCXY_L4.mapped.sorted.bam | head | awk 'BEGIN {FS="\t"} !(and($2,256) || and($2,2048)) { r1=substr($1, 1, 16); r1q=r1; gsub("A","F",r1q);gsub("T","F",r1q);gsub("G","F",r1q);gsub("C","F",r1q);gsub("N","#",r1q); sn=substr($1,18); r2=$10; r2q=$11; print sn,"\n" r1 "\n+\n" r1q >> "T0-ZCC-2_R1.fq"; print sn "\n" r2 "\n+\n" r2q >> "T0-ZCC-2_R2.fq"; }'` 

## I created this commond after learning AWK
### awk links
https://www.runoob.com/linux/linux-comm-awk.html   
https://www.runoob.com/w3cnote/awk-built-in-functions.html   
### Other links
https://github.com/samtools/samtools/issues/688   
https://onestopdataanalysis.com/convert-bam-to-fasta/   
https://baezortega.github.io/2018/09/22/command-line-manipulation-sequence-files/   

## The following commands also inspired me
`samtools view file.bam | awk 'BEGIN {FS="\t"} {print "@" $1 "\n" $10 "\n+\n" $11}' > file.fq`   
`samtools fastq -F 0x900 in.bam > all_reads.fq`




