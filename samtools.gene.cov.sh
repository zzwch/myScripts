samtools view -b -o accepted_hits.bam -q 60 accepted_hits.sam
samtools sort -o accepted_hits.sorted.bam accepted_hits.bam
samtools index accepted_hits.sorted.bam
rm -f accepted_hits.sam accepted_hits.bam
samtools depth -r chr1:4776465-4776470 accepted_hits.sorted.bam > gene.depth
