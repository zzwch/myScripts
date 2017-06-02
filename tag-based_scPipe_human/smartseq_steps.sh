#!/bin/bash
source ~/.bash_profile
pipelineDir="/home/lizongcheng/myScripts/smartseq_pipeline_human/"
workdir=$1
sample=$2
processn=$3
barcode_list=$4
count=`less $barcode_list | wc -l`

indir="$workdir/data_raw/"
outdir="$workdir/out_smartseq/"
pl_split="${pipelineDir}/scripts/0.1.smartseq_tang.paired2single.pl"
pl_trimtso="${pipelineDir}/scripts/0.2.trim_TSO_polyA.pl"
pl_trim="${pipelineDir}/scripts/0.3.QC_rm_primer-SE_truncatedReads.lzc.pl"
R1="${sample}*1.*f*q.gz"
R2="${sample}*2.*f*q.gz"
barcode_shift=0
prefix="${sample}_sc"

transcriptome_index="/data1/users/lizongcheng/genomes/Homo_sapiens/hg19_Tang/hg19_ERCC92_RGC.refGene"
genome="/data1/users/lizongcheng/genomes/Homo_sapiens/hg19_Tang/hg19_ERCC92_RGC"
refGene_GTF="/data1/users/lizongcheng/genomes/Homo_sapiens/hg19_Tang/hg19_ERCC92_RGC.refGene.gff"
py_umicount="${pipelineDir}/scripts/0.4.UMI_HTseq.py"

resdir="$workdir/out_results/"

f_split(){
	mkdir -p ${outdir}/${sample}/
	perl $pl_split ${indir}/$R1 ${indir}/$R2 $barcode_shift ${outdir}/${sample}/$prefix $barcode_list

}

f_tophat_htseq_copy(){
mkdir -p $resdir
for i in `seq 1 $count`; do 
  myqsub -r -n $prefix$i -c "
perl $pl_trimtso $outdir/${sample}/$prefix$i/$prefix$i.UMI.fq.gz $outdir/$sample/$prefix$i/$prefix$i.UMI.trim.fq.gz 0 \n
perl $pl_trim --indir $outdir/$sample/$prefix$i --outdir $outdir/$sample/$prefix$i --sample $prefix${i}.UMI.trim\n
tophat -p 8 --transcriptome-index $transcriptome_index --library-type fr-unstranded -o $outdir/$sample/$prefix$i/tophat_out $genome $outdir/$sample/$prefix$i/$prefix$i.UMI.trim.clean.fq.gz \n
samtools sort -n --output-fmt SAM -o $outdir/$sample/$prefix$i/tophat_out/accepted_hits.sort_name.sam $outdir/$sample/$prefix$i/tophat_out/accepted_hits.bam \n
htseq-count -s no -f sam -o $outdir/$sample/$prefix$i/tophat_out/accepted_hits.count.sam $outdir/$sample/$prefix$i/tophat_out/accepted_hits.sort_name.sam $refGene_GTF > $outdir/$sample/$prefix$i/tophat_out/accepted_hits.count.txt \n
python $py_umicount $outdir/$sample/$prefix$i/tophat_out/accepted_hits.count.sam $outdir/$sample/$prefix$i/tophat_out/accepted_hits.count.UMIcount.txt \n
cp $outdir/$sample/$prefix$i/tophat_out/align_summary.txt $resdir/$prefix$i.align_summary.txt \n
cp $outdir/$sample/$prefix$i/tophat_out/accepted_hits.count.UMIcount.txt $resdir/$prefix$i.UMIcount.txt \n
cp $outdir/$sample/$prefix$i/$prefix$i.UMI.trim.QC.log $resdir/" -d 1 -p 8 -m 10gb
  while [ `qstat | wc -l` -gt $processn ]
  do 
	sleep 60
  done
done
}


f_split;
f_tophat_htseq_copy;
