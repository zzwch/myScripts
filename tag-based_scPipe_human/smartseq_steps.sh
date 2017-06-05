#!/bin/bash
source ~/.bash_profile
pipelineDir="/home/lizongcheng/myScripts/smartseq_pipeline_human/"
workdir=$1
sample=$2
skipcode=$3
processn=$4
barcode_list=$5
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

gencodev19_index="/data1/users/lizongcheng/genomes/Homo_sapiens/hg19_Tang/hg19_ERCC92_RGC.gencodev19/gencode.v19.annotation"
gencodev19_GTF="/data1/users/lizongcheng/genomes/Homo_sapiens/hg19_Tang/hg19_ERCC92_RGC.gencodev19/gencode.v19.annotation.gff"

resdir="$workdir/out_results/"

f_split(){
	mkdir -p ${outdir}/${sample}/
	perl $pl_split ${indir}/$R1 ${indir}/$R2 $barcode_shift ${outdir}/${sample}/$prefix $barcode_list

}

f_tophat_htseq_both(){
mkdir -p $resdir/refGene/
mkdir -p $resdir/gencode/

for i in `seq 1 $count`; do 

qsub_trim="perl $pl_trimtso $outdir/${sample}/$prefix$i/$prefix$i.UMI.fq.gz $outdir/$sample/$prefix$i/$prefix$i.UMI.trim.fq.gz 0 \n
perl $pl_trim --indir $outdir/$sample/$prefix$i --outdir $outdir/$sample/$prefix$i --sample $prefix${i}.UMI.trim\n"
qsub_tophat_refgene="tophat -p 8 --transcriptome-index $transcriptome_index --library-type fr-unstranded -o $outdir/$sample/$prefix$i/tophat_out_refgene $genome $outdir/$sample/$prefix$i/$prefix$i.UMI.trim.clean.fq.gz \n"
qsub_htseq_refgene="samtools sort -n --output-fmt SAM -o $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.sort_name.sam $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.bam \n
htseq-count -s no -f sam -o $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.count.sam $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.sort_name.sam $refGene_GTF > $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.count.txt \n
python $py_umicount $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.count.sam $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.count.UMIcount.txt \n"
qsub_tophat_gencode="tophat -p 8 --transcriptome-index $gencodev19_index --library-type fr-unstranded -o $outdir/$sample/$prefix$i/tophat_out_gencode $genome $outdir/$sample/$prefix$i/$prefix$i.UMI.trim.clean.fq.gz \n"
qsub_htseq_gencode="samtools sort -n --output-fmt SAM -o $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.sort_name.sam $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.bam \n
htseq-count -s no -f sam -o $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.count.sam $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.sort_name.sam $gencodev19_GTF > $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.count.txt \n
python $py_umicount $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.count.sam $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.count.UMIcount.txt \n"
qsub_copy_refgene="cp $outdir/$sample/$prefix$i/tophat_out_refgene/align_summary.txt $resdir/refGene/$prefix$i.align_summary.txt \n
cp $outdir/$sample/$prefix$i/tophat_out_refgene/accepted_hits.count.UMIcount.txt $resdir/refGene/$prefix$i.UMIcount.txt \n"
qsub_copy_gencode="cp $outdir/$sample/$prefix$i/tophat_out_gencode/align_summary.txt $resdir/gencode/$prefix$i.align_summary.txt \n
cp $outdir/$sample/$prefix$i/tophat_out_gencode/accepted_hits.count.UMIcount.txt $resdir/gencode/$prefix$i.UMIcount.txt \n"

if [ $(($skipcode & 16)) -eq 16 ];then qsub_trim=''; fi
if [ $(($skipcode & 8)) -eq 8 ];then qsub_tophat_refgene=''; fi
if [ $(($skipcode & 4)) -eq 4 ];then qsub_tophat_gencode=''; fi
if [ $(($skipcode & 2)) -eq 2 ];then qsub_htseq_refgene=''; fi
if [ $(($skipcode & 1)) -eq 1 ];then qsub_htseq_gencode=''; fi

  myqsub -r -n $prefix$i -c "
$qsub_trim\n
$qsub_tophat_refgene\n
$qsub_tophat_gencode\n
$qsub_htseq_refgene\n
$qsub_htseq_gencode\n
$qsub_copy_refgene\n
$qsub_copy_gencode\n
cp $outdir/$sample/$prefix$i/$prefix$i.UMI.trim.QC.log $resdir/ \n
" -d 1 -p 8 -m 10gb
  while [ `qstat | wc -l` -gt $processn ]
  do 
	sleep 60
  done
done
}

if [ $(($skipcode & 32)) -eq 32 ];then 
  echo "skip split barcode!\n";
else
  f_split;
fi

f_tophat_htseq_both;
