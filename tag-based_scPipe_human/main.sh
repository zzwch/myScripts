pipelineDir="/home/lizongcheng/myScripts/smartseq_pipeline_human/"
# $1 - workspace directory
# $2 - file including a raw data file name per line
# $3 - skip code 
#      (trim, tophat_refgene, tophat_gencode, htseq_refgene, htseq_gencode;)
#      0  - 0b000000 no skip
#      1  - 0b000001 skip htseq_gencode
#      2  - 0b000010 skip htseq_refgene
#      4  - 0b000100 skip tophat_gencode
#      8  - 0b001000 skip tophat_refgene
#      16 - 0b010000 skip trim
#      32 - 0b100000 skip split
# $4 - qsub mission nums under processing
# $5 - barcode file
workdir=$1
files=$2

if [ $# -lt 3 ];then
  skipcode=$((0))
else
  skipcode=$3
fi

if [ $# -lt 4 ];then
  processn=40
else
  processn=$4
fi

if [ $# -lt 5 ];then
  barcode="$pipelineDir/96-8bp-barcode.txt"
else
  barcode=$5
fi



cat $files | while read sample
do 
  sh "${pipelineDir}/smartseq_steps.sh" $workdir $sample $skipcode $processn $barcode
done
