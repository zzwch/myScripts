pipelineDir="/home/lizongcheng/myScripts/smartseq_pipeline_human/"
# $1 - workspace directory
# $2 - file including a raw data file name per line
# $3 - qsub mission nums under processing
# $4 - barcode file
workdir=$1
files=$2

if [ $# -lt 3 ];then
  processn=40
else
  processn=$3
fi

if [ $# -lt 4 ];then
  barcode="$pipelineDir/96-8bp-barcode.txt"
else
  barcode=$4
fi
cat $files | while read sample
do 
  sh "${pipelineDir}/smartseq_steps.sh" $workdir $sample $processn $barcode
done
