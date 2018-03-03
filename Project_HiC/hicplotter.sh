#!/bin/bash
if [[ $# -lt 1 ]]; then
  echo "usage: `basename $0` -i in_dir -o out_prefix -s sample_name1,sample_name2,... -r resolution1,resolution2,... -c chr1,chr2,... "
  echo "This script will plot HiC-matrix usng HiCPlotter.py and *_iced_intraChrSums.matrix."
  echo "Produce plots into out_dir"
  exit
fi

while getopts "s:r:c:i:o:" arg
do
  case $arg in
    s)
      samples=$OPTARG
      ;;
    r)
      resolution=$OPTARG
      ;;
    c)
      chromosome=$OPTARG
      ;;
    i)
      in_dir=$OPTARG
      ;;
    o)
      out_prefix=$OPTARG
      ;;
    ?)
      echo "Unknow Argument"
      exit 1
      ;;
  esac
done

samples=${samples:?Please set -s sample_dirs!}
in_dir=${in_dir:?Please set -i in_dir!}
out_prefix=${out_prefix:?Please set -o out_prefix!}
resolution=${resolution:?Please set -r resolution!}
chromosome=${chromosome:-'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY'}
# strsplit by ','
sample_list=`echo $samples | awk -F',' 'BEGIN {OFS="\t"} {NF=NF; print $0}'`
samplename_list=`echo $samples | awk -F',' 'BEGIN {OFS=" "} {NF=NF; print $0}'`
resolution_list=`echo $resolution | awk -F',' 'BEGIN {OFS="\t"} {NF=NF; print $0}'`
chromosome_list=`echo $chromosome | awk -F',' 'BEGIN {OFS="\t"} {NF=NF; print $0}'`
echo $sample_list
echo $resolution_list
echo $chromosome_list
# processing
echo 'Processing start...'
mkdir -p `dirname $out_prefix`
#plot_prefix=`basename $out_prefix`
for j in $resolution_list;do 
  echo "  $j"
  samplefile_list=""
  standard_flag=0
  for i in $sample_list;do
    samplefile_list="$samplefile_list $in_dir/${i}_${j}_iced_intraChrNorm.matrix"
    if [[ $standard_flag -eq 0 ]];then
      bed_standard="$in_dir/${i}_${j}_abs.bed"
      #out_prefix="$i.intraChrNorm"
      standard_flag=1
      continue
    fi
  done
  for chr in $chromosome_list;do
    echo "    $chr"
    HiCPlotter.py -f $samplefile_list \
                  -n $samplename_list \
                  -bed $bed_standard \
                  -o $out_prefix -tri 1 -chr $chr -r $j \
                  -mm 2 -ptr 1 \
                  -c 1 -p 1
  done
done

