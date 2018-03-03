#!/bin/bash
if [[ $# -lt 1 ]]; then
  echo "usage: `basename $0` -s sample_dir1,sample_dir2,... -r resolution1,resolution2,... -o out_dir"
  echo "This script will do depth normalization for each chromosome to chromosome depth of sample_dir1."
  echo "Produce *_iced_intraChrNorm.matrix into out_dir"
  exit
fi

while getopts "s:r:o:" arg
do
  case $arg in
    s)
      samples=$OPTARG
      ;;
    r)
      resolution=$OPTARG
      ;;
    o)
      out_dir=$OPTARG
      ;;
    ?)
      echo "Unknow Argument"
      exit 1
      ;;
  esac
done

samples=${samples:?Please set -s sample_dirs!}
out_dir=${out_dir:?Please set -p out_dir!}
resolution=${resolution:?Please set -r resolution!}

# strsplit by ','
sample_list=`echo $samples | awk -F',' 'BEGIN {OFS="\t"} {NF=NF; print $0}'`
resolution_list=`echo $resolution | awk -F',' 'BEGIN {OFS="\t"} {NF=NF; print $0}'`
echo $sample_list
echo $resolution_list

# processing
echo 'Processing start...'
mkdir -p $out_dir
for j in $resolution_list;do 
  echo "$j"
  standard_flag=0
  for i in $sample_list;do
    sample_i=`basename $i`
    matrix_ji="$i/iced/$j/${sample_i}_${j}_iced.matrix"
    chrSum_ji="$i/iced/$j/${sample_i}_${j}_iced.intraChrSums"
    bed_ji="$i/raw/$j/${sample_i}_${j}_abs.bed"
    if [[ $standard_flag -eq 0 ]];then
      matrix_standard=$matrix_ji
      chrSum_standard=$chrSum_ji
      bed_standard=$bed_ji
      cp $matrix_standard $out_dir/${sample_i}_${j}_iced_intraChrNorm.matrix
      cp $bed_standard $out_dir/${sample_i}_${j}_abs.bed
      standard_flag=1
      echo "  Depth Normalization by $matrix_ji"
      continue
    fi
    echo "  Normalize Depth of $matrix_ji"
    awk 'ARGIND==1 {bed_ary[$4]=$1} ARGIND==2 {chr_ary[$1]=$2} ARGIND==3 {chr_ji_ary[$1]=$2} ARGIND==4 {if(bed_ary[$1]==bed_ary[$2]) print $1"\t"$2"\t"chr_ary[bed_ary[$1]]*$3/chr_ji_ary[bed_ary[$1]]}' $bed_standard $chrSum_standard $chrSum_ji $matrix_ji > $out_dir/${sample_i}_${j}_iced_intraChrNorm.matrix
  done
done

