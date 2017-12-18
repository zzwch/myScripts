#!/bin/bash
if [[ $# -lt 1 ]]; then
  echo "usage: `basename $0` -s sample_dir1,sample_dir2,... -r resolution1,resolution2,... -p pool_dir"
  echo "This script will produce sum of *_iced.matrix for each chromosome and Pooled matrix data"
  exit
fi

while getopts "s:r:p:" arg
do
  case $arg in
    s)
      samples=$OPTARG
      ;;
    r)
      resolution=$OPTARG
      ;;
    p)
      pool_dir=$OPTARG
      ;;
    ?)
      echo "Unknow Argument"
      exit 1
      ;;
  esac
done

samples=${samples:?Please set -s sample_dirs!}
pool_dir=${pool_dir:?Please set -p pool_dir!}
resolution=${resolution:?Please set -r resolution!}

if [ -d $pool_dir ]; then
  echo "-p pool_dir exists, Please delete it first."
  exit 1
fi

# strsplit by ','
sample_list=`echo $samples | awk -F',' 'BEGIN {OFS="\t"} {NF=NF; print $0}'`
resolution_list=`echo $resolution | awk -F',' 'BEGIN {OFS="\t"} {NF=NF; print $0}'`
echo $sample_list
echo $resolution_list

# processing
echo 'Processing start...'
echo 'Summing counts for each chromosome...'
for i in $sample_list;do
  echo $i
  sample_i=`basename $i`
  for j in $resolution_list;do
    echo "  $j"
    awk 'NR==FNR {ary[$4]=$1} NR>FNR {if(ary[$1]==ary[$2]) sum[ary[$1]]+=$3} END {for(key in sum) print key"\t"sum[key]}' \
        $i/raw/$j/${sample_i}_${j}_abs.bed $i/iced/$j/${sample_i}_${j}_iced.matrix | sort -g > $i/iced/$j/${sample_i}_${j}_iced.intraChrSums
    awk 'NR==FNR {ary[$4]=$1} NR>FNR {if(ary[$1]!=ary[$2]) sum[ary[$1]"\t"ary[$2]]+=$3} END {for(key in sum) print key"\t"sum[key]}' \
        $i/raw/$j/${sample_i}_${j}_abs.bed $i/iced/$j/${sample_i}_${j}_iced.matrix | sort -k1g,1 -k2g,2 > $i/iced/$j/${sample_i}_${j}_iced.interChrSums
  done
done

# pool data
echo "Pooling data ..."
mkdir -p $pool_dir
sample_pool=`basename $pool_dir`
for j in $resolution_list;do 
  echo "$j"
  mkdir -p $pool_dir/raw/$j
  mkdir -p $pool_dir/iced/$j
  sample_ji=''
  for i in $sample_list;do
    sample_i=`basename $i`
    sample_ji="$sample_ji $i/raw/$j/${sample_i}_${j}.matrix"
  done
  echo "  sort and merge triplet sparse matrix"
  cat $sample_ji | awk '{ary[$1"\t"$2]+=$3} END {for(key in ary) print key"\t"ary[key]}' | sort -k1g,1 -k2g,2 > $pool_dir/raw/$j/${sample_pool}_${j}.matrix 
  cp $i/raw/$j/${sample_i}_${j}_abs.bed $pool_dir/raw/$j/${sample_pool}_${j}_abs.bed
  echo "  iced normalization"
  python /home/Public/BioSoft/HiC-Pro_2.9.0/scripts/ice \
       --results_filename $pool_dir/iced/$j/${sample_pool}_${j}_iced.matrix \
       --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 \
       --remove-all-zeros-loci --output-bias 1 --verbose 1 $pool_dir/raw/$j/${sample_pool}_${j}.matrix
  echo "  Summing counts for each chromosome..."
  awk 'NR==FNR {ary[$4]=$1} NR>FNR {if(ary[$1]==ary[$2]) sum[ary[$1]]+=$3} END {for(key in sum) print key"\t"sum[key]}' \
        $pool_dir/raw/$j/${sample_pool}_${j}_abs.bed $pool_dir/iced/$j/${sample_pool}_${j}_iced.matrix | sort -g > $pool_dir/iced/$j/${sample_pool}_${j}_iced.intraChrSums
  awk 'NR==FNR {ary[$4]=$1} NR>FNR {if(ary[$1]!=ary[$2]) sum[ary[$1]"\t"ary[$2]]+=$3} END {for(key in sum) print key"\t"sum[key]}' \
        $pool_dir/raw/$j/${sample_pool}_${j}_abs.bed $pool_dir/iced/$j/${sample_pool}_${j}_iced.matrix | sort -k1g,1 -k2g,2 > $pool_dir/iced/$j/${sample_pool}_${j}_iced.interChrSums
done
