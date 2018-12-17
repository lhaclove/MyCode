echo "Shell 传递参数实例！";
echo "执行的文件名：$0";
echo "第一个参数为：$1";
echo "第二个参数为：$2";
echo "第三个参数为：$3";
if [$0=="PE"]
fi



prefix=$PWD
tail=
report=

if [ ! -d "./report/" ];then
mkdir ./report
fi

if [ ! -d "./align/" ];then
mkdir ./align
fi

if [ ! -d "./ballgown/" ];then
mkdir ./ballgown
fi



while read samplename;
do
fastp -i $prefix/raw/${samplename}.R1.clean.fastq.gz -I $prefix/raw/${samplename}.R2.clean.fastq.gz -o $prefix/raw/${samplename}.1.filt.fq.gz -O $prefix/raw/${samplename}.2.filt.fq.gz -h $tprefix/raw/${samplename}.html &&
hisat2 -p 60 -t --dta -x $REF439HTI -1 $prefix/raw/${samplename}.1.filt.fq.gz -2 $prefix/raw/${samplename}.2.filt.fq.gz -S $prefix/align/${samplename}.sam &&
samtools  sort -@60 $prefix/aligned/${samplename}.sam -o $prefix/aligned/${samplename}_sorted.bam &&
stringtie -e -B -p 60 -G $REF439GTF -o $prefix/ballgown/${samplename}/${samplename}_sorted.gtf -l ${samplename} $prefix/aligned/${samplename}_sorted.bam;
  done  < sample