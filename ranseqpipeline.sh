## define path to your seq-data(ptf)
## define path to your ref(ref)
ptf=/data1/rawdata/a7l/raw/
ref=/home/lhac/ref/V439


while read samplename;
 do
 fastp -i ${ptf}${samplename}.R1.clean.fastq.gz -I ${ptf}${samplename}.R2.clean.fastq.gz -o ${ptf}${samplename}_1.filt.fq.gz -O ${ptf}${samplename}_2.filt.fq.gz -h ${ptf}report/${samplename}.html
  done  < sample
