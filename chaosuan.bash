while read samplename;
 do
 fastp -i /data1/rawdata/a7l/raw/${samplename}.R1.clean.fastq -I /data1/rawdata/a7l/raw/${samplename}.R2.clean.fastq -o /data1/rawdata/a7l/raw/${samplename}_1.filt.fq -O /data1/rawdata/a7l/raw/${samplename}_2.filt.fq -h /data1/rawdata/a7l/report/${samplename}.html
  done  < sample

while read samplename;
do
/home/lhac/hisat2-2.1.0/hisat2 -p 60 -t --dta -x /home/lhac/ref/V439/Zea_mays.AGPv4.hisat2 -1 /data1/rawdata/a7/raw/${samplename}_1.filt.fq -2 /data1/rawdata/a7/raw/${samplename}_2.filt.fq -S /data1/rawdata/a7/aligned/${samplename}.sam
done  < sample

while read samplename;
do
bowtie2 -p 30 -x /home/lhac/ref/V439/Zea_mays.AGPv4.bowtie2 -1 /data1/rawdata/chip/raw/${samplename}_1.filt.fq -2 /data1/rawdata/chip/raw/${samplename}_2.filt.fq -S /data1/rawdata/chip/align/${samplename}.sam
  done  < sample


B519R-ChIPedA_L2_P501706_1.filt.fq


bowtie2 -p 30 -x /home/lhac/ref/V439/Zea_mays.AGPv4.bowtie2 -1 /data1/rawdata/chip/raw/B519R-ChIPedA_L2_P501706_1.filt.fq -2 /data1/rawdata/chip/raw/B519R-ChIPedA_L2_P501706_2.filt.fq -S /data1/rawdata/chip/raw/B519R-ChIPedA_L2_P501706.sam
bowtie2 -p 30 -x /home/lhac/ref/V439/Zea_mays.AGPv4.bowtie2 -1 /data1/rawdata/chip/raw/B519R-ChIPedA_L2_P501706_1.filt.fq -2 /data1/rawdata/chip/raw/B519R-ChIPedA_L2_P501706_2.filt.fq -S /data1/rawdata/chip/align/B519R-ChIPedA_L2_P501706.sam


macs14 -t "/data1/rawdata/chip/align/B519R-ChIPedA_L2_P501706.sam" --name "519R"  -g 2.4e10


/annotatePeaks.pl "/home/lhac/MACS-1.4.2/bin/519R_summits.bed" /home/lhac/ref/V439/Zea_mays.AGPv4.dna.toplevel.fa -gtf /home/lhac/ref/V439/Zea_mays.AGPv4.39.chr.gtf >/data1/rawdata/chip/515r
./annotatePeaks.pl "/home/lhac/MACS-1.4.2/bin/519R1_summits.bed" /home/lhac/ref/V439/Zea_mays.AGPv4.dna.toplevel.fa -gtf /home/lhac/ref/V439/Zea_mays.AGPv4.39.chr.gtf >/data1/rawdata/chip/515r1


##chip peak file need wig or bed file to visualization



while read samplename;
do
macs14 -t /data1/rawdata/chip/align/${samplename}.sam --format SAM --name /data1/rawdata/chip/peak/${samplename} -g 2.4e10  -B -S -p 1e6
  done  < sample
  while read samplename;
  macs2 callpeak -t ${samplename}.sam -n ${samplename} -g 2.4e10  -B -q 0.01
    done  < sample
	
	
	
	zsub -q blade1 -e "/public/home/liuhao/CHIP/A7/IGVTools/igvtools toTDF "/public/home/liuhao/CHIP/A7l/7-1-PE_control_lambda.bdg" "/public/home/liuhao/CHIP/A7l/7-1-PE_control_lambda.tdf" "/public/home/liuhao/CHIP/A7/IGVTools/v4.genome""

	
	while read samplename;  do  fastp -c -i /data1/rawdata/a7/raw/${samplename}.R1.clean.fastq -I /data1/rawdata/a7/raw/${samplename}.R2.clean.fastq -o /data1/rawdata/a7/raw/${samplename}_1.filt.fq -O /data1/rawdata/a7/raw/${samplename}_2.filt.fq -h /data1/rawdata/a7/report/${samplename}.html;   done  < sample
while read samplename;  do  fastp -c -i /data1/rawdata/a7/raw/${samplename}.R1.clean.fastq.1 -I /data1/rawdata/a7/raw/${samplename}.R2.clean.fastq.1 -o /data1/rawdata/a7/raw/${samplename}_1.filt.fq -O /data1/rawdata/a7/raw/${samplename}_2.filt.fq -h /data1/rawdata/a7/report/${samplename}.html;   done  < sample
while read samplename;  do  fastp -i /data1/rawdata/chip/raw/${samplename}.R1.clean.fastq -I /data1/rawdata/chip/raw/${samplename}.R2.clean.fastq -o /data1/rawdata/chip/raw/${samplename}_1.filt.fq -O /data1/rawdata/chip/raw/${samplename}_2.filt.fq -h /data1/rawdata/chip/report/${samplename}.html;   done  < sample
while read samplename; bowtie2 -p 30 -x /home/lhac/ref/V439/Zea_mays.AGPv4.bowtie2 -1 /data1/rawdata/chip/raw/${samplename}_1.filt.fq -2 /data1/rawdata/chip/raw/${samplename}_2.filt.fq -S /data1/rawdata/chip/align/${samplename}.sam;   done  < sample
while read samplename; bowtie2 -p 30 -x /home/lhac/ref/V439/Zea_mays.AGPv4.bowtie2 -1 /data1/rawdata/chip/raw/${samplename}_1.filt.fq -2 /data1/rawdata/chip/raw/${samplename}_2.filt.fq -S /data1/rawdata/chip/align/${samplename}.sam; done  < sample
while read samplename; do /home/lhac/hisat2-2.1.0/hisat2 -p 60 -t --dta -x /home/lhac/ref/V439/Zea_mays.AGPv4.hisat2 -1 /data1/rawdata/a7/raw/${samplename}_1.filt.fq -2 /data1/rawdata/a7/raw/${samplename}_2.filt.fq -S /data1/rawdata/a7/aligned/${samplename}.sam && samtools  sort -@60 /data1/rawdata/a7/aligned/${samplename}.sam -o /data1/rawdata/a7/aligned/${samplename}_sorted.bam && /home/lhac/stringtie -e -B -p 60 -G /public/home/liuhao/ref/Zea_mays.AGPv4.39.chr.gtf -o /data1/rawdata/a7/aligned/ballgown/${samplename}/${samplename}_sorted.gtf -l ${samplename} /data1/rawdata/a7/aligned/${samplename}_sorted.bam;  done  < sample
while read samplename; do bowtie2 -p 30 -x /home/lhac/ref/V439/Zea_mays.AGPv4.bowtie2 -1 /data1/rawdata/chip/raw/${samplename}_1.filt.fq -2 /data1/rawdata/chip/raw/${samplename}_2.filt.fq -S /data1/rawdata/chip/align/${samplename}.sam;   done  < sample
while read samplename; do macs14 -t /data1/rawdata/chip/align/${samplename}.sam --format SAM --name /data1/rawdata/chip/peak/${samplename} -g 2.4e10  -B -S -p 1e6;   done  < sample


annotatePeaks.pl "/public/home/liuhao/CHIP/A7/519R-merge-a-b-c-narrow.bed" /public/home/liuhao/ref/V439/Zea_mays.AGPv4lic/home/liuhao/ref/V439/Zea_mays.AGPv4.39.chr.gtf>"/public/home/liuhao/CHIP/A7/519R-merge-a-b-c-narrow.anno"