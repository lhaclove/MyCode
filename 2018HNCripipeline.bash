while read samplename;
do
fastp -i ${samplename}_1.fq -I ${samplename}_2.fq -c -o /data1/${samplename}.1.filt.fq -O /data1/${samplename}.2.filt.fq
./flash /data1/${samplename}.1.filt.fq /data1/${samplename}.2.filt.fq -o ${samplename} -d /data1/
cat /data1/${samplename}.extendedFrags.fastq | grep GGAATTGCGCAGGAAGCAAA | grep CGTCCAGCGCGGCGAAG > /data1/${samplename}.rv
cat /data1/${samplename}.extendedFrags.fastq | grep TTTGCTTCCTGCGCAATTCC | grep CTTCGCCGCGCTGGACG > /data1/${samplename}.fw
python makeRC.py /data1/${samplename}.rv > /data1/${samplename}.rv.rc
cat /data1/${samplename}.fw  /data1/${samplename}.rv.rc > /data1/${samplename}.all.fw
sort /data1/${samplename}.all.fw |uniq -c > /data1/${samplename}.all.fw.sort
python findprimerv5.py /data1/${samplename}.all.fw.sort  200 1 >/data1/${samplename}.200.1F.xls
python findprimerv5.py /data1/${samplename}.all.fw.sort  200 2 >/data1/${samplename}.200.2F.xls
sort /data1/${samplename}.200.1F.xls >/data1/${samplename}.200.1F.xls
sort /data1/${samplename}.200.2F.xls >/data1/${samplename}.200.2F.xls

  done  < sample
  