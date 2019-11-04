#!/public/home/liuhao/tool/anaconda2/envs/python34/bin/python3.4
import sys
from Bio import pairwise2


def fedit(target, query, header):
	flagedit = False
	aln = pairwise2.align.globalms(query, target, 2, -1, -10, -0.5)
	alntarget = aln[0][0]
	alnquery = aln[0][1]
	if alntarget.find("-") > 47:
		editsite = alntarget.find("-") - 44
		gaps = alntarget.count("-")
		header = header+"\t+"+str(editsite)+":-"+str(gaps)
		flagedit = True
	if alnquery.find("-") > 47:
		editsite = alnquery.find("-") - 44
		gaps = alnquery.count("-")
		header = header+"\t+"+str(editsite)+":+"+str(gaps)
		flagedit = True
	if not flagedit:
		header = header+"\tNE"
		print(header)
	else:
		print(header)
		print(alnquery)
		print(alntarget)


def DNA_com(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()


def DNA_rv(sequence):
    sequence = sequence.upper()
    return sequence[::-1]


if (len(sys.argv) != 3):
	print("need two argument: filename and threshold")
	sys.exit(0)


infile = sys.argv[1]

lenthres = sys.argv[2]
fw = ('ACAAACACTTTGCTTCCTGCGCAATTCC', 'ACAAGGGTTTTGCTTCCTGCGCAATTCC', 'ACATCTCTTTTGCTTCCTGCGCAATTCC', 'ACCCTGACTTTGCTTCCTGCGCAATTCC', 'ACGGGAACTTTGCTTCCTGCGCAATTCC', 'AGACAGTATTTGCTTCCTGCGCAATTCC',
      'AGCTTCCATTTGCTTCCTGCGCAATTCC', 'AGGTTGGGTTTGCTTCCTGCGCAATTCC', 'ATGAGGCTTTTGCTTCCTGCGCAATTCC', 'ATTGTACGTTTGCTTCCTGCGCAATTCC', 'CAGATGAATTTGCTTCCTGCGCAATTCC', 'CCTTAGATTTTGCTTCCTGCGCAATTCC')
fw1 = ("AATGCACATTTGCTTCCTGCGCAATTCC", "CCCGATAATTTGCTTCCTGCGCAATTCC", "AGGAATCCTTTGCTTCCTGCGCAATTCC", "CGATTGATTTTGCTTCCTGCGCAATTCC", "GAGTAGACTTTGCTTCCTGCGCAATTCC", "AACGTACCTTTGCTTCCTGCGCAATTCC",
       "GGGATGGTTTTGCTTCCTGCGCAATTCC", "ACTGCCGATTTGCTTCCTGCGCAATTCC", "TCATAGACTTTGCTTCCTGCGCAATTCC", "CAGGATACTTTGCTTCCTGCGCAATTCC", "TCGCGCTTTTTGCTTCCTGCGCAATTCC", "AAAGTTGCTTTGCTTCCTGCGCAATTCC")
rv = ('CTTCGCCGCGCTGGACGTATGCGACT', 'CTTCGCCGCGCTGGACGTTGGGATAT', 'CTTCGCCGCGCTGGACGCGAGAACAT', 'CTTCGCCGCGCTGGACGCTCGATTTG',
      'CTTCGCCGCGCTGGACGGGTTGCATG', 'CTTCGCCGCGCTGGACGACAGGAATG', 'CTTCGCCGCGCTGGACGATTAGGTGG', 'CTTCGCCGCGCTGGACGACGAACGGG')
a7cri = "TTGCTTCCTGCGCAATTCCTGCGGCAGCGCCTTTCTCTCTTTCAATGGACAGCGACTACGTTGCCAGCCTGCTCATGGGCATGGGCTCATCCGCCCCTGCGCTCGACTTCGCCGCGCTGGACG"
a7lcri = "TTGCTTCCTGCGCAATTCCTGCGGCAGCGCCTTTGTTTGTTTCAATGGACAGCGACTACATTGCCAGCCTGCTCATGGGCTCATCCGCCCCTGCACTCAGCTTCGCCGCGCTGGACG"
f = open(infile)
samplerow = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
flag = False

for line in f:
   
    seq=line
    if(seq.find("@"))>=0:
        continue
    if(seq.find("#"))>=0:
        continue
    if(seq.find("+"))>=0:
        continue
    # print count
    if len(seq) < int(lenthres):
        continue
    flagfind = 0
    m = 0
    n = 0
    end = ""
    start = ""
    for rvp in rv:
        if(seq.find(rvp) >= 0):
            head = ">"+str(samplerow[n])
            flagfind = flagfind+1
            end = seq.find(rvp)+17
        n = n+1
    for fwp in fw1:
        m = m+1
        if flagfind == 1:
            if(seq.find(fwp) >= 0):
                flagfind = flagfind+1
                head = head+str(m)+str(" 2F")
                start = seq.find(fwp)+9
        m=0
    for fwp in fw:
        m = m+1
        if flagfind == 1:
            if(seq.find(fwp) >= 0):
                flagfind = flagfind+1
                head = head+str(m)+str(" 1F")
                start = seq.find(fwp)+9
    if flagfind == 2:
        # print start,end
        seq = seq[start:end]
        # f = open("tmp.tmp",'w+')
        # f.write(head+"\t"+count+"\t"+seq
        # print len(seq)
        # print (head+"\t"+count+"\t"+seq)
        if seq.find("TCGACTTCGCCGCGC") >= 0:
            head = head+"\tABP7PCR\t"
            # print (head+"\t"+seq)
            fedit(a7cri, seq, head)
        if seq.find("TCAGCTTCGCCGCGC") >= 0:
            head = head+"\t7LikPCR\t"
            # print (head+"\t"+seq)
            fedit(a7lcri, seq, head)

    if flagfind == 0:
        m = 0
        n = 0
        end = ""
        start = ""
        seq = DNA_com(DNA_rv(seq))
        for rvp in rv:
            if(seq.find(rvp)>=0):
                head=">"+str(samplerow[n])
                flagfind=flagfind+1
                end=seq.find(rvp)+17
                n=n+1
        for fwp in fw1:
            m=m+1
            if flagfind == 1:
                if(seq.find(fwp)>=0):
                    flagfind=flagfind+1
                    head=head+str(m)+str(" 2F")
                    start=seq.find(fwp)+9
        m=0
        for fwp in fw:
            m=m+1
            if flagfind == 1:
                if(seq.find(fwp)>=0):
                    flagfind=flagfind+1
                    head=head+str(m)+str(" 1F")
                    start=seq.find(fwp)+9
        if flagfind == 2:
            seq=seq[start:end]
            if seq.find("TCGACTTCGCCGCGC")>=0:
                head=head+"\tABP7PCR\t"
                fedit(a7cri,seq,head)
            if seq.find("TCAGCTTCGCCGCGC")>=0:
                head=head+"\t7LikPCR\t"
                fedit(a7lcri,seq,head)

f.close()

