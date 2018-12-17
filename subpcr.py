fw=('ACAAACACTTTGCTTCCTGCGCAATTCC','ACAAGGGTTTTGCTTCCTGCGCAATTCC','ACATCTCTTTTGCTTCCTGCGCAATTCC','ACCCTGACTTTGCTTCCTGCGCAATTCC','ACGGGAACTTTGCTTCCTGCGCAATTCC','AGACAGTATTTGCTTCCTGCGCAATTCC','AGCTTCCATTTGCTTCCTGCGCAATTCC','AGGTTGGGTTTGCTTCCTGCGCAATTCC','ATGAGGCTTTTGCTTCCTGCGCAATTCC','ATTGTACGTTTGCTTCCTGCGCAATTCC','CAGATGAATTTGCTTCCTGCGCAATTCC','CCTTAGATTTTGCTTCCTGCGCAATTCC') 
rv=('CTTCGCCGCGCTGGACGTATGCGACT','CTTCGCCGCGCTGGACGTTGGGATAT','CTTCGCCGCGCTGGACGCGAGAACAT','CTTCGCCGCGCTGGACGCTCGATTTG','CTTCGCCGCGCTGGACGGGTTGCATG','CTTCGCCGCGCTGGACGACAGGAATG','CTTCGCCGCGCTGGACGATTAGGTGG','CTTCGCCGCGCTGGACGACGAACGGG')
a7cri="GCCTGCTCATGGGCATGGGC"
a7lcri="GCCTGCTCATGGGCTCATCC"
f = open("/Users/lhac/Downloads/oss-browser-darwin-x64/subPCR.txt")
samplerow=('A','B','C','D','E','F','G','H')
flag=False
for line in f:

    seq=line.strip('\n')
    flagfind=0
    m=0
    n=0
    for rvp in rv:
        if(seq.find(rvp)>=0):
            head=">"+str(samplerow[n])
            flagfind=flagfind+1
            end=seq.find(rvp)+26
        n=n+1
    for fwp in fw:
        m=m+1
        if flagfind == 1:
            if(seq.find(fwp)>=0):
                flagfind=flagfind+1
                head=head+str(m)
                start=seq.find(fwp)
    if flagfind == 2 and start < end:
        #print start,end
        seq=seq[start:end]
        if len(seq)>50:
            print seq
                #blast(seq,a7lcri)
        
        
        


f.close()
