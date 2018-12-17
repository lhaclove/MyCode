
def score(a,b):#scoring function
    score=0
    lst=['AC','GT','CA','TG']
    if a==b:
        score +=2
    elif a+b in lst:
        score += -5
    else:
        score += -7
    return score
 
def blast(seq1,seq2):
    #Basic Local Alignment Search Tool
    l1 = len(seq1)
    l2 = len(seq2)
    GAP =-5     #-5 for any gap
    scores =[]
    point =[]
    
    for j in range(l2+1):
        if j == 0:
            line1=[0]
            line2=[0]
            for i in range(1,l1+1):
                line1.append(GAP*i)
                line2.append(2)
        else:
            line1=[]
            line2=[]
            line1.append(GAP*j)
            line2.append(3)
        scores.append(line1)
        point.append(line2)
    
    #fill the blank of scores and point
    for j in range(1,l2+1):
        letter2 = seq2[j-1]
        for i in range(1,l1+1):
            letter1 = seq1[i-1]
            diagonal_score = score(letter1, letter2) + scores[j-1][i-1]
            left_score = GAP + scores[j][i-1]
            up_score = GAP + scores[j-1][i]
            max_score = max(diagonal_score, left_score, up_score)
            scores[j].append(max_score)
            
            if scores[j][i] == diagonal_score:
                point[j].append(1)
            elif scores[j][i] == left_score:
                point[j].append(2)
            else:
                point[j].append(3)
                
    #trace back
    alignment1=''
    alignment2=''
    i = l2
    j = l1
    #print 'scores =',scores[i][j]
    while True:
        if point[i][j] == 0:
            break
        elif point[i][j] == 1:
            alignment1 += seq1[j-1]
            alignment2 += seq2[i-1]
            i -= 1
            j -= 1
        elif point[i][j] == 2:
            alignment1 += seq1[j-1]
            alignment2 += '-'
            j -= 1
        else:
            alignment1 += '-'
            alignment2 += seq2[i-1]
            i -= 1
            
    #reverse alignment
    alignment1 = alignment1[::-1]
    alignment2 = alignment2[::-1]
    print 'The best alignment:'
    print alignment1
    print alignment2


fw=('ACAAACACTTTGCTTCCTGCGCAATTCC','ACAAGGGTTTTGCTTCCTGCGCAATTCC','ACATCTCTTTTGCTTCCTGCGCAATTCC','ACCCTGACTTTGCTTCCTGCGCAATTCC','ACGGGAACTTTGCTTCCTGCGCAATTCC','AGACAGTATTTGCTTCCTGCGCAATTCC','AGCTTCCATTTGCTTCCTGCGCAATTCC','AGGTTGGGTTTGCTTCCTGCGCAATTCC','ATGAGGCTTTTGCTTCCTGCGCAATTCC','ATTGTACGTTTGCTTCCTGCGCAATTCC','CAGATGAATTTGCTTCCTGCGCAATTCC','CCTTAGATTTTGCTTCCTGCGCAATTCC') 
rv=('CTTCGCCGCGCTGGACGTATGCGACT','CTTCGCCGCGCTGGACGTTGGGATAT','CTTCGCCGCGCTGGACGCGAGAACAT','CTTCGCCGCGCTGGACGCTCGATTTG','CTTCGCCGCGCTGGACGGGTTGCATG','CTTCGCCGCGCTGGACGACAGGAATG','CTTCGCCGCGCTGGACGATTAGGTGG','CTTCGCCGCGCTGGACGACGAACGGG')
a7cri="GCCTGCTCATGGGCATGGGC"
a7lcri="GCCTGCTCATGGGCTCATCC"
f = open("/Users/lhac/Downloads/oss-browser-darwin-x64/subPCR.txt")
samplerow=('A','B','C','D','E','F','G','H')
flag=False
for line in f:
    #line=line.strip('\n').lstrip()
    
    #sep=line.strip('\n').lstrip().split(' ')
    #count=sep[0]
    count="11"
    seq=line
    #print count
    if int(count) < 10:
        #print "count <10"
        continue
    
    flagfind=0
    m=0
    n=0
    end=""
    start=""
    for rvp in rv:
        if(seq.find(rvp)>=0):
            head=">"+str(samplerow[n])
            flagfind=flagfind+1
            end=seq.find(rvp)+26
            print end
        n=n+1
    for fwp in fw:
        m=m+1
        if flagfind == 1:
            if(seq.find(fwp)>=0):
                flagfind=flagfind+1
                head=head+str(m)
                start=seq.find(fwp)
    if flagfind == 2:
        #print start,end
        #seq=seq[start:end]
        #print len(seq)
        if seq.find("TCGACTTCGCCGCGC")>=0:
            
            #head=head+"\tABP7PCR"
            if seq.find(a7lcri)>=0:
                head=head+"\tA7-NE\t"+count
                print head
                print seq
            else:
                head=head+"\tA7-Edited\t"+count
                print head
                print seq
                #blast(seq,a7cri)
        if seq.find("TCAGCTTCGCCGCGC")>=0:
            #head=head+"\tABP7LikePCR"
            if seq.find(a7lcri)>=0:
                head=head+"\tA7L-NE\t"+count
                print head
                print seq
            else:
                head=head+"\tA7L-Edited\t"+count
                print head
                print seq
                #blast(seq,a7lcri)
        
        
        


f.close()

