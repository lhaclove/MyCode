filebase="D:\MyNasCould\code\zma00010.txt"
id="zma00010"
genelist=""
flag=False
with open(filebase,"r") as f:
    genes=f.readlines()
    for ll in genes:
        #print(ll)
        #genelist=genelist+"\t"+ll.split()
        if (ll.startswith ('GENE')):
            flag = True
        if (ll.startswith ('COMPOUND')):
            flag = False
            
        if flag:
            genelist=genelist+ll
            #print(ll)
    with open(id+"_gene.txt","w") as out:
        out.write(genelist)
        out.close( )
