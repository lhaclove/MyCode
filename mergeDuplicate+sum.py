dic = dict()

with open(r"E:\LPF\2017\KGL-ratio-1.gene2effect.np","r") as f:
    data = f.readlines()  
    for line in data:
        line.strip().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '').strip()
        key=line.split("\t")[0]
        value=line.split("\t")[1]
        if key in dic.keys():
            dic[key].append(value)
        else:
            dic[key] = [value]
    print(dic)
f6 = open(r"E:\LPF\2017\KGL-ratio-1.gene2effect.final.np",'w+')
for key in dic:
    #f6.write(key+'\t'+ ' '.join(str(i).replace('\n',"") for i in dic[key])+'\n')
    sum=0
    for i in dic[key]:
        sum+=abs(float(i))
    f6.write(key+'\t'+ str(sum)+'\n')
f6.close()