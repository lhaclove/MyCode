dic = dict()

with open("E:\\ZXH\\test\\gene_effect1.txt","r") as f:
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
f6 = open("E:\\ZXH\\test\\dict_th",'w+')
for key in dic:
    f6.write(key+'\t'+ ' '.join(str(i).replace('\n',"") for i in dic[key])+'\n')
f6.close()