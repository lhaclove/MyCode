import re
import urllib
import sys
import os
import time
import shutil
from os.path import join, getsize
from bs4 import BeautifulSoup
import socket
with open("out.txt","w+") as out:
    with open("E:\ZXH\database\新建文件夹\pathway_zma.txt","r") as f:
        data = f.readlines()  
        for line in data:
            genelist=""
            id=line.split("\t")[0]
            print("processing "+id)
            genelist=""
            flag=False
            with open(id+"_gene.txt","r") as f:
                genes=f.readlines()
                for ll in genes:
                    if(ll.startswith("GENE")):
                        genelist=genelist+ll.split()[1]+"\t"
                    else:
                        genelist=genelist+ll.split()[0]+"\t"
                out.write(id+"\t"+genelist+"\n")
out.close()
