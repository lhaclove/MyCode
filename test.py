import re
import urllib
import sys
import os
import time
import shutil
from os.path import join, getsize
from bs4 import BeautifulSoup
import socket

with open("E:\ZXH\database\新建文件夹\pathway_zma.txt","r") as f:
    data = f.readlines()  
    for line in data:
        id=line.split("\t")[0]
        print("processing "+id)
        genelist=""
        flag=False
        with open(id+".txt","r") as f:
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