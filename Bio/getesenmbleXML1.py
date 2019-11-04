#encoding='utf-8'
import re
import urllib.request

import sys
import os
import time
import shutil
from os.path import join, getsize
from bs4 import BeautifulSoup

def getHtml(url):
   
    User_Agent='Mozilla/5.0 (Windows NT 5.2) AppleWebKit/534.30 (KHTML, like Gecko) Chrome/12.0.742.122 Safari/534.30'
    req=urllib.request.Request(url)
    req.add_header('User-agent',User_Agent)
    page = urllib.request.urlopen(req)
    html = page.read().decode('utf-8')
    #html = html.decode('unicode', 'ignore')
    return html


myQueryURL = "https://www.uniprot.org/uniprot/"
myQueryList=['Q8RX37']
#f=open("id",'r')
#myQueryList= f.readlines()


pattern = re.compile(r'{.+?}')

pattern1 = re.compile(r'OrderedLocusNames=')
pattern2 = re.compile(r'SubName: Full=')
pattern3 = re.compile(r';')
if(1):
    for n in myQueryList:
        out=""
        outAC=""
        outDE=""
        outGN=""
        filename = n.strip() + ".txt"
        myQuery=myQueryURL+n.strip()+".txt"
        #print(myQuery)
        myPage = getHtml(myQuery)
        cont=str(myPage)
        all=cont.split("\n")
        for line in all:
            #if line[0:2]=="AC":
                #outAC=outAC+line.replace("AC   ","")
            if line[0:2]=="DE":
                outDE=outDE+line.replace("DE   ","")
            if line[0:2]=="GN":
                outGN=outGN+line.replace("GN   ","")
        out=n.strip()+"\t"+outDE+"\t"+outGN
        out.replace("SubName:","")  
        out.replace("Full=","")  
        out.replace("SubName: Full=","") 
        out.replace("Name=","") 
        out.replace(";",",")
        out = re.sub(pattern, '', out)
        out = re.sub(pattern1, '', out)
        out = re.sub(pattern2, '', out)
        out = re.sub(pattern3, ',', out)
        print (out)