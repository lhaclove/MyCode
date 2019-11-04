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
    html = page.read()
    
    #html = html.decode('unicode', 'ignore')
    return html

myQueryURL = "http://www.uniprot.org/uniprot/?query="
#myQueryList=['Zm00001d001820']
f=open("id",'r')
myQueryList= f.readlines()
for n in myQueryList:
    myQuery=myQueryURL+n
    myPage = getHtml(myQuery)
    #print(myPage)
    soup = BeautifulSoup(myPage,"html.parser")
    for pname in soup.find_all(name='div',attrs={"class":"long"}):
        pname=pname.get_text()
        print (n.strip()+"\t"+pname)
        #print (pname)
    






    

