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

def getXml(url,fn):
   
    User_Agent='Mozilla/5.0 (Windows NT 5.2) AppleWebKit/534.30 (KHTML, like Gecko) Chrome/12.0.742.122 Safari/534.30'
    req=urllib.request.Request(url)
    req.add_header('User-agent',User_Agent)
    page = urllib.request.urlretrieve(url,fn)
    return page


myQueryURL = "https://www.uniprot.org/uniprot/"
#myQueryList=['O48741']
f=open("id",'r')

myQueryList= f.readlines()
if(1):
    for n in myQueryList:
        filename = n.strip() + ".txt"
        myQuery=myQueryURL+n.strip()+".txt"
        print(myQuery)
        myPage = getXml(myQuery,filename)
if(0):
    for n in myQueryList:
        filename = n.strip() + ".xml"
        soup = BeautifulSoup(filename, 'xml')
       





    

