import re
import urllib.request
import sys
import os
import time
import shutil
from os.path import join, getsize
from bs4 import BeautifulSoup
import socket
socket.setdefaulttimeout(10)



url="http://rest.kegg.jp/get/"
with open("E:\ZXH\database\新建文件夹\pathway_zma.txt","r") as f:
    data = f.readlines()  
    for line in data:
        id=line.split("\t")[0]
        uri = url + id
        try:
            #req=urllib.request.Request(uri)
            #page = urllib.request.urlopen(req)
            #html = page.read()
            urllib.request.urlretrieve(uri,id+".txt")
            print("processing "+id)
        except:
            continue