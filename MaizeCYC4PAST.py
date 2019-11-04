import re
import urllib
import sys
import os
import time
import shutil
from os.path import join, getsize
from bs4 import BeautifulSoup
import socket


rootdir = r'E:\ZXH\database\pathway\MaizeCyc'  
      
for parent,dirnames,filenames in os.walk(rootdir):    #三个参数：分别返回1.父目录 2.所有文件夹名字（不含路径） 3.所有文件名字
    #for dirname in  dirnames:                       #输出文件夹信息
        #print ("parent is:" + parent)
        #print  ("dirname is" + dirname)
   
    for filename in filenames:                        
        #print ("parent is" + parent)
        #print ("filename is:" + filename)
        #print ("the full name of the file is:" + os.path.join(parent,filename) )#输出文件路径信息
        with open(os.path.join(parent,filename),"r") as f:
            #print(filename.split(".")[0],end="")
            data = f.readlines()  
            out = filename.split(".")[0]+"\t"+data[0].replace("\n","")
            for line in data:
                if(line.startswith("GBWI")):
                    #gene += "\t"
                    out1 = out+"\t"+line.split("\t")[1]
                    print(out1)
       