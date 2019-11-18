
#encoding='utf-8'
import re
import urllib.request

import sys
import os
import time
import shutil
from os.path import join, getsize
from bs4 import BeautifulSoup
myQueryURL = "http://www.t66y.com/thread0806.php?fid=5&type=2"
HomePage = "http://www.t66y.com/"
def getHtml(url):
   
    User_Agent='Mozilla/5.0 (Windows NT 5.2) AppleWebKit/534.30 (KHTML, like Gecko) Chrome/12.0.742.122 Safari/534.30'
    req=urllib.request.Request(url)
    req.add_header('User-agent',User_Agent)
    page = urllib.request.urlopen(req)
    html = page.read()
    
    #html = html.decode('unicode', 'ignore')
    return html

myPage = getHtml(myQueryURL)
soup = BeautifulSoup(myPage,"html.parser")
table = soup.find_all('h3')
for tab in table:
    link = tab.find_all('a')
    for url in link:
        url = url.get('href')
        print(HomePage+url)