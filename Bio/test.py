from urllib import request

aa = "https://www.uniprot.org/uniprot/A0A1P8AT43.xml"
request.urlretrieve(aa, "python.xml")