
from Bio import Entrez as ncbi
from Bio import SeqIO
import pandas as pd
import os
from time import sleep
import subprocess as sub
import re



path = './ftp.ncbi.nih.gov/genomes/all'


dirs = sub.check_output('ls %s'%path,shell=True).split()


for directory in dirs:
    retcode_genbank = sub.call('ls %s|grep gbff.gz$'%os.path.join(path,directory),shell=True)
    if not retcode_genbank:
        genbank_zip = sub.check_output('ls %s|grep gbff.gz$'%os.path.join(path,directory),shell=True).split()
        if len(genbank_zip) > 1:
            print "There are too many gbff-s in %s! Teminate!"%os.path.join(path,directory)
            sys.exit(1)
        else:
            # unpack from the list ...
            genbank_zip = genbank_zip[0]
        # unzip ...
        retcode = sub.call('gunzip %s'%os.path.join(path,directory,genbank_zip),shell=True)
        if retcode:
            print "Unable to unpack gbff.gz: %s"%os.path.join(path,directory,genbank_zip)
            sys.exit(1)

print "All genbank archives in subfolders of %s, were successfully unpacked!!!"%path





























