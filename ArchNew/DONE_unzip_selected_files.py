from Bio import Entrez as ncbi
from Bio import SeqIO
import pandas as pd
import os
from time import sleep
import subprocess as sub
import re

# path = './ftp.ncbi.nih.gov/genomes/all'
path = './ftp.ncbi.nlm.nih.gov/genomes/all'
# ftp.ncbi.nlm.nih.gov

dirs = sub.check_output('ls %s'%path,shell=True).split()

# Files we are hunting for:
# _genomic.fna.gz (genomic fasta)
# _genomic.gbff.gz (GenBank with or without annotation)
# ...



for directory in dirs:
    path_dir = os.path.join(path,directory)
    retcode_genbank = sub.call('ls %s|grep genomic.gbff.gz$'%path_dir,shell=True)
    retcode_fasta =   sub.call('ls %s|grep genomic.fna.gz$'%path_dir,shell=True)
    # sub.call returns the return code of the shell command inside, it does not return actual shell command result ...
    # if grep fails to find stuff sub.call yields 1 (NO SUCCESS) ... 
    #
    if (retcode_genbank==0)and(retcode_fasta==0):
        genbank_zip = sub.check_output('ls %s|grep genomic.gbff.gz$'%path_dir,shell=True).split()
        fasta_zip =   sub.check_output('ls %s|grep genomic.fna.gz$'%path_dir,shell=True).split()
        if (len(genbank_zip) > 1)or(len(fasta_zip)>1):
            print "There are too many gbff or fna files in %s! Teminate!"%path_dir
            sys.exit(1)
        else:
            # unpack from the list of a single element ...
            genbank_zip, = genbank_zip
            fasta_zip, = fasta_zip
        # unzip ...
        gb_retcode = sub.call('gunzip %s'%os.path.join(path_dir,genbank_zip),shell=True)
        fa_retcode = sub.call('gunzip %s'%os.path.join(path_dir,fasta_zip),shell=True)
        if (gb_retcode or fa_retcode):
            print "Unable to unpack fna or gbff file in %s"%path_dir
            sys.exit(1)
    else:
        print "No archive either for gbff or fna file in %s"%path_dir

#
#
print "All genomic gbff and fna archives in subfolders of %s, were successfully unpacked!!!"%path





























