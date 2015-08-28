import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
# import numpy as np
import pandas as pd
# import itertools
from functools import partial

from multiprocessing import Pool


aacids = list('CMFILVWYAGTSNQDEHRKP')

def get_aausage_proteome(seqrec):
    # seqrec = db[seqrec_id]
    features = seqrec.features
    proteome = []
    for feature in features:
        qualifiers = feature.qualifiers
        if (feature.type == 'CDS')and('translation' in qualifiers):
            proteome.append(qualifiers['translation'][0])
    #return the results ...
    proteome = ''.join(proteome)
    prot_len = float(len(proteome))
    aa_freq = tuple(proteome.count(aa)/prot_len for aa in aacids)
    #
    return aa_freq

def analyse_genome(db,seqrec_id):
    seqrec = db[seqrec_id]
    aa_freq = get_aausage_proteome(seqrec)
    # gc = SeqUtils.GC(seqrec.seq)
    gc = 0.0
    id = seqrec.id
    return (id,gc) + aa_freq

def analyse_genome_serial(seqrec):
    aa_freq = get_aausage_proteome(seqrec)
    # gc = SeqUtils.GC(seqrec.seq)
    gc = 0.0
    id = seqrec.id
    return (id,gc) + aa_freq


if __name__ == "__main__":
    path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
    #
    gb_fname = 'condensed.raw.gb'
    gbs = SeqIO.parse(os.path.join(path,gb_fname),'genbank')
    #
    print "100% serial implementation ..."
    result = ( analyse_genome_serial(seqrec) for seqrec in gbs )
    #
    print "file outputting ..."
    with open("res.serial.xxx","w") as fp:
        for res in result:
            fp.write(str(res)+'\n')





##########################
#   it took 17 minutes !!!!!!!!!!!!!!!!!
############################










