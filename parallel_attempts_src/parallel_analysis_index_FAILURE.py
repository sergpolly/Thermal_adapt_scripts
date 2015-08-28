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
import time

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
    seqrec = db.get(seqrec_id)
    aa_freq = get_aausage_proteome(seqrec)
    # gc = SeqUtils.GC(seqrec.seq) # nucleotides are not available in some cases, so skip that ...
    gc = 0.0
    id = seqrec.id
    return (id,gc) + aa_freq


if __name__ == "__main__":
    start = time.time()
    path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
    # dbx = SeqIO.index_db(os.path.join(path,'genbank.idx'))
    dbx = SeqIO.index(os.path.join(path,"condensed.raw.gb"),"genbank")
    dat = pd.read_csv(os.path.join(path,"env_catalog_compgenome.dat"))
    stop = time.time()
    print "Indexing is done! time ",(stop-start)
    #
    #
    def do_work(seqrecid):
        return analyse_genome(dbx,seqrecid)
    #
    work = list(dat.GenomicID)
    #
    print "launching processes, to do %d pieces of work ..."%len(work)
    start = time.time()
    pool = Pool(processes=4)
    result = pool.map(do_work, work)
    stop = time.time()
    print "parallel work is done! time ",(stop-start)
    #
    print "file outputting ..."
    with open("res.index.xxx","w") as fp:
        for res in result:
            fp.write(str(res)+'\n')



