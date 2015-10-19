import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
from functools import partial
import time

from multiprocessing import Pool

aacids = list('CMFILVWYAGTSNQDEHRKP')

# def get_aausage_proteome(seqrec):
#     # seqrec = db[seqrec_id]
#     features = seqrec.features
#     proteome = []
#     for feature in features:
#         qualifiers = feature.qualifiers
#         if (feature.type == 'CDS')and('translation' in qualifiers):
#             proteome.append(qualifiers['translation'][0])
#     #return the results ...
#     proteome = ''.join(proteome)
#     prot_len = float(len(proteome))
#     aa_freq = tuple(proteome.count(aa)/prot_len for aa in aacids)
#     #
#     return (int(prot_len),) + aa_freq

def analyse_genome(asm_suffix,asm_access,path):
    #
    asm_suffix = asm_suffix.replace(' ','_')
    #
    asm_root_name = '_'.join([asm_access,asm_suffix])
    asm_dir = os.path.join(path, asm_root_name)
    gb_fname = os.path.join(asm_dir,'_'.join([asm_root_name,'genomic.gbff']))
    fa_fname = os.path.join(asm_dir,'_'.join([asm_root_name,'genomic.fna']))
    #  we already know that these files exist (moreso, we don't need *.fna - as all nucs are present in the genbanks)
    SeqRecs = SeqIO.parse(gb_fname,format='genbank')
    #
    # start extraction of genome and proteome ...
    proteome = ''
    genome = ''
    for SeqRec in SeqRecs:
        first_feature =  SeqRec.features[0]
        # first feature is usually description of the SOURCE ...
        if (first_feature.type == 'source') and ('plasmid' in first_feature.qualifiers) :
            pass
            # skipping such a SeqRec, because it's simply a plasmid ...
        else:
            # count CDS features only ...
            proteome += ''.join(feature.qualifiers['translation'][0]
                                    for feature in SeqRec.features
                                    if (feature.type=='CDS'
                                            and 'translation' in feature.qualifiers))
            genome += str(SeqRec.seq)
        #
    # *1.0 to make it float and avoid 1/3=0
    aa_freq = tuple(proteome.count(aa)*1.0/len(proteome) for aa in aacids)
    gc = SeqUtils.GC(genome)
    return (asm_access,gc,len(proteome),) + aa_freq

if __name__ == "__main__":
    #
    #
    root_path = "."
    path = os.path.join(root_path,"ftp.ncbi.nlm.nih.gov/genomes/all")
    dat = pd.read_csv(os.path.join(root_path,'summary_organisms_interest.dat'))
    #
    #
    #
    def do_work(asm_tuple):
        asm_suffix,asm_access = asm_tuple
        return analyse_genome(asm_suffix,asm_access,path=path)
    #
    the_work = list(dat[['asm_name','assembly_accession']].itertuples(index=False))
    #
    print "launching processes, to do %d pieces of the work ..."%len(the_work)
    pool = Pool(processes=12)
    results = pool.map(do_work, the_work)
    #
    print "file outputting ..."
    with open(os.path.join(root_path,"proteome_arch.dat"),"w") as fp:
        fp.write("assembly_accession,GC,ProtLen,"+",".join(aacids)+"\n")
        for result in results:
            fp.write(','.join(str(item) for item in result)+'\n')




























