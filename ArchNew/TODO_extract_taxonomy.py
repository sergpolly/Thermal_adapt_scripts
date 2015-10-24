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


def analyse_genome(asm_suffix,asm_access,path):
    # combine suffix and access id to form asm_root name ...
    get_asm_root = lambda sfx,idx: idx+'_'+sfx.replace(' ','_')
    asm_root = get_asm_root(asm_suffix,asm_access)
    # assembly directory and GenBank file in it ...
    gb_fname = os.path.join(path, asm_root, asm_root+'_genomic.gbff')
    #  we already know that these files exist (moreso, we don't need *.fna - as all nucs are present in the genbanks)
    SeqRecs = SeqIO.parse(gb_fname,format='genbank')
    #
    # all we need this time around is the taxonomy ...
    minimalTaxonomySet = set(';'.join(SeqRec.annotations['taxonomy']) for SeqRec in SeqRecs if 'taxonomy' in SeqRec.annotations)
    if not minimalTaxonomySet:
        print "There is no taxonomy entries for %s assembly"%asm_root
        sys.exit(1)
    elif len(minimalTaxonomySet) > 1:
        print "Assembly %s got more than 1 taxonomy descriptions!"%asm_root
        # turned out 37 out of 263 assmeblies have more than 1 taxonomy lineage ...
        # sys.exit(1)
        # there is not much we can do, but return all unique taxonomy lineages ...
        taxon_to_return = ':'.join(minimalTaxonomySet)
    else:
        taxon_to_return, = minimalTaxonomySet
    # returnring ...
    return (asm_access,taxon_to_return)

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
    with open(os.path.join(root_path,"arch_taxonomy_interest.dat"),"w") as fp:
        fp.write("assembly_accession,tax_lineages\n")
        for result in results:
            fp.write(','.join(str(item) for item in result)+'\n')




























