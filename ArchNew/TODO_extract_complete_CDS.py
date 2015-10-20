import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
import numpy as np
import pandas as pd
import cairi
from multiprocessing import Pool



if __name__ == "__main__":
    #
    root_path = "."
    path = os.path.join(root_path,"ftp.ncbi.nlm.nih.gov/genomes/all")
    dat = pd.read_csv(os.path.join(root_path,'summary_organisms_interest.dat'))
    #
    # path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
    # # dbx = SeqIO.index_db(os.path.join(path,'genbank.idx'))
    # dbx = SeqIO.index_db(os.path.join(path,"subset.idx"))
    # # subset.idx is a previously indexed db for the single genbank file with complete genomes only (the one ~5GB).
    # dat = pd.read_csv(os.path.join(path,"env_catalog_compgenome.dat"))
    # # #
    def do_work(asm_tuple):
        asm_suffix,asm_access = asm_tuple
        #
        asm_suffix = asm_suffix.replace(' ','_')
        asm_root_name = '_'.join([asm_access,asm_suffix])
        asm_dir = os.path.join(path, asm_root_name)
        gb_fname = os.path.join(asm_dir,'_'.join([asm_root_name,'genomic.gbff']))
        # fa_fname = os.path.join(asm_dir,'_'.join([asm_root_name,'genomic.fna']))
        #  we already know that these files exist (moreso, we don't need *.fna - as all nucs are present in the genbanks)
        SeqRecs = SeqIO.parse(gb_fname,format='genbank')
        #
        # seqrec = dbx[seqrecid]
        # data = {"GenomicID":[],"fid":[],"pid":[],"cDNA":[],"protein":[],"product":[],"table":[],"status":[]}
        tmp_results = (cairi.extract_genes_features(seqrec,seqrec.seq) for seqrec in SeqRecs)
        to_return = pd.concat((pd.DataFrame(res) for res in tmp_results),ignore_index=True)
        to_return['assembly_accession'] = asm_access #the whole column of identical elements ...
        # pd.DataFrame(dat).to_csv(os.path.join(outpath,out(seqrecid)),index=False)
        return to_return
    #
    work = list(dat[['asm_name','assembly_accession']].itertuples(index=False))
    #
    processes = int(sys.argv[1])
    pool = Pool(processes=processes)
    print "launching %d processes, to do %d pieces of work ..."%(processes,len(work))
    results = pool.map(do_work, work)
    #
    print "parallel calculations are OVER!!!"
    print ""
    #
    final_df = pd.concat(results,ignore_index=True)
    print "file outputting ..."
    final_df.to_csv("complete_arch_CDS.dat",index=False)

# THIS ARE THE GENBANKS THAT GOT STUCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ['NC_006526.2', 'NC_010645.1', 'NC_015436.1']

# sr1 = dbx['NC_006526.2']
# sr2 = dbx['NC_010645.1']
# sr3 = dbx['NC_015436.1']

# fid = 11013 # it's a~4000 amino acid protein global alignment that it getting stuck ...

# trying to recover cDNA @ genbank NC_003888.3 feature 12136 ...






