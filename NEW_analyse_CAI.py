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


RIBO_LIMIT = 24


path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
dat = pd.read_csv(os.path.join(path,"complete_CDS.dat"))


# first: identify ribosomal proteins ...

# here is our heuristic way to check if it's a ribosomal protein or not, given corresponding gene's product description ...
ribo = re.compile("ribosomal.+protein",re.I)
ribo_check = lambda line: bool(ribo.search(line)) if not('transferase' in line) else False

dat['ribosomal'] = dat['product'].apply(ribo_check)


# based on these identified proteins, then calculate CAI ....

# group the data by the GenomicId ...
orgs = dat.groupby('GenomicID')

genom_id = orgs.groups.keys()

ribo_counts = [(idx,orgs.get_group(idx)['ribosomal'].nonzero()[0].size) for idx in genom_id]

ribo_counts = pd.DataFrame(ribo_counts,columns=['GenomicID','ribo_count'])

pid_cai_list = []
for idx in genom_id:
    cds_dat = orgs.get_group(idx)
    ribo_cds = cds_dat[cds_dat['ribosomal']]['cDNA'] # cDNA of ribosomal proteins ...
    codon_usage = cairi.count_codons(ribo_cds)
    codon_index = cairi.generate_codon_index(codon_usage,genetic_table=list(cds_dat['table'])[0]) # fix that ...
    pid_cai = ((pid,cairi.cai_for_gene(sequence,codon_index)) for pid,sequence in cds_dat[['pid','cDNA'])
    pid_cai_list.append( pd.DataFrame(pid_cai,columns=['pid','CAI']) )

#
pid_cai_df = pd.concat(pid_cai_list) 
# ttt = ["30S ribosomal subunit protein S9", "ribosomal-protein-alanine acetyltransferase", "Ribosomal protein L33", "ribosomal subunit interface protein", "Ribosomal protein S10", "ribosomal 5S rRNA E-loop binding protein Ctc/L25/TL5", "ribosomal-protein-alanine acetyltransferase", "16S ribosomal RNA methyltransferase KsgA/Dim1 family protein", "30S ribosomal proteinS16", "Acetyltransferases including N-acetylases of ribosomal proteins"]


































