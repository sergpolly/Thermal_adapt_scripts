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
import matplotlib.pyplot as plt
from Bio.SeqUtils import CodonUsage
from Bio import Data


def back_translate(protein,SynCodon,CodonWeight):
    # use extracted codon bias to back-translate given proteins ... 
    # np.random.choice([0,1,2],size=100000,p=[1/6.,1/3.,0.5])
    # conditional is here to avoid 'illegal' amino acids, while preserving alignments ('XXX' codon insertion ...)
    return ''.join(np.random.choice(SynCodon[aa],p=CodonWeight[aa]) if aa in SynCodon else 'XXX' for aa in protein)
# # use extracted codon bias to back-translate given proteins ... 
# # np.random.choice([0,1,2],size=100000,p=[1/6.,1/3.,0.5])
# # conditional is here to avoid 'illegal' amino acids, while preserving alignments ('XXX' codon insertion ...)
# back_translate = lambda prot: ''.join(np.random.choice(SynonymousCodons[aa],p=codon_weights[aa]) if aa in SynonymousCodons else 'XXX' for aa in prot)
#

#MAYBE WE DONT NEED THAT ANYMORE ...
# # STUPID FIX TO AVOID OLDER PANDAS HERE ...
# # PYTHONPATH seems to be ignored by th ipython ...
# sys.path.insert(1,"/home/venevs/.local/lib/python2.7/site-packages/")
# # this is needed just to avoid BUG in pandas (float indexing related: https://github.com/pydata/pandas/issues/5824)
# # when tsking quantile(q=0.75) ...
# import scipy.stats as stat


if __name__ == "__main__":
    PROT_COUNT = 1000
    path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
    dat = pd.read_csv(os.path.join(path,"complete_CDS.dat"))
    # GenomicID,cDNA,fid,pid,product,protein,status,table
    # group the data by the GenomicId ...
    orgs = dat.groupby('GenomicID')
    genom_ids = orgs.groups.keys()
    #############
    #
    shuffled_cdna_list = []
    #
    work = genom_ids 
    #
    def do_work(idx):
        # IMPORTANT!!! 
        np.random.seed()
        #
        cds_dat = orgs.get_group(idx)
        #
        genetic_table = cds_dat['table'].iloc[0]
        genetic_code = Data.CodonTable.unambiguous_dna_by_id[genetic_table]
        SynonymousCodons = dict([(aa,[]) for aa in genetic_code.protein_alphabet.letters])
        # SynonymousCodons['STOP'] = genetic_code.stop_codons # STOP codons are excluded from analysis ...
        for codon,aa in genetic_code.forward_table.iteritems():
            SynonymousCodons[aa].append(codon)
        #
        #
        # prot_cds_rnd = cds_dat['cDNA'].sample(PROT_COUNT) # cDNA sample proteins ...
        prot_cds_rnd = cds_dat['cDNA'] # let's use ALL cDNA to get the codon bias ...
        codon_usage = cairi.count_codons(prot_cds_rnd)
        #
        # generate codon weights based on codon counts ...
        codon_weights = {}
        for aa in SynonymousCodons:
            aa_codon_usage = [ codon_usage[codon] for codon in SynonymousCodons[aa] ]
            total_codons = sum(aa_codon_usage)
            # normalize ...
            codon_weights[aa] = np.true_divide(aa_codon_usage,float(total_codons))
        #
        # now rewrite (back translate) protein sequences keeping their index ...
        cdna_shuffled =  ( (ix,pid,back_translate(protein,SynonymousCodons,codon_weights))  for ix,pid,protein in cds_dat[['pid','protein']].itertuples() )
        cdna_shuffled = pd.DataFrame(cdna_shuffled,columns=['ix','pid','cDNA_rnd'])
        cdna_shuffled = cdna_shuffled.set_index(keys='ix')
        #
        #
        # shuffled_cdna_list.append( cdna_shuffled )
        return cdna_shuffled
    #
    #
    print "launching processes, to do %d pieces of work ..."%len(work)
    pool = Pool(processes=16)
    results = pool.map(do_work, work)
    #
    print "file outputting ..."
    #
    shuffled_cdna_df = pd.concat(results)
    dat_cDNA_rnd = pd.merge(dat, shuffled_cdna_df, how='left')
    dat_cDNA_rnd.to_csv(os.path.join(path,"complete_CDS_Rnd.dat"),index=False)



# #############
# # SERIAL ...
# #############
# for idx in genom_ids:
#     cds_dat = orgs.get_group(idx)
#     #
#     genetic_table = cds_dat['table'].iloc[0]
#     genetic_code = Data.CodonTable.unambiguous_dna_by_id[genetic_table]
#     SynonymousCodons = dict([(aa,[]) for aa in genetic_code.protein_alphabet.letters])
#     # SynonymousCodons['STOP'] = genetic_code.stop_codons # STOP codons are excluded from analysis ...
#     for codon,aa in genetic_code.forward_table.iteritems():
#         SynonymousCodons[aa].append(codon)
#     #
#     #
#     prot_cds_rnd = cds_dat['cDNA'].sample(PROT_COUNT) # cDNA sample proteins ...
#     codon_usage = cairi.count_codons(prot_cds_rnd)
#     #
#     # generate codon weights based on codon counts ...
#     codon_weights = {}
#     for aa in SynonymousCodons:
#         aa_codon_usage = [ codon_usage[codon] for codon in SynonymousCodons[aa] ]
#         total_codons = sum(aa_codon_usage)
#         # normalize ...
#         codon_weights[aa] = np.true_divide(aa_codon_usage,float(total_codons))
#     #
#     # now rewrite (back translate) protein sequences keeping their index ...
#     cdna_shuffled =  ( (ix,pid,back_translate(protein,SynonymousCodons,codon_weights))  for ix,pid,protein in cds_dat[['pid','protein']].itertuples() )
#     cdna_shuffled = pd.DataFrame(cdna_shuffled,columns=['ix','pid','cDNA_rnd'])
#     cdna_shuffled = cdna_shuffled.set_index(keys='ix')
#     #
#     #
#     shuffled_cdna_list.append( cdna_shuffled )
# # GenomicID,cDNA,fid,pid,product,protein,status,table
# ######################
# #############
# # SERIAL ...
# #############


#
# dat_with_cai_trop = pd.merge(dat_with_cai, org_cai_df, how='left', on='GenomicID')
#
# dat_with_cai = dat.join(pid_cai_df,lsuffix='',rsuffix='_wnans')
# # then simple check ...
# # all instances, where (pid != pid_wnans) must be NULL ...
# if dat_with_cai.pid_wnans[dat_with_cai.pid!=dat_with_cai.pid_wnans].isnull().all():
#     pass
# else:
#     print "ACHTUNG!!! All pid_wnans items whose (pid_wnans!=pid), must be NULL. Check"
#
# ########### let's try joining the 'org_cai_df' to the dat_with_cai as well, so that we'd be able to easily grab Trans.Optimized
# ########### organisms ...
# dat_with_cai_trop = pd.merge(dat_with_cai, org_cai_df, how='left', on='GenomicID')
# # apparently 'join' is a legacy procedure, so using 'merge' is encouraged instead!
# # http://stackoverflow.com/questions/10114399/pandas-simple-join-not-working
#
# # output CDS info with the calculated CAI ...
# dat_with_cai_trop[['GenomicID','cDNA','fid','pid','product','protein','status','table','ribosomal','CAI','TrOp']].to_csv(os.path.join(path,"complete_CDS_CAI_DNA.dat"),index=False)
#


#
# ##############################################################################
# genetic_code = Data.CodonTable.unambiguous_dna_by_id[11]
# SynonymousCodons = dict([(aa,[]) for aa in genetic_code.protein_alphabet.letters])
# # SynonymousCodons['STOP'] = genetic_code.stop_codons # STOP codons are excluded from analysis ...
# for codon,aa in genetic_code.forward_table.iteritems():
#     SynonymousCodons[aa].append(codon)
# # now exclude amino acids coded by a single codons: M,W for standard genetic code ...
# unambiguous_aacids = [aa for aa,codons in SynonymousCodons.iteritems() if len(codons)<2]
# for aa in unambiguous_aacids:
#     # if there is no ambiguity, there is no optimality ...
#     _ = SynonymousCodons.pop(aa)
#
#









