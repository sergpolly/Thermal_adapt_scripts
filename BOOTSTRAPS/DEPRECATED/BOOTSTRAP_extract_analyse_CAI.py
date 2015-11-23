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

#MAYBE WE DONT NEED THAT ANYMORE ...
# # STUPID FIX TO AVOID OLDER PANDAS HERE ...
# # PYTHONPATH seems to be ignored by th ipython ...
# sys.path.insert(1,"/home/venevs/.local/lib/python2.7/site-packages/")
# # this is needed just to avoid BUG in pandas (float indexing related: https://github.com/pydata/pandas/issues/5824)
# # when tsking quantile(q=0.75) ...
# import scipy.stats as stat

# RIBO_LIMIT = 24
PROT_LIMIT = 15


path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
dat = pd.read_csv(os.path.join(path,"complete_CDS.dat"))

plot_path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/plots')

# first: identify ribosomal proteins ...

# # here is our heuristic way to check if it's a ribosomal protein or not, given corresponding gene's product description ...
# ribo = re.compile("ribosomal.+protein",re.I)
# ribo_check = lambda line: bool(ribo.search(line)) if not('transferase' in line) else False

# either random or some other - non-ribo group of proteins ...
prot = re.compile("protein",re.I)
prot_check = lambda line: bool(prot.search(line))
dat['the_prot'] = dat['product'].apply(prot_check)


# based on these identified proteins, then calculate CAI ....

# group the data by the GenomicId ...
orgs = dat.groupby('GenomicID')

genom_id = orgs.groups.keys()

the_prot_counts = [(idx,orgs.get_group(idx)['the_prot'].nonzero()[0].size) for idx in genom_id]

prot_cai_info = pd.DataFrame(the_prot_counts,columns=['GenomicID','the_prot_count'])

print "Reference cDNA detected, min: %d, max %d"%(prot_cai_info.the_prot_count.min(),prot_cai_info.the_prot_count.max())

# some lists to describe organism's CAI distribution features ...
percentile = []
median = []
mean = []
sigma = []
idx_for_prot = []
prot_count_for_df = [] 
#
pid_cai_list = []
for idx,prot_count in prot_cai_info.itertuples(index=False):
    if prot_count >= PROT_LIMIT:
        cds_dat = orgs.get_group(idx)
        prot_cds = cds_dat[cds_dat['the_prot']]['cDNA'] # cDNA of protsomal proteins ...
        codon_usage = cairi.count_codons(prot_cds)
        codon_index = cairi.generate_codon_index(codon_usage,genetic_table=list(cds_dat['table'])[0]) # fix that ...
        # we need to track index from 'dat', as there are some stupid duplications ...
        pid_cai = pd.DataFrame(((dat_idx,pid,cairi.cai_for_gene(sequence,codon_index)) for dat_idx,pid,sequence in cds_dat[['pid','cDNA']].itertuples()),columns=['dat_idx','pid','CAI'])
        pid_cai = pid_cai.set_index(keys='dat_idx')
        # characterize CAI distribution for a given organism ...
        local_mean = pid_cai['CAI'].mean()
        local_median = pid_cai['CAI'].median()
        local_sigma = pid_cai['CAI'].std()
        mean.append(local_mean)
        median.append(local_median)
        sigma.append(local_sigma)
        idx_for_prot.append(idx)
        prot_count_for_df.append(prot_count)
        #
        local_prot_indexes = cds_dat['the_prot'].nonzero()[0]
        local_prot = pid_cai.iloc[local_prot_indexes].reset_index(drop=True) 
        # let's also check our t.o. score
        qH_all = pid_cai['CAI'].quantile(q=0.75)
        qL_rib = local_prot['CAI'].quantile(q=0.25)
        percentile.append( bool(qL_rib >= qH_all) )
        #
        # # OPTIONAL HISTOGRAM PLOTTING ...
        # # # # let's also plot histograms ...
        # xxx = dat_with_cai_trop['CAI'][dat_with_cai_trop.GenomicID=='NC_014500.1']
        # plt.clf()
        # plt.hist(xxx,range=(0,1),bins=100,color='blue',alpha=1.0)
        # # plt.hist(local_prot['CAI'],range=(0,1),bins=25,color='red',alpha=0.8)
        # # plt.title("%s, CAI median: %.2f, CoV %.3f, t.o. %s"%(idx,local_median,local_sigma/local_mean,str(qL_rib >= qH_all)))
        # # plt.savefig(os.path.join(plot_path,idx+".pdf"))
        #
        pid_cai_list.append( pid_cai )
# ttt = ["30S protsomal subunit protein S9", "protsomal-protein-alanine acetyltransferase", "protsomal protein L33", "protsomal subunit interface protein", "protsomal protein S10", "protsomal 5S rRNA E-loop binding protein Ctc/L25/TL5", "protsomal-protein-alanine acetyltransferase", "16S protsomal RNA methyltransferase KsgA/Dim1 family protein", "30S protsomal proteinS16", "Acetyltransferases including N-acetylases of protsomal proteins"]


org_cai_descr = {"GenomicID":idx_for_prot,"prot_count":prot_count_for_df,"TrOp":percentile,"median_cai":median,"mean_cai":mean,"sigma_cai":sigma}
org_cai_df = pd.DataFrame(org_cai_descr)

pid_cai_df = pd.concat(pid_cai_list)
#
# # before any mergings  ...
# ###########################################
# # MERGE BY THE INDEX  .... TO BE CONTINUED ...
# ###########################################
# # 1) merging 
# yyy = dat.join(pid_cai_df,lsuffix='',rsuffix='_wnans')#
# # 2) merging orther way ...
# xxx = pd.concat([dat,pid_cai_df],axis=1)
# #
# indexes = (xxx.CAI != yyy.CAI).nonzero()[0]
# # beware  (np.nan==np.nan) is False  ...
# # so there are ~1200 indexes ...
# # TO BE CONTINUED ...
# # merging is done, outputtting and that's it ...
dat_with_cai = dat.join(pid_cai_df,lsuffix='',rsuffix='_wnans')
# then simple check ...
# all instances, where (pid != pid_wnans) must be NULL ...
if dat_with_cai.pid_wnans[dat_with_cai.pid!=dat_with_cai.pid_wnans].isnull().all():
    pass
else:
    print "ACHTUNG!!! All pid_wnans items whose (pid_wnans!=pid), must be NULL. Check"

########### let's try joining the 'org_cai_df' to the dat_with_cai as well, so that we'd be able to easily grab Trans.Optimized
########### organisms ...
dat_with_cai_trop = pd.merge(dat_with_cai, org_cai_df, how='left', on='GenomicID')
# apparently 'join' is a legacy procedure, so using 'merge' is encouraged instead!
# http://stackoverflow.com/questions/10114399/pandas-simple-join-not-working

# output CDS info with the calculated CAI ...
dat_with_cai_trop[['GenomicID','cDNA','fid','pid','product','protein','status','table','the_prot','CAI','TrOp']].to_csv(os.path.join(path,"BOOT_complete_CDS_CAI_DNA.dat"),index=False)
# ['GenomicID', 'cDNA', 'fid', 'pid', 'product', 'protein', 'status', 'table', 'the_prot', 'pid_wnans', 'CAI']
# ['GenomicID', 'cDNA', 'fid', 'pid', 'product', 'protein', 'status', 'table', 'the_prot', 'CAI']

# #
# # some characterization plotting ...
# plt.clf()
# org_cai_trop = org_cai_df[org_cai_df["TrOp"]]
# org_cai_notrop = org_cai_df[~org_cai_df["TrOp"]]
# trop_dots = plt.plot(org_cai_trop.median_cai,np.true_divide(org_cai_trop.sigma_cai,org_cai_trop.mean_cai),'ro',label='translational optimization')
# notrop_dots = plt.plot(org_cai_notrop.median_cai,np.true_divide(org_cai_notrop.sigma_cai,org_cai_notrop.mean_cai),'bo',alpha=0.8,label='No translational optimization')
# ax = plt.gca()
# ax.set_title("Organism level CAI: t.o. criteria comparison (Margalit vs ours)")
# ax.set_xlabel("median CAI")
# ax.set_ylabel("CAI coefficient of variation") # using plain sigma works worse ...
# ax.legend(loc='best')
# plt.savefig(os.path.join(path,"org_cai_to.pdf"))
# #
# # 
# plt.clf()
# size_dot = lambda x: 10 if 50<x<60 else 120
# plt.scatter(x=org_cai_df.median_cai,y=np.true_divide(org_cai_df.sigma_cai,org_cai_df.mean_cai),s=org_cai_df.ribo_count.apply(size_dot),c="blue",edgecolor=None)
# ax = plt.gca()
# ax.set_title("Organism level CAI: effect of # ribosomal proteins (no effect)")
# ax.set_xlabel("median CAI")
# ax.set_ylabel("CAI coefficient of variation") # using plain sigma works worse ...
# plt.savefig(os.path.join(path,"org_cai_ribonum.pdf"))































