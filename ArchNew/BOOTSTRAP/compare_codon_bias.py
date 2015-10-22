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

RIBO_LIMIT = 24
PROT_COUNT = 100


path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
dat = pd.read_csv(os.path.join(path,"complete_CDS.dat"))

plot_path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/plots')



#####################################################################
# here is our heuristic way to check if it's a ribosomal protein or not, given corresponding gene's product description ...
ribo = re.compile("ribosomal.+protein",re.I)
ribo_check = lambda line: bool(ribo.search(line)) if not('transferase' in line) else False
dat['ribosomal'] = dat['product'].apply(ribo_check)
#####################################################################




# group the data by the GenomicId ...
orgs = dat.groupby('GenomicID')

genom_id = orgs.groups.keys()

ribo_counts = [(idx,orgs.get_group(idx)['ribosomal'].nonzero()[0].size) for idx in genom_id]
ribo_cai_info = pd.DataFrame(ribo_counts,columns=['GenomicID','ribo_count'])

#############
#
#
cix_prot = {}
cix_ribo = {}
#
#
#############
for idx,ribo_count in ribo_cai_info.itertuples(index=False):
    #
    cds_dat = orgs.get_group(idx)
    prot_cds_rnd = cds_dat['cDNA'].sample(PROT_COUNT) # cDNA sample proteins ...
    codon_usage_rnd = cairi.count_codons(prot_cds_rnd)
    codon_index_rnd = cairi.generate_codon_index(codon_usage_rnd,genetic_table=list(cds_dat['table'])[0]) # fix that ...
    cix_prot[idx] = codon_index_rnd
    #
    if ribo_count >= RIBO_LIMIT:
        pass
        ribo_cds = cds_dat[cds_dat['ribosomal']]['cDNA'] # cDNA of ribosomal proteins ...
        codon_usage = cairi.count_codons(ribo_cds)
        codon_index = cairi.generate_codon_index(codon_usage,genetic_table=list(cds_dat['table'])[0]) # fix that ...
        cix_ribo[idx] = codon_index

######################


###########################################################################################################
#  some testing 
############################
from Bio.SeqUtils import CodonUsage
import math
from Bio import Data
import matplotlib.gridspec as gridspec

good_gid = 'NC_004606.1'
bad_gid = 'NC_007332.1'


# # print those RSCUs ...
# for aa in SynonymousCodons:
#     aa_codons = SynonymousCodons[aa]
#     print "amino acid %s:"%aa
#     for codon in aa_codons:
#         print "    %s %.3f %.3f"%(codon, cix_prot[gid][codon], cix_ribo[gid][codon])


# # ugly plot for RSCU ...
# colors = ['red','green','blue','black','magenta','red','green','blue','black','magenta','red','green','blue','black','magenta','red','green','blue','black','magenta']
# i=0
# plt.clf()
# for aa in SynonymousCodons:
#     aa_codons = SynonymousCodons[aa]
#     x,y = [],[]
#     # print "amino acid %s:"%aa
#     for codon in aa_codons:
#         x.append(cix_ribo[gid][codon])
#         y.append(cix_prot[gid][codon])
#         # print "    %s %.3f %.3f"%(codon, cix_prot[gid][codon], cix_ribo[gid][codon])
#     plt.plot(x,y,'o--',color=colors[i])
#     i+=1




cix = cix_ribo[good_gid]
##############################################################################
genetic_code = Data.CodonTable.unambiguous_dna_by_id[11]
SynonymousCodons = dict([(aa,[]) for aa in genetic_code.protein_alphabet.letters])
# SynonymousCodons['STOP'] = genetic_code.stop_codons # STOP codons are excluded from analysis ...
for codon,aa in genetic_code.forward_table.iteritems():
    SynonymousCodons[aa].append(codon)
# now exclude amino acids coded by a single codons: M,W for standard genetic code ...
unambiguous_aacids = [aa for aa,codons in SynonymousCodons.iteritems() if len(codons)<2]
for aa in unambiguous_aacids:
    # if there is no ambiguity, there is no optimality ...
    _ = SynonymousCodons.pop(aa)

###################################
## TRATATATATA!!!!!!!!!!
## CODON BIAS VISUALIZATION SIMILAR TO THE ONE HERE http://www.ijbs.com/v08p0093.pdf
## OR HERE: http://www.biomedcentral.com/1471-2164/13/276/figure/F5?highres=y

# we need SynonymousCodons and The RSCU indexes (some kind of relative codon usage thing ...)
# amino_acids = SynonymousCodons.keys()
# codons = [SynonymousCodons[aa] for aa in amino_acids] ...
colors = ['red','green','blue','white','magenta','yellow','cyan','red','green','blue','black','magenta','red','green','blue','black','magenta','red','green','blue','black','magenta']


bars_num = len(SynonymousCodons)
indices = np.arange(bars_num)
width = 0.8

plt.clf()
fig = plt.figure(figsize=(16,9))
gs = gridspec.GridSpec(4,1)
# plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, **fig_kw)
# fig,axs = plt.subplots(nrows=2, ncols=1,figsize=(16,9))
# datax = axs[0]
# legax = axs[1]
datax = plt.subplot(gs[:3, :])
legax = plt.subplot(gs[3, :]) #create the second subplot, that MIGHT be there


# def codon_color(codon):
#     gc = codon.count("G")+codon.count("C")



for idx,aa in zip(indices,sorted(SynonymousCodons)):
    aa_codons = SynonymousCodons[aa]
    # normalize those RSCUs ...
    norm = sum(cix[codon] for codon in aa_codons)
    norm = len(aa_codons)/norm
    bottom,leg_bottom = 0,6
    for jdx,codon in enumerate(aa_codons):
        # color = 
        datax.bar(idx,cix[codon]*norm,width=width,color=colors[jdx],bottom=bottom)
        legax.bar(idx,0.5,color=colors[jdx],bottom=leg_bottom,edgecolor='none')
        legax.text(idx+width/2.0,leg_bottom+0.25,codon,horizontalalignment='center',verticalalignment='center')
        bottom += cix[codon]*norm
        leg_bottom -= 0.5
    #
legax.axis('off')
datax.set_ylabel('RSCU')
datax.set_ylim(0,6.2)
# plt.title('Scores by group and gender')
datax.set_xticks(indices+width/2.)
datax.set_xticklabels( sorted(SynonymousCodons) )
fig.savefig(r"/home/venevs/Dropbox (UMASS MED - BIB)/xxx.png")# plt.yticks(np.arange(0,81,10))
# plt.legend( (p1[0], p2[0]), ('Men', 'Women') )







#MISTAKE ............






































