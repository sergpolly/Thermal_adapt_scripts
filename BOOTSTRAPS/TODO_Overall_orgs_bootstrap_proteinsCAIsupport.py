import pandas as pd
import os
import subprocess as sub
import re
import sys
from Bio import SeqUtils
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st


root_path = "."
# # archaea/genomes
# data_path = os.path.join(root_path,"ftp.ncbi.nih.gov/genomes/all")
# asm_path = os.path.join(root_path,"archaea_links")
# result_path = os.path.join("..","archaea250")

# we'll calculate protein level CAI and store em in individual files ...
path_CAI = os.path.join(root_path,'genomes_with_CAI')

# # read table with organisms and their taxonid ...
dat = pd.read_csv(os.path.join(root_path,"assembly_description.dat"))

aacids = sorted(list('CMFILVWYAGTSNQDEHRKP'))


def init_iter_dat():
    iteration_dat = {}
    # we'll populate it with GC,topt, aacids - everything that's needed for the thremal adaptation trends calculation ...
    iteration_dat['GC'] = []
    iteration_dat['topt'] = []
    for aa in aacids:
        iteration_dat[aa] = []
    return iteration_dat


# consider organisms with the detailed protein description available ...

##### excluding Halophiles and empty entries ...
indcs = (dat.subdivision!='Halobacteria')&(dat.subdivision!='Nanohaloarchaea')
#####
valid_dat_subset = dat[indcs]
# taking into account only those with detailed protein info ...
valid_dat_subset = valid_dat_subset[valid_dat_subset['protein_details']]
# valid_dat_subset = dat[dat['protein_details']]
#####
# reset index just in case ...
valid_dat_subset = valid_dat_subset.reset_index()


PERCENTILE = 0.1
# FRACTION = 0.4
num_iterations = 30
slopes_generated = {}
# ribocheck = re.compile('ribosomal protein')
#####################################
for iteration in xrange(num_iterations):
    # sample_indicies = np.random.choice(valid_dat_subset.index,int(dat_size*FRACTION))
    # let's create local tiny copy of the dat DataFrame ...
    # we'll populate it with GC,topt, aacids - everything that's needed for the thremal adaptation trends calculation ...
    iteration_dat = init_iter_dat()
    print
    print "iteration %d"%iteration
    print
    for asm, topt in valid_dat_subset[['assembly','topt']].get_values():
        # open a file with the analysed organismal proteins...
        protein_fname = os.path.join(path_CAI,'%s_genes.dat'%asm)
        # load the data ...
        protein_dat = pd.read_csv(protein_fname)
        protein_dat_size = protein_dat.index.size
        subsample_size = int(protein_dat.index.size*PERCENTILE)
        ####################
        # now instead of taking proeins with high CAI, will be taking random subset of proteins ...
        prot_sample_indicies = np.random.choice(protein_dat.index,subsample_size)
        prot_subsample = protein_dat.loc[prot_sample_indicies]
        # print "taking %d random proteins from some organism: average cai (total vs subs): %.2f vs %.2f"%(subsample_size,protein_dat['cai'].mean(),prot_subsample['cai'].mean())
        cai_proteome = ''.join(prot_subsample['prot_seq'])
        cai_proteome_len = float(len(cai_proteome))
        cai_genome = ''.join(prot_subsample['gene_seq'])
        #
        iteration_dat['topt'].append(topt)
        iteration_dat['GC'].append(SeqUtils.GC(cai_genome))
        #
        for aa in aacids:
            # ACHTUNG ACHTUNG ACHTUNG ACHTUNG ...
            # multiply by 100.0 or not ?!
            iteration_dat[aa].append(100.0*cai_proteome.count(aa)/cai_proteome_len)
    # now, as soon as the iteration_dat is ready to go, let's just get amino acid trends out of it ...
    slopes_generated[iteration] = []
    for aa in aacids:
        a,b,r,pval,_ = st.linregress(iteration_dat['topt'],iteration_dat[aa])
        slopes_generated[iteration].append(a)



# print counter,'out of ',totcount
slopes_generated = pd.DataFrame(slopes_generated,index=aacids)
slopes_generated.to_csv('BOOTSTRAPED_prots_allORGS_supportCAI_arch.dat')


# plt.clf(); plt.imshow(slopes_generated.get_values(),interpolation='nearest'); plt.show()

















































