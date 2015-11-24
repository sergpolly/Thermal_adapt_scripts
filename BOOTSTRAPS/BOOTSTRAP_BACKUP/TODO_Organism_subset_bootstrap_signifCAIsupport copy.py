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

PERCENTILE = 0.1
FRACTION = 0.4
num_iterations = 50
dat_size = valid_dat_subset.index.size

slopes_generated = {}

for iteration in xrange(num_iterations):
    sample_indicies = np.random.choice(valid_dat_subset.index,int(dat_size*FRACTION))
    # get the subsample here ...
    subsample = valid_dat_subset.loc[sample_indicies]
    # check:
    print
    print "check ..."
    print subsample.index.size,subsample.index.get_values()
    print subsample['topt'].min(),subsample['topt'].max(),subsample['topt'].mean()
    # let's create local tiny copy of the dat DataFrame ...
    # we'll populate it with GC,topt, aacids - everything that's needed for the thremal adaptation trends calculation ...
    iteration_dat = init_iter_dat()
    for asm, topt in subsample[['assembly','topt']].get_values():
        # open a file with the analysed organismal proteins...
        protein_fname = os.path.join(path_CAI,'%s_genes.dat'%asm)
        # load the data ...
        protein_dat = pd.read_csv(protein_fname)
        ####################
        percPERCtot = protein_dat['cai'].quantile(q=(1.0-PERCENTILE))
        significant_proteins = protein_dat[protein_dat['cai']>=percPERCtot]
        cai_proteome = ''.join(significant_proteins['prot_seq'])
        cai_proteome_len = float(len(cai_proteome))
        cai_genome = ''.join(significant_proteins['gene_seq'])
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


slopes_generated = pd.DataFrame(slopes_generated,index=aacids)
slopes_generated.to_csv('BOOTSTRAPED_ORGS_signifCAI10_arch.dat')











