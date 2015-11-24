import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
import numpy as np
from functools import partial
import time
from scipy import stats as st
import copy 
from multiprocessing import Pool

# for plotting ...
import matplotlib.pyplot as plt
import matplotlib as mpl
font = {'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
mpl.rc('font', **font)
# data lopading ...
#
# reset the seed  here ... (?!) not thread safe (?!)
np.random.seed()


root_path = os.path.expanduser('~')
bact_path = os.path.join(root_path,'GENOMES_BACTER_RELEASE69/genbank')
arch_path = os.path.join(root_path,'GENOMES_ARCH_SEP2015')

# SOME ARCHAEAL DATA ...
arch        = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest.dat'))
arch_nohalo = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest_no_halop.dat'))

# SOME BACTERIAL DATA ...
# complete genomes only ...
bact        = pd.read_csv(os.path.join(bact_path,'env_catalog_compgenome.dat'))

###############################################
# complete_CDS_CAI_DNA.dat same thing ...
arch_cai_fname = os.path.join(arch_path,"complete_arch_CDS_CAI_DNA.dat")
bact_cai_fname = os.path.join(bact_path,"complete_CDS_CAI_DNA.dat")

arch_cai = pd.read_csv(arch_cai_fname)
bact_cai = pd.read_csv(bact_cai_fname)

bact_cai_by_org = bact_cai.groupby('GenomicID')
arch_cai_by_org = arch_cai.groupby('assembly_accession')


aacids = sorted('ACDEFGHIKLMNPQRSTVWY')




FRACTION = 0.4
ITERATIONS = 25
PERCENTILE = 0.1
# #########################################################
# #########################################################
# #########################################################
# PERCENTILE = 0.1
# FRACTION = 0.4
# num_iterations = 50
# dat_size = valid_dat_subset.index.size
# slopes_generated = {}
# for iteration in xrange(num_iterations):
#     sample_indicies = np.random.choice(valid_dat_subset.index,int(dat_size*FRACTION))
#     # get the subsample here ...
#     subsample = valid_dat_subset.loc[sample_indicies]
#     # check:
#     print
#     print "check ..."
#     print subsample.index.size,subsample.index.get_values()
#     print subsample['topt'].min(),subsample['topt'].max(),subsample['topt'].mean()
# #########################################################
# #########################################################
# #########################################################



# generate 'dataset':
# select CDS translations (protein sequences) for a given criteria ...
# criteria includes CAI top 10%, all proteins, TrOp, noTrOp, ribosomal ...
# 2 types of criteria: organismal  and  CDS-level ...
def get_random_slopeset(all_cds,dat,uid_key,cds_criteria='all',org_criteria='random',random_trop_fraction=0.5,topt='OptimumTemperature'):
    #
    # ACHTUNG !!!
    # the way we treated organisms with too little ribosomal proteins (~<24), makes it hard
    # to distinguish between non-TrOp and the ones the former ones.
    # This can cause downstream 'nan' evaluations as a result ...
    # To avoid troubles: treat the TrOp==np.nan differently !
    def extract_aausage(uid,criteria,cds_trans_key='protein'):
        local_cds = all_cds.get_group(uid)
        # select protein according to the criteria ...
        if criteria == 'cai':
            # cai10 for cai related criteria ...
            cai10 = local_cds['CAI'].quantile(q=1.0 - PERCENTILE)
            selected_aa = local_cds[local_cds['CAI'] >= cai10][cds_trans_key]
        elif criteria == 'ribo':
            selected_aa = local_cds[local_cds['ribosomal']][cds_trans_key]
        elif criteria == 'cai_noribo':
            # cai10 for cai related criteria ...
            cai10 = local_cds['CAI'].quantile(q=1.0 - PERCENTILE)
            selected_aa = local_cds[(local_cds['CAI'] >= cai10)&(~local_cds['ribosomal'])][cds_trans_key]
        elif criteria == 'all':
            selected_aa = local_cds[cds_trans_key]
        else:
            raise ValueError('CDS criteria must be either cai,ribo,cai_noribo or all!')
        #
        selected_aa = ''.join(selected_aa)
        total_aacount = float(len(selected_aa))
        # return freuqencies of aa usage
        get_aacount_str = lambda seq: np.asarray([seq.count(aa) for aa in aacids])
        return get_aacount_str(selected_aa)*100.0/total_aacount
    #
    # for each dataset get some slopes info ...
    def get_slopes(dat,topt='OptimumTemperature'):
        # exp_T = dat[dat[topt]>=50][aacids].mean()
        # exp_M = dat[(dat[topt]>=20)&(dat[topt]<=40)][aacids].mean()
        # exp_A = dat[aacids].mean()
        exp_D = dat[aacids].apply(lambda x: x.cov(dat[topt])/dat[topt].var())
        # check exp_D, just in case ...
        exp_D_check = dat[aacids].apply(lambda x: st.linregress(dat[topt],x)[0])
        if ( np.abs(exp_D - exp_D_check) > 1e-7 ).any():
            raise ValueError('exp_D caluclation failed!')
        # returning ...
        return exp_D
    #######################
    # we are still required to if CAI can be determined for the organism ...
    # are there enuogh ribo proteins, or in other words - TrOp status mut be not null ...
    def is_not_null(idx):
        org_cds = all_cds.get_group(idx)
        return org_cds['TrOp'].notnull().all()

    ############################################################################
    # instead of calculating TrOp we shoud do random subsampling of organisms ...
    subsample_idx = np.random.choice(dat.index,int(dat.shape[0]*random_trop_fraction))
    # get the subsample here ...
    dat_subsample = copy.deepcopy(dat.loc[subsample_idx])
    # filter the nulls ...
    dat_subsample['notnull'] = [is_not_null(uid) for uid in dat_subsample[uid_key]]
    dat_subsample_filtered = dat_subsample[dat_subsample['notnull']]
    ##############
    # the only criteria that makes sense ...
    if org_criteria in ['trop','random']:
        the_aausage = [ extract_aausage(uid,criteria=cds_criteria) for uid in dat_subsample_filtered[uid_key] ]
        # it's a resources wasting, but the code is cleaner this way ...
        return_topt = dat_subsample_filtered[topt]
    else:
        raise ValueError("Organism criteria must be 'trop' or 'random'!")
    # transform the aausage 
    the_aausage = np.asarray(the_aausage)
    the_aausage = pd.DataFrame(the_aausage,columns=aacids)
    the_aausage[topt] = return_topt.reset_index(drop=True)
    #
    return get_slopes(the_aausage,topt=topt)
###################################################
# cds_crits = ['cai','ribo','cai_noribo','all']
# org_crits = ['random']
###################################################
# cds_criteria = 'cai'
cds_criteria = 'all'
org_criteria = 'random'
# for iteration in range(ITERATIONS):
arch_iter_slopes = []
bact_iter_slopes = []
for iteration in range(ITERATIONS):
    one_arch_iter_slope = get_random_slopeset(arch_cai_by_org,arch_nohalo,'assembly_accession',cds_criteria,org_criteria,FRACTION)
    one_bact_iter_slope = get_random_slopeset(bact_cai_by_org,bact,'GenomicID',cds_criteria,org_criteria,FRACTION)
    one_arch_iter_slope.name = 'iter%d'%iteration
    one_bact_iter_slope.name = 'iter%d'%iteration
    arch_iter_slopes.append(one_arch_iter_slope)
    bact_iter_slopes.append(one_bact_iter_slope)


arch_boostrat_slopes = pd.concat(arch_iter_slopes,axis=1)
bact_boostrat_slopes = pd.concat(bact_iter_slopes,axis=1)


arch_fname = "BOOTSTRAP_%s_org-%s_cds-%s.dat"%('arch',org_criteria,cds_criteria)
bact_fname = "BOOTSTRAP_%s_org-%s_cds-%s.dat"%('bact',org_criteria,cds_criteria)


# csv export's default header=True, index=True must work just fine ...
arch_boostrat_slopes.to_csv(arch_fname)
bact_boostrat_slopes.to_csv(bact_fname)



print "Bootstrap for CSD:%s ORG:%s is over! "%(cds_criteria,org_criteria)





