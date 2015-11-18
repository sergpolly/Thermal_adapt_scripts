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

from multiprocessing import Pool

# for plotting ...
import matplotlib.pyplot as plt
import matplotlib as mpl
font = {'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
mpl.rc('font', **font)
# data lopading ...



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


aacids = list('ACDEFGHIKLMNPQRSTVWY')




# generate 'dataset':
# select CDS translations (protein sequences) for a given criteria ...
# criteria includes CAI top 10%, all proteins, TrOp, noTrOp, ribosomal ...
# 2 types of criteria: organismal  and  CDS-level ...
def get_dataset(all_cds,dat,uid_key,cds_criteria='all',org_criteria='all',calculate_trop=False,topt='OptimumTemperature'):
    #
    # ACHTUNG !!!
    # the way we treated organisms with too little ribosomal proteins (~<24), makes it hard
    # to distinguish between non-TrOp and the ones the former ones.
    # This can cause downstream 'nan' evaluations as a result ...
    # To avoid troubles: treat the TrOp==np.nan differently !
    #
    #
    # # This is a bug potentially ...
    # pd.Series([np.nan, np.nan, np.nan]).all() = True
    # pd.Series([np.nan, np.nan, np.nan]).any() = False
    def get_one_trop(idx):
        org_cds = all_cds.get_group(idx)
        # check if TrOp ...
        # for a given organism(id) all TrOp values must be same
        trop_vals = org_cds['TrOp'].unique()
        assert trop_vals.size == 1
        # then just figure out TrOp value after unpacking ...
        trop, = trop_vals
        if pd.isnull(trop):
        	# special return - not enough ribosomal proteins ...
        	return 'none'
        if not trop:
            # False, return False 
            return 'false'
        elif trop == True:
            # if it's True just return ...
            return 'true'
        else:
            raise ValueError
    #
    #
    def extract_aausage(uid,criteria,cds_trans_key='protein'):
        local_cds = all_cds.get_group(uid)
        # select protein according to the criteria ...
        if criteria == 'cai':
            # cai10 for cai related criteria ...
            cai10 = local_cds['CAI'].quantile(q=1.0 - 0.1)
            selected_aa = local_cds[local_cds['CAI'] >= cai10][cds_trans_key]
        elif criteria == 'ribo':
            selected_aa = local_cds[local_cds['ribosomal']][cds_trans_key]
        elif criteria == 'cai_noribo':
            # cai10 for cai related criteria ...
            cai10 = local_cds['CAI'].quantile(q=1.0 - 0.1)
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
    # Calculate or check if we have the TrOp info for each organism ...
    if calculate_trop:
        dat['TrOp'] = [get_one_trop(idx) for idx in dat[uid_key]]
    else:
        assert 'TrOp' in dat.columns
    #
    #
    if org_criteria == 'all':
        return_aausage = [ extract_aausage(uid,criteria=cds_criteria) for uid in dat[dat['TrOp']!='none'][uid_key]]
        return_topt = dat[dat['TrOp']!='none'][topt]
    elif org_criteria == 'trop':
        return_aausage = [ extract_aausage(uid,criteria=cds_criteria) for uid in dat[dat['TrOp']=='true'][uid_key] ]
        return_topt = dat[dat['TrOp']=='true'][topt]
    else:
        raise ValueError('Organism criteria must be either all or trop!')
    # transform the aausage 
    return_aausage = np.asarray(return_aausage)
    return_aausage = pd.DataFrame(return_aausage,columns=aacids)
    return_aausage[topt] = return_topt.reset_index(drop=True)
    #
    return return_aausage




def plot_comparison(ax,x,y,xlab='x',ylab='y'):
    #
    ax.plot(x,y,'bo',visible=False)
    for i,aa in enumerate(aacids):
        ax.text(x[i], y[i], aa, transform=ax.transData)
    ###################################################
    a,b,r,pval,_ = st.linregress(x,y)
    #
    middle = 0.5*(x.min()+x.max())
    delta = 1.1*(x.max()-x.min())
    x_range = [middle-0.5*delta,middle+0.5*delta]
    #############################################
    fff = np.vectorize(lambda x: a*x + b)
    lin_fit, = ax.plot(x_range,fff(x_range),color='gray',linewidth=1.5,linestyle='-',zorder=100,label='linear fit, R=%.3f'%r)
    ###################################################
    ###################################################
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ###################################################
    ###################################################
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.legend((lin_fit,),('linear fit, R=%.3f'%r,),loc='best',frameon=False)
    plt.tight_layout(pad=0.4, h_pad=None, w_pad=None)
    ###################################################
    ###################################################
    # axis limits ...
    middle = 0.5*(x.min()+x.max())
    delta = 1.3*(x.max()-x.min())
    x_lims = [middle-0.5*delta,middle+0.5*delta]
    middle = 0.5*(y.min()+y.max())
    delta = 1.3*(y.max()-y.min())
    y_lims = [middle-0.5*delta,middle+0.5*delta]
    # ###################################################
    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    ###################################################
    # ax.locator_params(axis='both',tight=True)#nbins=10)
    # ax.locator_params(axis='x',nbins=5)
    ###################################################





def get_axis(xbins = 2, ybins = 3):
    plt.clf()
    # create grid of subplots on a large figure canvas
    # share x&y axes among all plots
    fig, ax = plt.subplots(ybins, xbins, figsize=(7.5,ybins*7.5/xbins))#, sharex=True, sharey=True)
    # # no space between subplots
    # fig.subplots_adjust(hspace=0.0, wspace=0.08)
    # # make some room for the axes' labels
    # l,b,r,t = 0.05,0.12,0.98,0.92
    # w, h = r-l, t-b
    # fig.subplots_adjust(bottom=b, left=l, right=r, top=t)
    # #
    # # returning axis ...
    return (fig,ax)



# aafreqs_topt = get_dataset(arch_cai_by_org,arch_nohalo,'assembly_accession',cds_criteria='all',org_criteria='all',calculate_trop=True,topt='OptimumTemperature')


cds_crits = ['cai','ribo','cai_noribo','all']
org_crits = ['trop','all']

combinations = [(cc,oo) for cc in cds_crits for oo in org_crits]
combinations = [[{'cds_criteria':cc,'org_criteria':oo} for cc in cds_crits] for oo in org_crits]
# # with trop
# combinations[0][0,2,3]
# # no trop ...
# combinations[1][0,1,3]

# ######################################
# for i,(for_trop,for_notrop) in enumerate(zip([0,2,3],[0,1,3])):
#     #
#     print i,combinations[0][for_trop]
#     #
#     kwargs = combinations[0][for_trop]
#     kwargs['calculate_trop'] = True
#     # {'cds_criteria':combinations[0][cc],'org_criteria':combinations[0][oo],'calculate_trop':True}
#     arch_slopes_x = get_slopes(get_dataset(arch_cai_by_org,arch_nohalo,'assembly_accession',**kwargs))
#     bact_slopes_y = get_slopes(get_dataset(bact_cai_by_org,bact,'GenomicID',**kwargs))
#     print arch_slopes_x,bact_slopes_y
#     # #
#     print i,combinations[1][for_notrop]
#     #
#     kwargs = combinations[1][for_notrop]
#     kwargs['calculate_trop'] = True
#     # kwargs = {'cds_criteria':combinations[1][cc],'org_criteria':combinations[1][oo],'calculate_trop':True}
#     arch_slopes_x = get_slopes(get_dataset(arch_cai_by_org,arch_nohalo,'assembly_accession',**kwargs))
#     bact_slopes_y = get_slopes(get_dataset(bact_cai_by_org,bact,'GenomicID',**kwargs))
#     print arch_slopes_x,bact_slopes_y
# #######################################


# # #TESTING ....
# #
# # kkk = {'cds_criteria': 'cai', 'org_criteria': 'trop', 'calculate_trop': True}
# kkk = {'cds_criteria': 'cai', 'org_criteria': 'all', 'calculate_trop': True}
# # xxx = get_dataset(arch_cai_by_org,arch_nohalo,'assembly_accession',**kkk)
# xxx = get_dataset(bact_cai_by_org,bact,'GenomicID',**kkk)
# #
# xxx[xxx['A'].isnull()] # there is 1 that evaluate to np.nan?!!!!!!!!!!!!!!!


# {'cds_criteria': 'ribo', 'org_criteria': 'all'}

# get slopes ...
def get_slopes(dat,topt='OptimumTemperature'):
    return np.asarray([st.linregress(dat[topt],dat[aa])[0] for aa in aacids])

# # draw the thing ...
fig,ax = get_axis()

for i,(for_trop,for_notrop) in enumerate(zip([0,2,3],[0,1,3])):
    #
    print i,combinations[0][for_trop]
    #
    kwargs = combinations[0][for_trop]
    kwargs['calculate_trop'] = True
    # {'cds_criteria':combinations[0][cc],'org_criteria':combinations[0][oo],'calculate_trop':True}
    arch_slopes_x = get_slopes(get_dataset(arch_cai_by_org,arch_nohalo,'assembly_accession',**kwargs))
    bact_slopes_y = get_slopes(get_dataset(bact_cai_by_org,bact,'GenomicID',**kwargs))
    with open("slopes_%s_%s.info"%(kwargs['cds_criteria'],kwargs['org_criteria']),'w') as fp:
        fp.write("aa,arch_a,bact_a\n")
        for aa,arch_a,bact_a in zip(aacids,arch_slopes_x,bact_slopes_y):
            fp.write("%s,%.4f,%.4f\n"%(aa,arch_a,bact_a))
    plot_comparison(ax[i,0],arch_slopes_x,bact_slopes_y,ylab='_'.join([str(_) for _ in kwargs.values()]),xlab='x')
    # #
    print i,combinations[1][for_notrop]
    #
    kwargs = combinations[1][for_notrop]
    kwargs['calculate_trop'] = True
    # kwargs = {'cds_criteria':combinations[1][cc],'org_criteria':combinations[1][oo],'calculate_trop':True}
    arch_slopes_x = get_slopes(get_dataset(arch_cai_by_org,arch_nohalo,'assembly_accession',**kwargs))
    bact_slopes_y = get_slopes(get_dataset(bact_cai_by_org,bact,'GenomicID',**kwargs))
    with open("slopes_%s_%s.info"%(kwargs['cds_criteria'],kwargs['org_criteria']),'w') as fp:
        fp.write("aa,arch_a,bact_a\n")
        for aa,arch_a,bact_a in zip(aacids,arch_slopes_x,bact_slopes_y):
            fp.write("%s,%.4f,%.4f\n"%(aa,arch_a,bact_a))
    plot_comparison(ax[i,1],arch_slopes_x,bact_slopes_y,ylab='_'.join([str(_) for _ in kwargs.values()]),xlab='x')



plt.show()


# # # # simple thing to extract 20 slopes per dataset ...
# # # with open('cai10_archaea.slopes_info','w') as fp:
# # #     fp.write("aa,a,b,r,pval\n")
# # #     for aa in aacids:
# # #         a,b,r,pval,_ = st.linregress(dat['topt'],dat[aa])
# # #         fp.write("%s,%.4f,%.4f,%.4f,%.e\n"%(aa,a,b,r,pval))









