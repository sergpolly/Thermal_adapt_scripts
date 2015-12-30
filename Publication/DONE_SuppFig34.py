import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
import numpy as np
from functools import partial
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




def CAI_hist(ax,cds_grouped,unique_id,title="",legend=True,xlims=None,ylims=None):
    #
    all_cds = cds_grouped.get_group(unique_id)
    # extract needed data ...
    all_CAI = np.asarray(all_cds['CAI'])
    ribo_CAI= np.asarray(all_cds[all_cds['ribosomal']]['CAI'])
    #
    # let's also plot histograms ...
    bins_all    = np.linspace(0,1,num=100)
    bins_ribo   = np.linspace(0,1,num=25)
    # all CAI
    all_hist    = ax.hist(all_CAI,bins=bins_all,color='blue',alpha=0.95,label='CDS all proteins',histtype='stepfilled',lw=0)
    # ribo CAI
    ribo_hist   = ax.hist(ribo_CAI,bins=bins_ribo,color='red',alpha=0.8,label='ribosomal proteins',histtype='stepfilled',lw=0)
    #
    ax.set_title(title)
    ax.set_xlabel('CAI')
    if legend:
        ax.legend(loc='upper right',frameon=False)#,markeredgecolor='none')
    #
    # set ylims xlims ....
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    #
    # return handles ...
    return (all_hist,ribo_hist)


def get_axis(xbins = 2, ybins = 1):
    plt.clf()
    # create grid of subplots on a large figure canvas
    # share x&y axes among all plots
    fig, ax = plt.subplots(ybins, xbins, figsize=(7.5,ybins*7.5/xbins), sharex=True, sharey=True)
    # no space between subplots
    fig.subplots_adjust(hspace=0.0, wspace=0.08)
    # make some room for the axes' labels
    l,b,r,t = 0.05,0.12,0.98,0.92
    w, h = r-l, t-b
    fig.subplots_adjust(bottom=b, left=l, right=r, top=t)
    #
    # returning axis ...
    return (fig,ax)




def TrOp_hist(ax,cds_grouped,unique_ids,title="",legend=True,xlims=None,ylims=None):
    # # This is a bug potentially ...
    # pd.Series([np.nan, np.nan, np.nan]).all() = True
    # pd.Series([np.nan, np.nan, np.nan]).any() = False
    def get_one_trop(idx):
        org_cds = cds_grouped.get_group(idx)
        # check if TrOp ...
        # for a given organism(id) all TrOp values must be same
        trop_vals = org_cds['TrOp'].unique()
        assert trop_vals.size == 1
        # then just figure out TrOp value after unpacking ...
        trop, = trop_vals
        if (not trop)or(pd.isnull(trop)):
            # False or Nan, return False 
            return False
        elif trop == True:
            # if it's True just return ...
            return trop
        else:
            raise ValueError
    #
    #
    cai_means_TrOp = [ cds_grouped.get_group(idx)['CAI'].mean() for idx in unique_ids if get_one_trop(idx) ]
    cai_means_noTrOp = [ cds_grouped.get_group(idx)['CAI'].mean() for idx in unique_ids if not get_one_trop(idx) ]
    # let's also plot histograms ...
    bins = np.linspace(0,1,num=50)
    # TrOp CAI
    trop_hist   = ax.hist(cai_means_TrOp,bins=bins,color='red',alpha=0.95,label='CUS',histtype='stepfilled',lw=0)
    # no TrOp CAI
    notrop_hist = ax.hist(cai_means_noTrOp,bins=bins,color='blue',alpha=0.8,label='non-CUS',histtype='stepfilled',lw=0)
    #
    ax.set_title(title)
    ax.set_xlabel('mean organismal CAI')
    if legend:
        ax.legend(loc='upper left',frameon=False)#,markeredgecolor='none')
    #
    # set ylims xlims ....
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    #
    # return handles ...
    return (trop_hist,notrop_hist)




# PLOT SUPP FIGURE 3 -> TROP vs NOTROP comparison ...

fig,ax = get_axis()
# choose a couple GenomicID using 'arch_nohalo' and/or 'bact' dataframes ...
arch1 = 'GCA_000204585.1' # good notrop example
arch2 = 'GCA_000013725.1'
# arch3 = 'GCA_000017945.1'
# arch3 = 'GCA_000446015.1'
# arch3 = 'GCA_000018365.1' # nice ...
# arch3 = 'GCA_000221185.1' # nice ...
arch3 = 'GCA_000009965.1' # quite nice ...
# arch3 = 'GCA_000017225.1' # nice
# arch3 = 'GCA_000006175.2'  #nice
# arch3 = 'GCA_000024185.1' # nice
# arch3 = 'GCA_000016125.1' # nice
get_trop = lambda grouped,idx: grouped.get_group(idx)['TrOp'].unique()[0]
gta = lambda idx: get_trop(arch_cai_by_org,idx)
get_org = lambda dat,idx,uid: dat[dat[uid]==idx]['organism_name'].iloc[0]
goa = lambda idx: get_org(arch_nohalo,idx,'assembly_accession')
#
CAI_hist(ax[0],arch_cai_by_org,arch1,legend=False,xlims=None,ylims=None,title=goa(arch1))
CAI_hist(ax[1],arch_cai_by_org,arch3,legend=True,xlims=None,ylims=None,title=goa(arch3))
# CAI_hist(ax[0],arch_cai_by_org,arch1,legend=False,xlims=None,ylims=None,title='non Translationally Optimized')
# CAI_hist(ax[1],arch_cai_by_org,arch3,legend=True,xlims=None,ylims=None,title='Translationally Optimized')
# plt.title("%s, CAI median: %.2f, CoV %.3f, t.o. %s"%(idx,local_median,local_sigma/local_mean,str(qL_rib >= qH_all)))
plt.savefig("SuppFig3.pdf")




# PLOT SUPP FIGURE 4 -> MEAN CAI vs OUR25/75 TROP METRIC ...



fig,ax = get_axis()
#
arch_ids = arch_nohalo['assembly_accession']
bact_ids = bact['GenomicID']
#
TrOp_hist(ax[0],arch_cai_by_org,arch_ids,title="Archaea",legend=True,xlims=None,ylims=None)
TrOp_hist(ax[1],bact_cai_by_org,bact_ids,title="Bacteria",legend=False,xlims=None,ylims=None)
#
plt.savefig("SuppFig4.pdf")




















