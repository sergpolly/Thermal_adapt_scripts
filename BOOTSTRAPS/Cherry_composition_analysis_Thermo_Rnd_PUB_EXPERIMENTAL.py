import pandas as pd
import os
import subprocess as sub
import re
import sys
from Bio import SeqUtils
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy import stats as st

import matplotlib as mpl
#
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# #
#
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{textcomp}',   # i need upright \micro symbols, but you need...
       # r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]
#
font = {#'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
rc('font', **font)
# # data loading ...



# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
root_path = os.path.expanduser('~')
bact_path = os.path.join(root_path,'GENOMES_BACTER_RELEASE69/genbank')
arch_path = os.path.join(root_path,'GENOMES_ARCH_SEP2015')

# this is BACTERIA script, so just do 
path = bact_path

uid = 'GenomicID'
topt_id = 'OptimumTemperature'
aacids = sorted(list('CMFILVWYAGTSNQDEHRKP'))
the_num_of_quantiles = 5

#######################
def get_one_trop(all_cds_grouped,idx):
    org_cds = all_cds_grouped.get_group(idx)
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





def get_quantiles_summary(cds_cai_dat,num_of_quantiles,R20_vec_compare,vec_cost):
    # we can use this 'qcut' function from pandas to divide our proteins by the quantiles ...
    category,bins = pd.qcut(cds_cai_dat['CAI'],q=num_of_quantiles,retbins=True,labels=False)
    # then we could iterate over proteins/cDNAs in these categories ...
    fivywrel_cat, r20_cat, cost_cat = [],[],[]
    for cat in range(num_of_quantiles):
        cds_cai_category = cds_cai_dat[category==cat]
        protein_length_distro = cds_cai_category['protein'].str.len()
        # average protein length per quantile as a stability measure ...
        average_length = protein_length_distro.mean()
        # total proteins length in quantile for AA freqs calculations ...
        total_length = protein_length_distro.sum()
        IVYWREL = sum(cds_cai_category['protein'].str.count(aa).sum() for aa in list('IVYWREL'))
        # IVYWREL = cds_cai_category['protein'].str.count('|'.join("IVYWREL")).sum() # tiny bit slower ...
        f_IVYWREL = float(IVYWREL)/float(total_length)
        # 20-vector for of amino acid composition ...
        aa_freq_20 = np.true_divide([cds_cai_category['protein'].str.count(aa).sum() for aa in aacids],float(total_length))
        # slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        _1,_2,R20,_4,_5 = stats.linregress(aa_freq_20, R20_vec_compare)
        # Akashi ...
        cost = np.dot(aa_freq_20,vec_cost)
        # storing info ...
        fivywrel_cat.append(f_IVYWREL)
        r20_cat.append(R20)
        cost_cat.append(cost)
    #returning ...
    return (fivywrel_cat,r20_cat,cost_cat)






def quantile_summary_storage_init(num_of_quantiles=5):
    # data structure to store info ...
    the_storage = {uid:[],
                topt_id:[],
                'TrOp':[]}
    for i in range(num_of_quantiles):
        the_storage['q%d'%i] = []
        the_storage['R20_q%d'%i] = []
        the_storage['Akashi_q%d'%i] = []
        # the_storage['ProtLen_q%d'%i] = []
    # returning ...
    return the_storage






# ['DbxRefs','Description','FeaturesNum','assembly_accession','GenomicLen','GenomicName','Keywords','NucsPresent','Organism_des',
# 'SourceDbxRefs','SourceOrganism','SourcePlasmid','SourceStrain','Taxonomy','BioProject','TaxonID','Organism_env',
# 'OptimumTemperature','TemperatureRange','OxygenReq','Habitat','Salinity','crit_NC','crit_WGS','crit_genlen',
# 'crit_features','crit_comp_genome','crit_plasmid']
# env_dat = pd.read_csv(os.path.join(path,"summary_organisms_interest.dat"))
env_dat = pd.read_csv(os.path.join(path,'env_catalog_compgenome.dat'))
#['assembly_accession','cDNA','fid','pid','product','protein','status','table','ribosomal','CAI','TrOp']
gen_dat = pd.read_csv(os.path.join(path,"complete_CDS_CAI_DNA_Rnd.dat"))
#
original_cai_dat = pd.read_csv(os.path.join(path,"complete_CDS_CAI_DNA.dat"))
#
# PROTEOME LEVEL AMINO ACID FREQUENCIES ...
# "proteome_all.dat"
# # file with the organisms of interest 
# dat_fname = os.path.join(bib2_scr_path,'catalog_with_accesion.dat')
# dat = pd.read_csv(dat_fname)
#
cost_vec_path = path
akashi = os.path.join(cost_vec_path,'akashi-cost.d')
argentina = os.path.join(cost_vec_path,'argentina-cost.d')
#
akashi_cost = pd.read_csv(akashi,header=None,sep=' ')
argentina_cost = pd.read_csv(argentina,header=None,sep=' ')
thermo_freq = pd.read_csv(os.path.join(path,'thermo.dat'),header=None,sep=' ')
#
akashi_cost.set_index(0,inplace=True)
argentina_cost.set_index(0,inplace=True)
thermo_freq.set_index(0,inplace=True)
#
akashi_cost.sort_index(inplace=True)
argentina_cost.sort_index(inplace=True)
thermo_freq.sort_index(inplace=True)


# 
#####################################################
#  QUANTILE STATISTICALL DATA EXTRACTION ....
#####################################################
gen_dat_org = gen_dat.groupby(uid)
original_gen_dat_org = original_cai_dat.groupby(uid)
#
# data structure to store info ...
stat_dat = quantile_summary_storage_init(num_of_quantiles=the_num_of_quantiles)
# data structure to ORIGINAL store info ...
original_stat_dat = quantile_summary_storage_init(num_of_quantiles=the_num_of_quantiles)
# start iterating over different organisms ...
for idx,topt in env_dat[[uid,topt_id]].itertuples(index=False):
    if topt < 1000:
        # halophiles already excluded ...
        cds_cai_dat =  gen_dat_org.get_group(idx) 
        original_cds_cai_dat = original_gen_dat_org.get_group(idx) 
        # is it a translationally optimized organism ?
        # after messing up codons, 0 organisms are going to be TrOp..
        # so, just use the original definition to check if TrOp affects the results ...
        trop_status = get_one_trop(original_gen_dat_org,idx)
        #
        if trop_status != 'none':
            # fill in the record for each quantile (random-shuffled) ...
            stat_dat[uid].append(idx)
            stat_dat[topt_id].append(topt)
            stat_dat['TrOp'].append(trop_status)
            counter = 0
            for fivywrel_qs,r20_qs,cost_qs in zip(*get_quantiles_summary(cds_cai_dat,the_num_of_quantiles,thermo_freq[1],akashi_cost[1])):
                stat_dat['q%d'%counter].append(fivywrel_qs)
                stat_dat['R20_q%d'%counter].append(r20_qs)
                stat_dat['Akashi_q%d'%counter].append(cost_qs)
                # stat_dat['ProtLen_q%d'%counter].append(protlen_qs)
                counter += 1
            ########################################
            # fill in the record for each quantile (ORIGINAL) ...
            original_stat_dat[uid].append(idx)
            original_stat_dat[topt_id].append(topt)
            original_stat_dat['TrOp'].append(trop_status)
            counter = 0
            for fivywrel_qs,r20_qs,cost_qs in zip(*get_quantiles_summary(original_cds_cai_dat,the_num_of_quantiles,thermo_freq[1],akashi_cost[1])):
                original_stat_dat['q%d'%counter].append(fivywrel_qs)
                original_stat_dat['R20_q%d'%counter].append(r20_qs)
                original_stat_dat['Akashi_q%d'%counter].append(cost_qs)
                # original_stat_dat['ProtLen_q%d'%counter].append(protlen_qs)
                counter += 1
#
cai_stats_quant = pd.DataFrame(stat_dat)
cai_stats_quant_TrOp = cai_stats_quant[cai_stats_quant['TrOp']=='true']
cai_stats_quant_noTrOp = cai_stats_quant[cai_stats_quant['TrOp']=='false']
#
original_cai_stats_quant = pd.DataFrame(original_stat_dat)
original_cai_stats_quant_TrOp = original_cai_stats_quant[original_cai_stats_quant['TrOp']=='true']
original_cai_stats_quant_noTrOp = original_cai_stats_quant[original_cai_stats_quant['TrOp']=='false']
#


###############################################
# PLOTTING ...
###############################################
# k1, k2, k3 = 'q%d'%i,'R20_q%d'%i, 'Akashi_q%d'%i
def quantile_plotter(dat,kx,fname,ax=None,savefig=False,title='',color='blue',lims=False,x_offset=0):
    #############################
    def lims_to_range(local_lims,ext_coeff = 1.1):
        # small conversion func to turnd data range limits
        # into slightly extended limist for plotting ...
        mean = 0.5*sum(local_lims)
        ext_half = ext_coeff*abs(local_lims[0] - mean)
        return (mean - ext_half, mean + ext_half)
    ###############################
    if ax is None:
        plt.clf()
        ax = plt.gca()
    else:
        pass
    ###############################
    x_toplot, y_toplot, err_plot = [],[],[]
    for i in range(the_num_of_quantiles):
        kxi,ylabel = kx(i,label=True)
        x_toplot.append(i+1+x_offset)
        y_toplot.append(dat[kxi].mean())
        err_plot.append(dat[kxi].std())
    # plotting ...
    ax.errorbar(x_toplot,y_toplot,yerr=err_plot,fmt='o',color=color, label=title, mew=0,ms=6)
    ax.set_xlim(0.5,5.5)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('CAI quantile')
    #
    #########################
    a,b,r,pval,_ = st.linregress(x_toplot,y_toplot)
    # label = flabel(r,pval)
    x_toplot = np.asarray(x_toplot)
    f_label = lambda r,pval: "linear fit, R=%.2f "%r + (r"$p=%.3f$"%pval if pval>=0.001 else r"$p<0.001$")
    ax.plot(x_toplot,a*x_toplot+b,'-',color=color,lw=2,label=f_label(r,pval))
    ######################
    # # ...
    # if title:
    #     ax.set_title(title)
    if lims:
        ax.set_ylim( lims_to_range(lims) )
    # make room for legends ...
    ymin,ymax = ax.get_ylim()
    y_span = ymax-ymin
    ax.set_ylim((ymin-0.5*y_span,ymax+0.5*y_span))
    #
    # save figures on demand ...
    if savefig:
        ax.legend(loc='best',ncol=2)
        plt.savefig(fname)
    else:
        return ax


def update_lims(dat,fkey,old_lims=False):
    # calculate the new limits anyways ...
    lims_min = min( dat[fkey(i)].mean()-dat[fkey(i)].std() for i in range(the_num_of_quantiles) )
    lims_max = max( dat[fkey(i)].mean()+dat[fkey(i)].std() for i in range(the_num_of_quantiles) )
    new_lims = (lims_min,lims_max)
    # update these limits if necessary ...
    if not old_lims:
        return new_lims
    else:
        updated_lims_min = min(new_lims[0],old_lims[0])
        updated_lims_max = max(new_lims[1],old_lims[1])
        return ( updated_lims_min, updated_lims_max )







def IVYWREL_key(i,label=False):
    return ('q%d'%i,'IVYWREL') if label else 'q%d'%i

def R20_key(i,label=False):
    return ('R20_q%d'%i,'R20 self-exp_T') if label else 'R20_q%d'%i

def Akashi_key(i,label=False):
    return ('Akashi_q%d'%i,'Akashi cost') if label else 'Akashi_q%d'%i

# def ProtLen_key(i,label=False):
#     return ('ProtLen_q%d'%i,'mean protein length') if label else 'ProtLen_q%d'%i




# unified limits for the IVYWREL ...
# IVYWREL_lims = update_lims(cai_stats_quant,IVYWREL_key)
# IVYWREL_lims = update_lims(cai_stats_quant_noTrOp,IVYWREL_key,old_lims=IVYWREL_lims)
# IVYWREL_lims = update_lims(cai_stats_quant_TrOp,IVYWREL_key,old_lims=IVYWREL_lims)
IVYWREL_lims = update_lims(cai_stats_quant_TrOp,IVYWREL_key)
# IVYWREL_lims = update_lims(original_cai_stats_quant,IVYWREL_key,old_lims=IVYWREL_lims)
# IVYWREL_lims = update_lims(original_cai_stats_quant_noTrOp,IVYWREL_key,old_lims=IVYWREL_lims)
# IVYWREL_lims = update_lims(original_cai_stats_quant_TrOp,IVYWREL_key,old_lims=IVYWREL_lims)
IVYWREL_lims = update_lims(original_cai_stats_quant_TrOp,IVYWREL_key,old_lims=IVYWREL_lims)
# unified limits for the R20 ...
# R20_lims = update_lims(cai_stats_quant,R20_key)
# R20_lims = update_lims(cai_stats_quant_noTrOp,R20_key,old_lims=R20_lims)
# R20_lims = update_lims(cai_stats_quant_TrOp,R20_key,old_lims=R20_lims)
R20_lims = update_lims(cai_stats_quant_TrOp,R20_key)
# R20_lims = update_lims(original_cai_stats_quant,R20_key,old_lims=R20_lims)
# R20_lims = update_lims(original_cai_stats_quant_noTrOp,R20_key,old_lims=R20_lims)
# R20_lims = update_lims(original_cai_stats_quant_TrOp,R20_key,old_lims=R20_lims)
R20_lims = update_lims(original_cai_stats_quant_TrOp,R20_key,old_lims=R20_lims)
# unified limits for the Akashi ...
# Akashi_lims = update_lims(cai_stats_quant,Akashi_key)
# Akashi_lims = update_lims(cai_stats_quant_noTrOp,Akashi_key,old_lims=Akashi_lims)
# Akashi_lims = update_lims(cai_stats_quant_TrOp,Akashi_key,old_lims=Akashi_lims)
Akashi_lims = update_lims(cai_stats_quant_TrOp,Akashi_key)
# Akashi_lims = update_lims(original_cai_stats_quant,Akashi_key,old_lims=Akashi_lims)
# Akashi_lims = update_lims(original_cai_stats_quant_noTrOp,Akashi_key,old_lims=Akashi_lims)
# Akashi_lims = update_lims(original_cai_stats_quant_TrOp,Akashi_key,old_lims=Akashi_lims)
Akashi_lims = update_lims(original_cai_stats_quant_TrOp,Akashi_key,old_lims=Akashi_lims)
# #########################
# ProtLen_lims = update_lims(cai_stats_quant_TrOp,ProtLen_key)
# ProtLen_lims = update_lims(original_cai_stats_quant_TrOp,ProtLen_key,old_lims=ProtLen_lims)
#########################


def get_fname_title(keyid,kingdom,shuff_stat,trop_stat):
    # generate fname right away ...
    fname = "%s_%s_qunatile_trend.%s.%s.png"%(keyid,kingdom,shuff_stat,trop_stat)
    # figure out the readable title ...
    the_kingdom = 'Archaea' if kingdom=='arch' else 'Bacteria'
    the_shuff_stat = 'Shuffled' if shuff_stat=='shuff' else 'Original'
    the_trop_stat = 'CUS' if trop_stat=='trop' else ('non-CUS' if trop_stat=='notrop' else 'All')
    title = "%s %s (%s)"%(the_shuff_stat,the_kingdom,the_trop_stat)
    return (fname,title)


# TROP ONLY ...
fname,title = get_fname_title('IVYWREL','bact','shuff','trop')
ax_IVYWREL = quantile_plotter(cai_stats_quant_TrOp,IVYWREL_key,fname,ax=None,savefig=False,title=title,lims=IVYWREL_lims,x_offset=-0.05)
fname,title = get_fname_title('IVYWREL','bact','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,IVYWREL_key,fname,ax=ax_IVYWREL,savefig=True,title=title,color='red',lims=IVYWREL_lims,x_offset=0.05)


fname,title = get_fname_title('R20','bact','shuff','trop')
ax_R20 = quantile_plotter(cai_stats_quant_TrOp,R20_key,fname,ax=None,savefig=False,title=title,lims=R20_lims,x_offset=-0.05)
fname,title = get_fname_title('R20','bact','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,R20_key,fname,ax=ax_R20,savefig=True,title=title,color='red',lims=R20_lims,x_offset=0.05)


fname,title = get_fname_title('Akashi','bact','shuff','trop')
ax_Akashi = quantile_plotter(cai_stats_quant_TrOp,Akashi_key,fname,ax=None,savefig=False,title=title,lims=Akashi_lims,x_offset=-0.05)
fname,title = get_fname_title('Akashi','bact','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,Akashi_key,fname,ax=ax_Akashi,savefig=True,title=title,color='red',lims=Akashi_lims,x_offset=0.05)





def quantile_plotter_II(dat,kx,fname,ax=None,savefig=False,title='',color='blue',lims=False,x_offset=0):
    #############################
    def lims_to_range(local_lims,ext_coeff = 1.1):
        # small conversion func to turnd data range limits
        # into slightly extended limist for plotting ...
        mean = 0.5*sum(local_lims)
        ext_half = ext_coeff*abs(local_lims[0] - mean)
        return (mean - ext_half, mean + ext_half)
    ###############################
    if ax is None:
        plt.clf()
        ax = plt.gca()
    else:
        pass
    ###############################
    x_toplot, y_toplot, err_plot = [],[],[]
    for i in range(the_num_of_quantiles):
        kxi,ylabel = kx(i,label=True)
        x_toplot.append(i+1+x_offset)
        # y_toplot.append(dat[kxi].mean())
        y_toplot.append(dat[kxi].median())
        err_plot.append(dat[kxi].std())
    # plotting ...
    #
    # ax.errorbar(x_toplot,y_toplot,yerr=err_plot,fmt='o',color=color, label=title, mew=0,ms=6)
    rects1 = ax.bar(x_toplot, y_toplot, width=0.35, color=color, yerr=err_plot, label=title,error_kw=dict(elinewidth=2,ecolor='black'))
    # error_kw=dict(elinewidth=2,ecolor='red')
    #
    #
    ax.set_xlim(0.5,5.5)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('CAI quantile')
    #
    # #########################
    # a,b,r,pval,_ = st.linregress(x_toplot,y_toplot)
    # # label = flabel(r,pval)
    # x_toplot = np.asarray(x_toplot)
    # f_label = lambda r,pval: "linear fit, R=%.2f "%r + (r"$p=%.3f$"%pval if pval>=0.001 else r"$p<0.001$")
    # ax.plot(x_toplot,a*x_toplot+b,'-',color=color,lw=2,label=f_label(r,pval))
    # ######################
    # # ...
    # if title:
    #     ax.set_title(title)
    if lims:
        ax.set_ylim( lims_to_range(lims) )
    # make room for legends ...
    ymin,ymax = ax.get_ylim()
    y_span = ymax-ymin
    ax.set_ylim((ymin-0.2*y_span,ymax))
    #
    # save figures on demand ...
    if savefig:
        ax.legend(loc='lower left',ncol=2)
        # plt.savefig(fname)
        plt.show()
    else:
        return ax



fname,title = get_fname_title('R20','bact','shuff','trop')
ax_R20 = quantile_plotter(cai_stats_quant_TrOp,R20_key,fname,ax=None,savefig=False,title=title,lims=R20_lims,x_offset=-0.05)
fname,title = get_fname_title('R20','bact','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,R20_key,fname,ax=ax_R20,savefig=True,title=title,color='red',lims=R20_lims,x_offset=0.05)


from scipy.stats import norm

# original_cai_stats_quant_TrOp
ccc = ['R20_q0', 'R20_q1', 'R20_q2', 'R20_q3', 'R20_q4']
# ccc = ['Akashi_q0', 'Akashi_q1', 'Akashi_q2', 'Akashi_q3', 'Akashi_q4']
xmin = original_cai_stats_quant_TrOp[ccc].min().min()
xmax = original_cai_stats_quant_TrOp[ccc].max().max()
bins = np.linspace(xmin,xmax,50)
plt.clf()
fig,ax = plt.subplots(nrows=5,ncols=2,sharex=True,figsize=(7,7))
for i,c in enumerate(ccc):
    # Fit a normal distribution to the data:
    dat = cai_stats_quant_TrOp[cai_stats_quant_TrOp['OptimumTemperature']<50][c]
    mu, std = norm.fit(dat)
    p = norm.pdf(bins, mu, std)
    ax[len(ccc)-i-1,0].hist(dat.as_matrix(),bins=bins,color='blue',alpha=0.5,linewidth=0,log=False)
    ax[len(ccc)-i-1,0].plot(bins, p, color='blue', linewidth=2)
    #
    #
    #
    #
    dat = original_cai_stats_quant_TrOp[original_cai_stats_quant_TrOp['OptimumTemperature']<50][c]
    mu, std = norm.fit(dat)
    p = norm.pdf(bins, mu, std)
    ax[len(ccc)-i-1,1].hist(dat.as_matrix(),bins=bins,color='red',alpha=0.5,linewidth=0,log=False)
    # ax[len(ccc)-i-1].hist(cai_stats_quant_TrOp[c],bins=bins,color='blue',alpha=0.7)
    ax[len(ccc)-i-1,1].plot(bins, p, color='red', linewidth=2)
plt.show()






####################################################################################
#  EXTRA ANALYSIS FOR TROY'S LINEAR MODEL STUFF ...
####################################################################################
# gen_dat_org = gen_dat.groupby(uid)
# original_gen_dat_org = original_cai_dat.groupby(uid)
#
# def quantile_summary_storage_init(num_of_quantiles=5):
    # data structure to store info ...
def storage_init():
    the_storage = {uid:     [],
                'prot_id':  [],
                'prot_len': [],
                topt_id:    [],
                'CAI':      [],
                'TrOp':     []}
    for aa in aacids:
        the_storage[aa] = []
    # returning ...
    return the_storage
#
#
# data structure to store info ...
unroll_stat_dat = storage_init()
# data structure to ORIGINAL store info ...
unroll_original_stat_dat = storage_init()
# start iterating over different organisms ...
for idx,topt in env_dat[[uid,topt_id]].itertuples(index=False):
    # halophiles already excluded ...
    cds_cai_dat =  gen_dat_org.get_group(idx) 
    original_cds_cai_dat = original_gen_dat_org.get_group(idx) 
    # is it a translationally optimized organism ?
    # after messing up codons, 0 organisms are going to be TrOp..
    # so, just use the original definition to check if TrOp affects the results ...
    trop_status = get_one_trop(original_gen_dat_org,idx)
    #
    if trop_status == 'true':
        #
        # # fill in the record for each protein ...
        # [uid, 'pid', 'protein', 'CAI', 'TrOp']
        # [uid, 'pid', 'protein', 'CAI', 'TrOp']
        # cds_cai_dat
        #
        # SHUFFLED CODONS INFORMATION GATHERING ...
        for pid,prot,cai in cds_cai_dat[['pid', 'protein', 'CAI']].itertuples(index=False):
            ########################
            # some calculations ...
            prot_len = len(prot)
            aa_freq_20 = np.true_divide([prot.count(aa) for aa in aacids],float(prot_len))
            # print aa_freq_20
            ########################
            # storing the data ...
            unroll_stat_dat[uid].append(idx)
            unroll_stat_dat['prot_id'].append(pid)
            unroll_stat_dat['prot_len'].append(prot_len)
            unroll_stat_dat[topt_id].append(topt)
            unroll_stat_dat['CAI'].append(cai)
            unroll_stat_dat['TrOp'].append(trop_status)
            for aa_id, aa in enumerate(aacids):
                unroll_stat_dat[aa].append(aa_freq_20[aa_id])
        #
        # ORIGINAL GENOMES INFORMATION GATHERING ...
        for pid,prot,cai in original_cds_cai_dat[['pid', 'protein', 'CAI']].itertuples(index=False):
            ########################
            # some calculations ...
            prot_len = len(prot)
            aa_freq_20 = np.true_divide([prot.count(aa) for aa in aacids],float(prot_len))
            # print aa_freq_20
            ########################
            # storing the data ...
            unroll_original_stat_dat[uid].append(idx)
            unroll_original_stat_dat['prot_id'].append(pid)
            unroll_original_stat_dat['prot_len'].append(prot_len)
            unroll_original_stat_dat[topt_id].append(topt)
            unroll_original_stat_dat['CAI'].append(cai)
            unroll_original_stat_dat['TrOp'].append(trop_status)
            for aa_id, aa in enumerate(aacids):
                unroll_original_stat_dat[aa].append(aa_freq_20[aa_id])
        #
#
flat_stats_shuffled = pd.DataFrame(unroll_stat_dat)
flat_stats_original = pd.DataFrame(unroll_original_stat_dat)
#

####################################################################################
# [uid, 'prot_id', 'prot_len', topt_id, 'CAI'] + aacids




















