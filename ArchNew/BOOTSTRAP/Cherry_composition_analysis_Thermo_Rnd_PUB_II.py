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


# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
root_path = os.path.expanduser('~')
bact_path = os.path.join(root_path,'GENOMES_BACTER_RELEASE69/genbank')
arch_path = os.path.join(root_path,'GENOMES_ARCH_SEP2015')

# this is Archaea script, so just do 
path = arch_path

uid = 'assembly_accession'
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
    fivywrel_cat, r20_cat, cost_cat, prot_len_cat = [],[],[],[]
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
        prot_len_cat.append(average_length)
    #returning ...
    return (fivywrel_cat,r20_cat,cost_cat,prot_len_cat)




def quantile_summary_storage_init(num_of_quantiles=5):
    # data structure to store info ...
    the_storage = {uid:[],
                topt_id:[],
                'TrOp':[]}
    for i in range(num_of_quantiles):
        the_storage['q%d'%i] = []
        the_storage['R20_q%d'%i] = []
        the_storage['Akashi_q%d'%i] = []
        the_storage['ProtLen_q%d'%i] = []
    # returning ...
    return the_storage






# ['DbxRefs','Description','FeaturesNum','assembly_accession','GenomicLen','GenomicName','Keywords','NucsPresent','Organism_des',
# 'SourceDbxRefs','SourceOrganism','SourcePlasmid','SourceStrain','Taxonomy','BioProject','TaxonID','Organism_env',
# 'OptimumTemperature','TemperatureRange','OxygenReq','Habitat','Salinity','crit_NC','crit_WGS','crit_genlen',
# 'crit_features','crit_comp_genome','crit_plasmid']
# env_dat = pd.read_csv(os.path.join(path,"summary_organisms_interest.dat"))
env_dat = pd.read_csv(os.path.join(path,'summary_organisms_interest_no_halop.dat'))
#['assembly_accession','cDNA','fid','pid','product','protein','status','table','ribosomal','CAI','TrOp']
gen_dat = pd.read_csv(os.path.join(path,"complete_arch_CDS_CAI_DNA_Rnd.dat"))
#
original_cai_dat = pd.read_csv(os.path.join(path,"complete_arch_CDS_CAI_DNA.dat"))
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
thermo_freq = pd.read_csv(os.path.join(path,'arch_thermo.dat'),header=None,sep=' ')
#
akashi_cost.set_index(0,inplace=True)
argentina_cost.set_index(0,inplace=True)
thermo_freq.set_index(0,inplace=True)
#
akashi_cost.sort_index(inplace=True)
argentina_cost.sort_index(inplace=True)
thermo_freq.sort_index(inplace=True)

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
    # halophiles already excluded ...
    cds_cai_dat = gen_dat_org.get_group(idx) 
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
        for fivywrel_qs,r20_qs,cost_qs,protlen_qs in zip(*get_quantiles_summary(cds_cai_dat,the_num_of_quantiles,thermo_freq[1],akashi_cost[1])):
            stat_dat['q%d'%counter].append(fivywrel_qs)
            stat_dat['R20_q%d'%counter].append(r20_qs)
            stat_dat['Akashi_q%d'%counter].append(cost_qs)
            stat_dat['ProtLen_q%d'%counter].append(protlen_qs)
            counter += 1
        ########################################
        # fill in the record for each quantile (ORIGINAL) ...
        original_stat_dat[uid].append(idx)
        original_stat_dat[topt_id].append(topt)
        original_stat_dat['TrOp'].append(trop_status)
        counter = 0
        for fivywrel_qs,r20_qs,cost_qs,protlen_qs in zip(*get_quantiles_summary(original_cds_cai_dat,the_num_of_quantiles,thermo_freq[1],akashi_cost[1])):
            original_stat_dat['q%d'%counter].append(fivywrel_qs)
            original_stat_dat['R20_q%d'%counter].append(r20_qs)
            original_stat_dat['Akashi_q%d'%counter].append(cost_qs)
            original_stat_dat['ProtLen_q%d'%counter].append(protlen_qs)
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
    ax.errorbar(x_toplot,y_toplot,yerr=err_plot,fmt='o',color=color, label=title)
    ax.set_xlim(0,6)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('CAI quantile')
    #
    #########################
    a,b,r,pval,_ = st.linregress(x_toplot,y_toplot)
    # label = flabel(r,pval)
    x_toplot = np.asarray(x_toplot)
    ax.plot(x_toplot,a*x_toplot+b,'-',color=color,lw=3,label="fit R=%.2f p=%.3f"%(r,pval))
    ######################
    # # ...
    # if title:
    #     ax.set_title(title)
    if lims:
        ax.set_ylim( lims_to_range(lims) )
    # save figures on demand ...
    if savefig:
        ax.legend(loc='best')
        plt.savefig(fname)
    else:
        return ax
######################################


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

def ProtLen_key(i,label=False):
    return ('ProtLen_q%d'%i,'mean protein length') if label else 'ProtLen_q%d'%i


# y limits for the plots ...
IVYWREL_lims = update_lims(cai_stats_quant_TrOp,IVYWREL_key)
IVYWREL_lims = update_lims(original_cai_stats_quant_TrOp,IVYWREL_key,old_lims=IVYWREL_lims)
#########################
R20_lims = update_lims(cai_stats_quant_TrOp,R20_key)
R20_lims = update_lims(original_cai_stats_quant_TrOp,R20_key,old_lims=R20_lims)
#########################
Akashi_lims = update_lims(cai_stats_quant_TrOp,Akashi_key)
Akashi_lims = update_lims(original_cai_stats_quant_TrOp,Akashi_key,old_lims=Akashi_lims)
#########################
ProtLen_lims = update_lims(cai_stats_quant_TrOp,ProtLen_key)
ProtLen_lims = update_lims(original_cai_stats_quant_TrOp,ProtLen_key,old_lims=ProtLen_lims)
#########################

def get_fname_title(keyid,kingdom,shuff_stat,trop_stat):
    # generate fname right away ...
    fname = "%s_%s_qunatile_trend.%s.%s.png"%(keyid,kingdom,shuff_stat,trop_stat)
    # figure out the readable title ...
    the_kingdom = 'Archaea' if kingdom=='arch' else 'Bacteria'
    the_shuff_stat = 'Shuffled' if shuff_stat=='shuff' else 'Original'
    the_trop_stat = 'Tr.Op.' if trop_stat=='trop' else ('non-Tr.Op.' if trop_stat=='notrop' else 'All')
    title = "%s %s (%s)"%(the_shuff_stat,the_kingdom,the_trop_stat)
    return (fname,title)


# # SHUFFLED ...
# fname,title = get_fname_title('IVYWREL','arch','shuff','notrop')
# quantile_plotter(cai_stats_quant_noTrOp,IVYWREL_key,fname,title,lims=IVYWREL_lims)
# fname,title = get_fname_title('R20','arch','shuff','notrop')
# quantile_plotter(cai_stats_quant_noTrOp,R20_key,fname,title,lims=R20_lims)
# fname,title = get_fname_title('Akashi','arch','shuff','notrop')
# quantile_plotter(cai_stats_quant_noTrOp,Akashi_key,fname,title,lims=Akashi_lims)
# #####################################################################################################
# fname,title = get_fname_title('IVYWREL','arch','shuff','all')
# quantile_plotter(cai_stats_quant,IVYWREL_key,fname,title,lims=IVYWREL_lims)
# fname,title = get_fname_title('R20','arch','shuff','all')
# quantile_plotter(cai_stats_quant,R20_key,fname,title,lims=R20_lims)
# fname,title = get_fname_title('Akashi','arch','shuff','all')
# quantile_plotter(cai_stats_quant,Akashi_key,fname,title,lims=Akashi_lims)
# #####################################################################################################
# fname,title = get_fname_title('IVYWREL','arch','shuff','trop')
# quantile_plotter(cai_stats_quant_TrOp,IVYWREL_key,fname,title,lims=IVYWREL_lims)
# fname,title = get_fname_title('R20','arch','shuff','trop')
# quantile_plotter(cai_stats_quant_TrOp,R20_key,fname,title,lims=R20_lims)
# fname,title = get_fname_title('Akashi','arch','shuff','trop')
# quantile_plotter(cai_stats_quant_TrOp,Akashi_key,fname,title,lims=Akashi_lims)
# #####################################################################################################


# # ORIGINAL ...
# fname,title = get_fname_title('IVYWREL','arch','original','notrop')
# quantile_plotter(original_cai_stats_quant_noTrOp,IVYWREL_key,fname,title,color='red',lims=IVYWREL_lims)
# fname,title = get_fname_title('R20','arch','original','notrop')
# quantile_plotter(original_cai_stats_quant_noTrOp,R20_key,fname,title,color='red',lims=R20_lims)
# fname,title = get_fname_title('Akashi','arch','original','notrop')
# quantile_plotter(original_cai_stats_quant_noTrOp,Akashi_key,fname,title,color='red',lims=Akashi_lims)
# #####################################################################################################
# fname,title = get_fname_title('IVYWREL','arch','original','all')
# quantile_plotter(original_cai_stats_quant,IVYWREL_key,fname,title,color='red',lims=IVYWREL_lims)
# fname,title = get_fname_title('R20','arch','original','all')
# quantile_plotter(original_cai_stats_quant,R20_key,fname,title,color='red',lims=R20_lims)
# fname,title = get_fname_title('Akashi','arch','original','all')
# quantile_plotter(original_cai_stats_quant,Akashi_key,fname,title,color='red',lims=Akashi_lims)
# #####################################################################################################
# fname,title = get_fname_title('IVYWREL','arch','original','trop')
# quantile_plotter(original_cai_stats_quant_TrOp,IVYWREL_key,fname,title,color='red',lims=IVYWREL_lims)
# fname,title = get_fname_title('R20','arch','original','trop')
# quantile_plotter(original_cai_stats_quant_TrOp,R20_key,fname,title,color='red',lims=R20_lims)
# fname,title = get_fname_title('Akashi','arch','original','trop')
# quantile_plotter(original_cai_stats_quant_TrOp,Akashi_key,fname,title,color='red',lims=Akashi_lims)
# #####################################################################################################




# TROP ONLY ...
fname,title = get_fname_title('IVYWREL','arch','shuff','trop')
ax_IVYWREL = quantile_plotter(cai_stats_quant_TrOp,IVYWREL_key,fname,ax=None,savefig=False,title=title,lims=IVYWREL_lims,x_offset=-0.05)
fname,title = get_fname_title('IVYWREL','arch','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,IVYWREL_key,fname,ax=ax_IVYWREL,savefig=True,title=title,color='red',lims=IVYWREL_lims,x_offset=0.05)


fname,title = get_fname_title('R20','arch','shuff','trop')
ax_R20 = quantile_plotter(cai_stats_quant_TrOp,R20_key,fname,ax=None,savefig=False,title=title,lims=R20_lims,x_offset=-0.05)
fname,title = get_fname_title('R20','arch','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,R20_key,fname,ax=ax_R20,savefig=True,title=title,color='red',lims=R20_lims,x_offset=0.05)


fname,title = get_fname_title('Akashi','arch','shuff','trop')
ax_Akashi = quantile_plotter(cai_stats_quant_TrOp,Akashi_key,fname,ax=None,savefig=False,title=title,lims=Akashi_lims,x_offset=-0.05)
fname,title = get_fname_title('Akashi','arch','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,Akashi_key,fname,ax=ax_Akashi,savefig=True,title=title,color='red',lims=Akashi_lims,x_offset=0.05)


fname,title = get_fname_title('ProtLen','arch','shuff','trop')
ax_ProtLen = quantile_plotter(cai_stats_quant_TrOp,ProtLen_key,fname,ax=None,savefig=False,title=title,lims=ProtLen_lims,x_offset=-0.05)
fname,title = get_fname_title('ProtLen','arch','original','trop')
quantile_plotter(original_cai_stats_quant_TrOp,ProtLen_key,fname,ax=ax_ProtLen,savefig=True,title=title,color='red',lims=ProtLen_lims,x_offset=0.05)






















