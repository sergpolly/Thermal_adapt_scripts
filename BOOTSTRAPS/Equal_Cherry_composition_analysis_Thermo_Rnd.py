import pandas as pd
import os
import subprocess as sub
import re
import sys
from Bio import SeqUtils
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')


# ['DbxRefs','Description','FeaturesNum','GenomicID','GenomicLen','GenomicName','Keywords','NucsPresent','Organism_des',
# 'SourceDbxRefs','SourceOrganism','SourcePlasmid','SourceStrain','Taxonomy','BioProject','TaxonID','Organism_env',
# 'OptimumTemperature','TemperatureRange','OxygenReq','Habitat','Salinity','crit_NC','crit_WGS','crit_genlen',
# 'crit_features','crit_comp_genome','crit_plasmid']
env_dat = pd.read_csv(os.path.join(path,"env_catalog_compgenome.dat"))

#['GenomicID','cDNA','fid','pid','product','protein','status','table','ribosomal','CAI','TrOp']
gen_dat = pd.read_csv(os.path.join(path,"complete_CDS_CAI_DNA_Rnd_Equal.dat"))


org_dat = pd.read_csv(os.path.join(path,"proteome_all.dat"))

# PROTEOME LEVEL AMINO ACID FREQUENCIES ...
# "proteome_all.dat"

# # file with the organisms of interest 
# dat_fname = os.path.join(bib2_scr_path,'catalog_with_accesion.dat')
# dat = pd.read_csv(dat_fname)
aacids = sorted(list('CMFILVWYAGTSNQDEHRKP'))


cost_vec_path = path
akashi = os.path.join(cost_vec_path,'akashi-cost.d')
argentina = os.path.join(cost_vec_path,'argentina-cost.d')

akashi_cost = pd.read_csv(akashi,header=None,sep=' ')
argentina_cost = pd.read_csv(argentina,header=None,sep=' ')
thermo_freq = pd.read_csv(os.path.join(path,'thermo.dat'),header=None,sep=' ')

akashi_cost.set_index(0,inplace=True)
argentina_cost.set_index(0,inplace=True)
thermo_freq.set_index(0,inplace=True)

akashi_cost.sort_index(inplace=True)
argentina_cost.sort_index(inplace=True)
thermo_freq.sort_index(inplace=True)


#
gen_dat_org = gen_dat.groupby('GenomicID')
# genom_id = orgs.groups.keys() # env_dat['GenomicID'] ...
# gen_dat_grouped.get_group(idx)
#
# how to get quantile ...
# q75 = pid_cai['CAI'].quantile(q=0.75)
#
#
num_of_quantiles = 5
#
stat_dat = {'GenomicID':[],
            'OptimumTemperature':[],
            'TrOp':[],
            'GC':[],
            'ProtLen':[]}
for i in range(num_of_quantiles):
    stat_dat['q%d'%i] = []
    stat_dat['R20_q%d'%i] = []
    stat_dat['Akashi_q%d'%i] = []
#
#
for idx,topt in env_dat[['GenomicID','OptimumTemperature']].itertuples(index=False):
    cds_cai_dat = gen_dat_org.get_group(idx) 
    # is it a translationally optimized organism ?
    all,any = cds_cai_dat['TrOp'].all(),cds_cai_dat['TrOp'].any()
    if all == any:
        trans_opt = all
    else:  #any != all
        print "%s@T=%f: Something wrong is happening: TrOp flag is not same for all ..."%(idx,topt)
    # THIS IS just a stupid precaution measure, in case we messed something upstream ...
    # not that stupid after all, because NaN is behaving badly here ...
    if cds_cai_dat['TrOp'].notnull().all():
        #
        # we can use this 'qcut' function from pandas to divide our proteins by the quantiles ...
        category,bins = pd.qcut(cds_cai_dat['CAI'],q=num_of_quantiles,retbins=True,labels=False)
        #
        stat_dat['GenomicID'].append(idx)
        stat_dat['OptimumTemperature'].append(topt)
        stat_dat['TrOp'].append(trans_opt)
        stat_dat['GC'].append(org_dat[org_dat.GenomicID==idx]['GC'])
        stat_dat['ProtLen'].append(org_dat[org_dat.GenomicID==idx]['ProtLen'])
        #
        # then we could iterate over proteins/cDNAs in these categories ...
        for cat in range(num_of_quantiles):
            cds_cai_category = cds_cai_dat[category==cat]
            total_length = cds_cai_category['protein'].str.len().sum()
            IVYWREL = sum(cds_cai_category['protein'].str.count(aa).sum() for aa in list('IVYWREL'))
            # IVYWREL = cds_cai_category['protein'].str.count('|'.join("IVYWREL")).sum() # tiny bit slower ...
            f_IVYWREL = float(IVYWREL)/float(total_length)
            # 20-vector for of amino acid composition ...
            aa_freq_20 = np.true_divide([cds_cai_category['protein'].str.count(aa).sum() for aa in aacids],float(total_length))
            # slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            _1,_2,R20,_4,_5 = stats.linregress(aa_freq_20, thermo_freq[1])
            # Akashi ...
            cost = np.dot(aa_freq_20,akashi_cost[1])
            # appending ...
            #
            #
            stat_dat['q%d'%cat].append(f_IVYWREL)
            stat_dat['R20_q%d'%cat].append(R20)
            stat_dat['Akashi_q%d'%cat].append(cost)
#
#
#
cai_stats_quant = pd.DataFrame(stat_dat)
#
cai_stats_quant_TrOp = cai_stats_quant[cai_stats_quant.TrOp]
cai_stats_quant_noTrOp = cai_stats_quant[~cai_stats_quant.TrOp]


plt.clf()
bins = np.linspace(-0.05,0.05,50)
# plt.hist(list(cai_stats_quant_TrOp.q4 - cai_stats_quant_TrOp.q1),bins=bins,color='blue')
plt.hist(list(cai_stats_quant.q4 - cai_stats_quant.q1),bins=bins,color='red',alpha=0.8)#,cumulative=True)
# plt.show()
plt.savefig("IVYWREL_quantile_hist_Eq.png")


plt.clf()
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.q1,'bo',alpha=0.8)
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.q4,'ro',alpha=0.8)
# plt.show()
plt.savefig("IVYWREL_dots_compare_Eq.png")
###############
 # plt.show()


plt.clf()
bins = np.linspace(-0.15,0.15,50)
# bins=50
plt.hist(list(cai_stats_quant[cai_stats_quant.OptimumTemperature<=50].R20_q4 - cai_stats_quant[cai_stats_quant.OptimumTemperature<=50].R20_q1),bins=bins,color='black',cumulative=False)
# plt.hist(list(cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature<=50].R20_q4 - cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature<=50].R20_q1),bins=bins,color='blue',cumulative=False)
# plt.hist(list(cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature<=50].R20_q4 - cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature<=50].R20_q1),bins=bins,color='red',alpha=0.8,cumulative=False)
plt.xlabel('$R^{4}_{T} - R^{1}_{T}$')
# plt.hist(list(cai_stats_quant_TrOp.R20_q4 - cai_stats_quant_TrOp.R20_q1),bins=bins,color='blue')
# plt.hist(list(cai_stats_quant_noTrOp.R20_q4 - cai_stats_quant_noTrOp.R20_q1),bins=bins,color='red',alpha=0.8)
plt.show()


# R20 ...
plt.clf()
delta_R20 = cai_stats_quant['R20_q4'] - cai_stats_quant['R20_q1']
plt.plot(cai_stats_quant[cai_stats_quant.OptimumTemperature>50].GC,delta_R20[cai_stats_quant.OptimumTemperature>50],'ro')
plt.plot(cai_stats_quant[cai_stats_quant.OptimumTemperature<=50].GC,delta_R20[cai_stats_quant.OptimumTemperature<=50],'bo')
plt.ylabel('$R^{4}_{T} - R^{1}_{T}$')
plt.xlabel('GC')
plt.show()

# IVYWREL
plt.clf()
delta_R20 = cai_stats_quant['q4'] - cai_stats_quant['q1']
plt.plot(cai_stats_quant[cai_stats_quant.OptimumTemperature>50].GC,delta_R20[cai_stats_quant.OptimumTemperature>50],'ro')
plt.plot(cai_stats_quant[cai_stats_quant.OptimumTemperature<=50].GC,delta_R20[cai_stats_quant.OptimumTemperature<=50],'bo')
plt.ylabel('$IVYWREL_{q4} - IVYWREL_{q1}$')
plt.xlabel('GC')
plt.show()

# Akashi 
plt.clf()
delta_R20 = cai_stats_quant['Akashi_q4'] - cai_stats_quant['Akashi_q1']
plt.plot(cai_stats_quant[cai_stats_quant.OptimumTemperature>50].GC,delta_R20[cai_stats_quant.OptimumTemperature>50],'ro')
plt.plot(cai_stats_quant[cai_stats_quant.OptimumTemperature<=50].GC,delta_R20[cai_stats_quant.OptimumTemperature<=50],'bo')
plt.ylabel('$Akashi_{q4} - Akashi_{q1}$')
plt.xlabel('GC')
plt.show()


# for the real ribosomal calculated CAI only ...
plt.clf()
delta_R20 = cai_stats_quant['R20_q4'] - cai_stats_quant['R20_q1']
plt.plot(cai_stats_quant[cai_stats_quant.TrOp].GC,delta_R20[cai_stats_quant.TrOp],'ro')
plt.plot(cai_stats_quant[~cai_stats_quant.TrOp].GC,delta_R20[~cai_stats_quant.TrOp],'bo')
plt.ylabel('$R^{4}_{T} - R^{1}_{T}$')
plt.xlabel('GC')
plt.show()
# for the real ribosomal calculated CAI only ...
plt.clf()
delta_R20 = cai_stats_quant['q4'] - cai_stats_quant['q1']
plt.plot(cai_stats_quant[cai_stats_quant.TrOp].GC,delta_R20[cai_stats_quant.TrOp],'ro')
plt.plot(cai_stats_quant[~cai_stats_quant.TrOp].GC,delta_R20[~cai_stats_quant.TrOp],'bo')
plt.ylabel('$IVYWREL_{q4} - IVYWREL_{q1}$')
plt.xlabel('GC')
plt.show()
# for the real ribosomal calculated CAI only ...
plt.clf()
delta_R20 = cai_stats_quant['Akashi_q4'] - cai_stats_quant['Akashi_q1']
plt.plot(cai_stats_quant[cai_stats_quant.TrOp].GC,delta_R20[cai_stats_quant.TrOp],'ro')
plt.plot(cai_stats_quant[~cai_stats_quant.TrOp].GC,delta_R20[~cai_stats_quant.TrOp],'bo')
plt.ylabel('$Akashi_{q4} - Akashi_{q1}$')
plt.xlabel('GC')
plt.show()




##################################
plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    kkk = k3
    #
    plt.errorbar([i+1,],cai_stats_quant[cai_stats_quant.OptimumTemperature>0][kkk].mean(),yerr=cai_stats_quant[cai_stats_quant.OptimumTemperature>0][kkk].std(),fmt='o')
plt.xlim(0,6)
plt.ylabel(kkk if kkk!=k1 else 'IVYWREL')
plt.xlabel('CAI quantile')
plt.show()
##################################
plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k2].mean(),yerr=cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k2].std(),fmt='o')
plt.xlim(0,6)
plt.show()




# R20 is flat on average (strange bi-modality?!) 
#       | meso thermo
# ------+-------------
# TrOp  |  NA  NA
# noTrOp|  ~~+   ~~-


# Akashi is flat on average (strange local minimum at middle CAI quantile)
#       | meso thermo
# ------+-------------
# TrOp  |  NA   NA
# noTrOp|  ~   ~


# IVYWREL is declining on average (?!)
#       | meso thermo
# ------+-------------
# TrOp  |  NA   NA
# noTrOp|  --   --








































