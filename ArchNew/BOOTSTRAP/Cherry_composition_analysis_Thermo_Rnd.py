import pandas as pd
import os
import subprocess as sub
import re
import sys
from Bio import SeqUtils
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
path = "."


# ['DbxRefs','Description','FeaturesNum','assembly_accession','GenomicLen','GenomicName','Keywords','NucsPresent','Organism_des',
# 'SourceDbxRefs','SourceOrganism','SourcePlasmid','SourceStrain','Taxonomy','BioProject','TaxonID','Organism_env',
# 'OptimumTemperature','TemperatureRange','OxygenReq','Habitat','Salinity','crit_NC','crit_WGS','crit_genlen',
# 'crit_features','crit_comp_genome','crit_plasmid']
env_dat = pd.read_csv(os.path.join(path,"summary_organisms_interest.dat"))

#['assembly_accession','cDNA','fid','pid','product','protein','status','table','ribosomal','CAI','TrOp']
gen_dat = pd.read_csv(os.path.join(path,"complete_arch_CDS_CAI_DNA_Rnd.dat"))

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
thermo_freq = pd.read_csv(os.path.join(path,'arch_thermo.dat'),header=None,sep=' ')

akashi_cost.set_index(0,inplace=True)
argentina_cost.set_index(0,inplace=True)
thermo_freq.set_index(0,inplace=True)

akashi_cost.sort_index(inplace=True)
argentina_cost.sort_index(inplace=True)
thermo_freq.sort_index(inplace=True)


#
gen_dat_org = gen_dat.groupby('assembly_accession')
# genom_id = orgs.groups.keys() # env_dat['assembly_accession'] ...
# gen_dat_grouped.get_group(idx)
#
# how to get quantile ...
# q75 = pid_cai['CAI'].quantile(q=0.75)
#
#
num_of_quantiles = 5
#
stat_dat = {'assembly_accession':[],
            'OptimumTemperature':[],
            'TrOp':[]}
for i in range(num_of_quantiles):
    stat_dat['q%d'%i] = []
    stat_dat['R20_q%d'%i] = []
    stat_dat['Akashi_q%d'%i] = []
#
#
for idx,topt in env_dat[['assembly_accession','OptimumTemperature']].itertuples(index=False):
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
        stat_dat['assembly_accession'].append(idx)
        stat_dat['OptimumTemperature'].append(topt)
        stat_dat['TrOp'].append(trans_opt)
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
plt.xlabel("IVYWREL(HExp)-IVYWREL(LExp)")
# plt.show()
plt.savefig("IVYWREL_quantile_hist_arch.png")


plt.clf()
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.q1,'bo',alpha=0.8)
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.q4,'ro',alpha=0.8)
plt.xlabel('Temperature')
plt.ylabel('IVYWREL(HE:red;LE:blue)')
# plt.show()
plt.savefig("IVYWREL_dots_compare_arch.png")



plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k1].mean(),yerr=cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k1].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k1)
plt.xlabel('CAI quantile')
plt.savefig("IVYWREL_arch_qunatile_trend_Shuff.png")

plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k2].mean(),yerr=cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k2].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k2)
plt.xlabel('CAI quantile')
plt.savefig("R20_arch_qunatile_trend_Shuff.png")

plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k3].mean(),yerr=cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature>0][k3].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k3)
plt.xlabel('CAI quantile')
plt.savefig("Akashi_arch_qunatile_trend_Shuff.png")

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








































