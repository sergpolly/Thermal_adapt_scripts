import pandas as pd
import os
import subprocess as sub
import re
import sys
from Bio import SeqUtils
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

path = "."


# ['DbxRefs','Description','FeaturesNum','GenomicID','GenomicLen','GenomicName','Keywords','NucsPresent','Organism_des',
# 'SourceDbxRefs','SourceOrganism','SourcePlasmid','SourceStrain','Taxonomy','BioProject','TaxonID','Organism_env',
# 'OptimumTemperature','TemperatureRange','OxygenReq','Habitat','Salinity','crit_NC','crit_WGS','crit_genlen',
# 'crit_features','crit_comp_genome','crit_plasmid']
env_dat = pd.read_csv(os.path.join(path,"summary_organisms_interest.dat"))
# summary_organisms_interest.dat

taxon_dat = pd.read_csv(os.path.join(path,"arch_taxonomy_interest.dat"))
check_halo = lambda tax_class: any(_ in tax_class for _ in ('Halobacteria','Nanohaloarchaea'))
taxon_dat['halo'] = taxon_dat['tax_lineages'].apply(lambda lins: any( check_halo(lin.split(';')) for lin in lins.split(':') ) )

#['GenomicID','cDNA','fid','pid','product','protein','status','table','ribosomal','CAI','TrOp']
gen_dat = pd.read_csv(os.path.join(path,"complete_arch_CDS_CAI_DNA.dat"))

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


#################################################
# after we processed all organism's genomes, we would need to 
# pull out proteins with top XXX percent CAI, to analyse their amino acid compositions ...
#################################################################
#################################################################
# (1) for each asm - open protein detail file, sort by CAI and analyse proteins with top 10% CAI ...
# (2) output analysis results to external file ...
#################################################################
#################################################################
# def get_aausage_proteome(seqrec):
#     # seqrec = db[seqrec_id]
#     features = seqrec.features
#     proteome = []
#     for feature in features:
#         qualifiers = feature.qualifiers
#         if (feature.type == 'CDS')and('translation' in qualifiers):
#             proteome.append(qualifiers['translation'][0])
#     #return the results ...
#     proteome = ''.join(proteome)
#     prot_len = float(len(proteome))
#     aa_freq = tuple(proteome.count(aa)/prot_len for aa in aacids)
#     #
#     return (int(prot_len),) + aa_freq

# def analyse_genome(db,seqrec_id):
#     seqrec = db[seqrec_id]
#     pl_aa_freq = get_aausage_proteome(seqrec)
#     gc = SeqUtils.GC(seqrec.seq)
#     id = seqrec.id
    # return (id,gc) + pl_aa_freq
# PERCENTILE = 0.1


# accounted_GC = []
# aafs = {}
# for aa in aacids:
#     aafs[aa] = []

# genome_length = []
# proteome_length = []
# and for each assembley it goes ...
######################
# fname = os.path.join(path_CAI,'%s_genes.dat'%asm)
######################
#
#
gen_dat_org = gen_dat.groupby('assembly_accession')
# genom_id = orgs.groups.keys() # env_dat['GenomicID'] ...
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
env_dat_tax = pd.merge(env_dat,taxon_dat,on='assembly_accession')
#
for idx,topt,halo in env_dat_tax[['assembly_accession','OptimumTemperature','halo']].itertuples(index=False):
    # excluding halophiles ...
    if not halo:
        cds_cai_dat = gen_dat_org.get_group(idx) 
        # is it a translationally optimized organism ?
        all_org,any_org = cds_cai_dat['TrOp'].all(),cds_cai_dat['TrOp'].any()
        if all_org == any_org:
    	    trans_opt = all_org
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
plt.hist(list(cai_stats_quant_TrOp.q4 - cai_stats_quant_TrOp.q1),bins=bins,color='blue')
plt.hist(list(cai_stats_quant_noTrOp.q4 - cai_stats_quant_noTrOp.q1),bins=bins,color='red',alpha=0.8)
plt.xlabel("IVYWREL(TrOp:blue,noTrOp:red)")
# plt.show()
plt.savefig("IVYWREL_arch_hist_TrueCodon.png")


plt.clf()
bins = np.linspace(-0.15,0.15,50)
# bins=50
plt.hist(list(cai_stats_quant[cai_stats_quant.OptimumTemperature<=50].R20_q4 - cai_stats_quant[cai_stats_quant.OptimumTemperature<=50].R20_q1),bins=bins,color='black',cumulative=False)
plt.hist(list(cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature<=50].R20_q4 - cai_stats_quant_noTrOp[cai_stats_quant_noTrOp.OptimumTemperature<=50].R20_q1),bins=bins,color='blue',cumulative=False)
plt.hist(list(cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature<=50].R20_q4 - cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature<=50].R20_q1),bins=bins,color='red',alpha=0.8,cumulative=False)
plt.xlabel('$R^{4}_{T} - R^{1}_{T}$')
plt.ylabel("hist AllMeso:black; noTrOpMeso:blue; TrOpMeso:red")
# plt.hist(list(cai_stats_quant_TrOp.R20_q4 - cai_stats_quant_TrOp.R20_q1),bins=bins,color='blue')
# plt.hist(list(cai_stats_quant_noTrOp.R20_q4 - cai_stats_quant_noTrOp.R20_q1),bins=bins,color='red',alpha=0.8)
plt.savefig("R20_arch_hist_TrueCodon.png")



plt.clf()
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.q1,'bo',alpha=0.8)
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.q4,'ro',alpha=0.8)
plt.xlabel("Temperature")
plt.ylabel("IVYWREL(HE:red;LE:blue)")
# plt.show()
plt.savefig("IVYWREL_T_dots_TrueCodon.png")
# #

plt.clf()
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.R20_q1,'bo',alpha=0.8)
plt.plot(cai_stats_quant.OptimumTemperature,cai_stats_quant.R20_q4,'ro',alpha=0.8)
plt.xlabel("Temperature")
plt.ylabel("R20(HE:red;LE:blue)")
plt.savefig("R20_T_dots_TrueCodon.png")

############################
# GC is lacking from the file env_dat is reffering to  ...
# FIX IT ...
###########################
# plt.clf()
# plt.plot(cai_stats_quant.GC,cai_stats_quant.R20_q1,'bo',alpha=0.8)
# plt.plot(cai_stats_quant.GC,cai_stats_quant.R20_q4,'ro',alpha=0.8)
# plt.show()


# plt.clf()
# for i in range(num_of_quantiles):
#     k1 = 'q%d'%i
#     k2 = 'R20_q%d'%i
#     k3 = 'Akashi_q%d'%i
#     #
#     plt.plot([i+1,]*cai_stats_quant.shape[0],cai_stats_quant[k1],alpha=0.7)

# plt.xlim(0,6)
# plt.show()



plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature>0][k1].mean(),yerr=cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature>0][k1].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k1)
plt.xlabel('CAI quantile')
plt.savefig("IVYWREL_arch_qunatile_trend_TrueCodon_TrOp.png")
# plt.show()


plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature>0][k2].mean(),yerr=cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature>0][k2].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k2)
plt.xlabel('CAI quantile')
plt.savefig("R20_arch_qunatile_trend_TrueCodon_TrOp.png")
# plt.show()


plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature>0][k3].mean(),yerr=cai_stats_quant_TrOp[cai_stats_quant_TrOp.OptimumTemperature>0][k3].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k3)
plt.xlabel('CAI quantile')
plt.savefig("Akashi_arch_qunatile_trend_TrueCodon_TrOp.png")
# plt.show()



######################################################################################


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
plt.savefig("IVYWREL_arch_qunatile_trend_TrueCodon_TrOp.png")
# plt.show()


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
plt.savefig("R20_arch_qunatile_trend_TrueCodon_noTrOp.png")
# plt.show()


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
plt.savefig("Akashi_arch_qunatile_trend_TrueCodon_noTrOp.png")
# plt.show()


########################################

plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant[cai_stats_quant.OptimumTemperature>0][k1].mean(),yerr=cai_stats_quant[cai_stats_quant.OptimumTemperature>0][k1].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k1)
plt.xlabel('CAI quantile')
plt.savefig("IVYWREL_arch_qunatile_trend_TrueCodon_ALL.png")
# plt.show()


plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant[cai_stats_quant.OptimumTemperature>0][k2].mean(),yerr=cai_stats_quant[cai_stats_quant.OptimumTemperature>0][k2].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k2)
plt.xlabel('CAI quantile')
plt.savefig("R20_arch_qunatile_trend_TrueCodon_ALL.png")
# plt.show()


plt.clf()
for i in range(num_of_quantiles):
    k1 = 'q%d'%i
    k2 = 'R20_q%d'%i
    k3 = 'Akashi_q%d'%i
    #
    plt.errorbar([i+1,],cai_stats_quant[cai_stats_quant.OptimumTemperature>0][k3].mean(),yerr=cai_stats_quant[cai_stats_quant.OptimumTemperature>0][k3].std(),fmt='o')

plt.xlim(0,6)
plt.ylabel(k3)
plt.xlabel('CAI quantile')
plt.savefig("Akashi_arch_qunatile_trend_TrueCodon_ALL.png")


# R20 grows on average, 
#       | meso thermo
# ------+-------------
# TrOp  |  ++  ~+
# noTrOp|  +   ~


# Akashi is declining on average 
#       | meso thermo
# ------+-------------
# TrOp  |  --   --
# noTrOp|  ~-   ~-


# IVYWREL is declining on average 
#       | meso thermo
# ------+-------------
# TrOp  |  --   ~-
# noTrOp|  ~-    -




# After reshuffling given CAI, everything becomes flat and 0-centered with much narrower distributions 






# #
# #     # move on ...
# #     # quantiles calculation ...
# #     q20,q40,q60,q80 = cds_cai_dat['CAI'].quantile(q=[0.2,0.4,0.6,0.8])
# #     #
# # q1_idx =                          (cds_cai_dat['CAI']<=q20)
# # q2_idx = (q20<cds_cai_dat['CAI'])&(cds_cai_dat['CAI']<=q40)
# # q3_idx = (q40<cds_cai_dat['CAI'])&(cds_cai_dat['CAI']<=q60)
# # q4_idx = (q60<cds_cai_dat['CAI'])&(cds_cai_dat['CAI']<=q80)
# # q5_idx = (q80<cds_cai_dat['CAI'])
# # # q3_idx = q40<cds_cai_dat['CAI']<=q60
# # # q4_idx = q60<cds_cai_dat['CAI']<=q80
# # # q5_idx = q80<cds_cai_dat['CAI']


# # ['q1', 'q2', 'q3', 'q4', 'q5']



# for IndexId,prot_det in dat[['IndexId','protein_details']].get_values():
#     ###################
#     # in the case of Bacteria, we know for sure that a single accession number refers to related
#     if prot_det:
#         # open a file with the analysed organismal proteins...
#         protein_fname = os.path.join(path_CAI,'%s_genes.dat'%IndexId)
#         # load the data ...
#         protein_dat = pd.read_csv(protein_fname)
#         # "prot_id,cai,gene_product,gene_seq,prot_seq" are the columns ...
#         # we'll be taking proteins with top PERCENTILE% CAI in the list ...
#         # protein total number is the first dimension of the table here ...
#         number_of_proteins,_ = protein_dat.shape
#         accounted_proteins = int(number_of_proteins*PERCENTILE)
#         # top PERCENTILE proteins will be considered for analysis ...
#         accounted_data = protein_dat.sort(columns='cai',ascending=False)[:accounted_proteins]
#         # analyse that stuff ...
#         cai_proteome = ''.join(accounted_data['prot_seq'])
#         cai_proteome_len = float(len(cai_proteome))
#         cai_genome = ''.join(accounted_data['gene_seq'])
#         #
#         accounted_GC.append(SeqUtils.GC(cai_genome))
#         #
#         for aa in aacids:
#             aafs[aa].append(cai_proteome.count(aa)/cai_proteome_len)
#     else:
#         # no protein details exist, no corresponding file exist at all...
#         accounted_GC.append(0.0)
#         for aa in aacids:
#             aafs[aa].append(0.0)




# dat['GC'] = accounted_GC
# for aa in aacids:
#     dat[aa] = aafs[aa]
# dat.to_csv('cai5_bacter.dat',index=False)





















































