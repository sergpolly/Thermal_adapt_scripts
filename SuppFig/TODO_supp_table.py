import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
from functools import partial
import time

from multiprocessing import Pool

topt = 'OptimumTemperature'


# SPECIAL CASE: archaeal assemblies from GenBank@NCBI before any filtering 
# see second (2) paragraph in the Archaeal Wiki on GitHub:
# original file is saved as "assembly_summary_arch_ncbi.txt" and is the part of the Script Collection
dropbox_path = os.path.join( os.path.expanduser('~'),'Dropbox (UMASS MED - BIB)')
path_assembly_summary = os.path.join( dropbox_path, 'Thermal_adapt_scripts', 'ArchNew' )
arch_genbank_original = os.path.join(path_assembly_summary,'assembly_summary_arch_ncbi.txt')
asm_sum = pd.read_csv(arch_genbank_original,sep='\t')
original_number_archaea_assembly = asm_sum.shape[0]
print "# of archaeal assemblies in the GenBank DB at NCBI,", original_number_archaea_assembly

therm_ref = pd.read_csv("termoprop.txt",sep='\t')
therm_www = pd.read_csv("thermoprop.csv",header=None,names=['taxonid','name','ogt','tmin','tmax','trange','origin','href'],na_values=['\N','0','\"\"'])
therm_www_true = therm_www[therm_www['ogt'].notnull()]


# just check is 2 thermo files refers to the same set of archaea ...
print "The following 2 lines must be True, otherwise termoprop.txt and thermoprop.csv are reffering to different archaea ..."
_1 = therm_ref.name.isin(therm_www_true.name).all()
_2 = therm_www_true.name.isin(therm_ref.name).all()
print _1
print _2
if _1 and _2:
    print "Everything seems OK!"
else:
    print "check files content!"
    sys.exit(1)

print "Same for taxonid-s ..."
therm_ref.taxonid.isin(therm_www_true.taxonid).all()
therm_www_true.taxonid.isin(therm_ref.taxonid).all()
########################################################################################


#
path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
#
root_path = os.path.expanduser('~')
bact_path = os.path.join(root_path,'GENOMES_BACTER_RELEASE69/genbank')
arch_path = os.path.join(root_path,'GENOMES_ARCH_SEP2015')
# SOME ARCHAEAL DATA ...
arch        = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest.dat'))
arch_nohalo = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest_no_halop.dat'))
number_arch_topt_genome = arch.shape[0]
print "# of archaeal genomes with sufficient annotation and OGT, ",number_arch_topt_genome
number_arch_topt_genome_nohalo = arch_nohalo.shape[0]
print "# of archaeal genomes with sufficient annotation and OGT (Halophiles excluded), ",number_arch_topt_genome_nohalo
# complete_arch_CDS_CAI_DNA.dat is needed to analyse TrOp organisms etc ...



# SOME BACTERIAL DATA ...
# complete genomes only ...
bact        = pd.read_csv(os.path.join(bact_path,'env_catalog_compgenome.dat'))
number_bact_topt_compgenome = bact.shape[0]
print "# of bacterial genomes with full annotations and OGT (analyzed genomes), ",number_bact_topt_compgenome


# TOTAL 
# total # of genomes analyzed ...
print "Total # of genomes analysed: archaea without Halophiles, bacteria complete genomes only = %d"%\
											(number_arch_topt_genome_nohalo + number_bact_topt_compgenome)
print
###################################################################################################################
print
print
print "We are matching thermoprop.csv's href to arch dataframe using both taxonid-s and organism names in both tables:"
print "There are 2 taxonid-s available for arch df: taxid and species_taxid, and 1 taxonid on the thermoprop side."
print "There is 1-on-1 organism name correspondance between arch and thermoprop."
atax = arch.taxid.isin(therm_www_true.taxonid)
asptax = arch.species_taxid.isin(therm_www_true.taxonid)
aname = arch.organism_name.isin(therm_www_true.name)
print " Can we match thermoprop with arch dataframe? Answer: "+str((atax|asptax|aname).all())
print
#########################################################
print
print "Matching procedure is in action!"
matched_indices = []
for taxidx, sptaxidx, aname in arch[['taxid','species_taxid','organism_name']].itertuples(index=False):
    _1 = therm_www_true[ therm_www_true.taxonid == taxidx ].index.tolist()
    _2 = therm_www_true[ therm_www_true.taxonid == sptaxidx ].index.tolist()
    _3 = therm_www_true[ therm_www_true.name == aname ].index.tolist()
    #
    if ((len(_1)<=1)and(len(_2)<=1)and(len(_3)<=1))and(bool(_1) or bool(_2) or bool(_3)):
        idx1, = _1 if bool(_1) else [None,]
        idx2, = _2 if bool(_2) else [None,]
        idx3, = _3 if bool(_3) else [None,]
        # print idx1,idx2,idx3
    else:
        print "Achtung! %s %s %s arch guy isn't found!"%(str(taxidx),str(sptaxidx),aname)
        sys.exit(1)
    # small histogramming along the way 
    match_results = dict((idx,0) for idx in set([idx1,idx2,idx3]))
    for idx in [idx1,idx2,idx3]:
        match_results[idx] += 1 
    ############################################################################
    # some rather artificial and empirical choice mechanisms is to follow ...
    ############################################################################
    if None in match_results:
        # print idx1,idx2,idx3
        if len(match_results)==2:
            _,the_index = sorted(match_results)
            # print the_index,': ',idx1,idx2,idx3
        else:
            _,the_index,_ = sorted(match_results)
            # print the_index,': ',idx1,idx2,idx3
        matched_indices.append(the_index)
    else:
        # print idx1,idx2,idx3
        assert 1<=len(match_results)<=2
        the_index = sorted(match_results,key=lambda _: match_results[_],reverse=True)[0]
        # print the_index,':   ',idx1,idx2,idx3
        matched_indices.append(the_index)


#######################################
# SAME FOR ARCH_NOHALO ...
print
print "Do similar analysis for archaea without halophiles ..."
matched_indices_nohalo = []
for taxidx, sptaxidx, aname in arch_nohalo[['taxid','species_taxid','organism_name']].itertuples(index=False):
    _1 = therm_www_true[ therm_www_true.taxonid == taxidx ].index.tolist()
    _2 = therm_www_true[ therm_www_true.taxonid == sptaxidx ].index.tolist()
    _3 = therm_www_true[ therm_www_true.name == aname ].index.tolist()
    #
    if ((len(_1)<=1)and(len(_2)<=1)and(len(_3)<=1))and(bool(_1) or bool(_2) or bool(_3)):
        idx1, = _1 if bool(_1) else [None,]
        idx2, = _2 if bool(_2) else [None,]
        idx3, = _3 if bool(_3) else [None,]
        # print idx1,idx2,idx3
    else:
        print "Achtung! %s %s %s arch guy isn't found!"%(str(taxidx),str(sptaxidx),aname)
        sys.exit(1)
    # small histogramming along the way 
    match_results = dict((idx,0) for idx in set([idx1,idx2,idx3]))
    for idx in [idx1,idx2,idx3]:
        match_results[idx] += 1 
    ############################################################################
    # some rather artificial and empirical choice mechanisms is to follow ...
    ############################################################################
    if None in match_results:
        # print idx1,idx2,idx3
        if len(match_results)==2:
            _,the_index = sorted(match_results)
            # print the_index,': ',idx1,idx2,idx3
        else:
            _,the_index,_ = sorted(match_results)
            # print the_index,': ',idx1,idx2,idx3
        matched_indices_nohalo.append(the_index)
    else:
        # print idx1,idx2,idx3
        assert 1<=len(match_results)<=2
        the_index = sorted(match_results,key=lambda _: match_results[_],reverse=True)[0]
        # print the_index,':   ',idx1,idx2,idx3
        matched_indices_nohalo.append(the_index)


arch['href'] = therm_www_true['href'].loc[matched_indices].tolist() # this line correctly matched what we need ...
arch_nohalo['href'] = therm_www_true['href'].loc[matched_indices_nohalo].tolist() # this line correctly matched what we need ...



####################################################################################################################

# print
# print "Thermophiles ..."
# ############################
# # how many Thermophiles (OGT>50) ...
# number_arch_topt_genome_nohalo_THERM = arch_nohalo[arch_nohalo[topt] > 50].shape[0]
# print "# of archaeal Thermophiles with sufficient annotation and OGT, ",number_arch_topt_genome_nohalo_THERM
# ####################
# number_bact_topt_compgenome_THERM = bact[bact[topt] > 50].shape[0]
# print "# of bacterial Thermophiles with full annotations and OGT, ",number_bact_topt_compgenome_THERM
# ############################


###############################################
# complete_CDS_CAI_DNA.dat same thing ...

arch_cai_fname = os.path.join(arch_path,"complete_arch_CDS_CAI_DNA.dat")
bact_cai_fname = os.path.join(bact_path,"complete_CDS_CAI_DNA.dat")

arch_cai = pd.read_csv(arch_cai_fname)
bact_cai = pd.read_csv(bact_cai_fname)

bact_cai_org = bact_cai.groupby('GenomicID')
arch_cai_org = arch_cai.groupby('assembly_accession')


# # This is a bug potentially ...
# pd.Series([np.nan, np.nan, np.nan]).all() = True
# pd.Series([np.nan, np.nan, np.nan]).any() = False
def get_one_trop(idx, all_cai):
    org_cds = all_cai.get_group(idx)
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


# num_TrOp_arch_nohalo    = sum( get_one_trop(idx,arch_cai_org) for idx in arch_nohalo['assembly_accession'] )
# num_TrOp_arch           = sum( get_one_trop(idx,arch_cai_org) for idx in arch['assembly_accession'])
# num_TrOp_bact           = sum( get_one_trop(idx,bact_cai_org) for idx in bact['GenomicID'])
# print "# of Archaeal species with Translational Optimization, ",num_TrOp_arch
# print "# of Archaeal species with Translational Optimization (excluding Halophiles), ",num_TrOp_arch_nohalo
# print "# of Bacterial species with Translational Optimization, ",num_TrOp_bact

col_TrOp_arch_nohalo    = [get_one_trop(idx,arch_cai_org) for idx in arch_nohalo['assembly_accession'] ]
col_TrOp_arch           = [get_one_trop(idx,arch_cai_org) for idx in arch['assembly_accession'] ]
col_TrOp_bact           = [get_one_trop(idx,bact_cai_org) for idx in bact['GenomicID'] ]
##################################
arch['CUS']          = col_TrOp_arch 
arch_nohalo['CUS']   = col_TrOp_arch_nohalo
bact['CUS']          = col_TrOp_bact







####################################################################
####################################################################
#
#    PREPARING ARCHAEAL TABLE ...
#
####################################################################
####################################################################
arch['Halophilic'] = ~arch['assembly_accession'].isin(arch_nohalo['assembly_accession'])
# join taxon id-s to free some space ...
arch['taxonid'] = arch['species_taxid'].apply(str)+','+arch['taxid'].apply(str)

arch_cols = ['AssemblyID',
        'Organism name',
        'OGT',
        # 'OptimumTemperature',
        'taxonid',
        # 'bioproject',
        'CUS',
        'Halophilic',
        'FTP ref',
        'OGT link']
#########################################
f_rename = {'OptimumTemperature':'OGT','assembly_accession':'AssemblyID','organism_name':'Organism name','href':'OGT link','ftp_path':'FTP ref'}
# formatters
f_truth = lambda _:'Y'if _ else 'N'
f_identity = lambda _:_
# f_ext = lambda _: str(_[:100])
f_fomat = {'href':f_identity, #lambda _:"\\tiny{%s}"%_,
            'CUS':f_identity,
            'Halophilic':f_identity}
###########################
pd.set_option('max_colwidth',2000)
##############################
with open('arch_table.htm','w') as fp:
    # fp.write("\documentclass[8pt]{extreport}\n")
    # fp.write("\usepackage{booktabs}\n")
    # fp.write("\usepackage{longtable}\n")
    # fp.write("\usepackage[top=1in, bottom=1.25in, left=0.3in, right=1.25in]{geometry}\n")
    # fp.write("\\begin{document}\n")
    # fp.write("{\small \n")
    # to_latex output ...
    arch.rename(columns=f_rename).to_html(fp, columns=arch_cols, col_space=20,
                header=True, index=False, na_rep='NaN',
                formatters=f_fomat, float_format=None, sparsify=None,
                index_names=True, bold_rows=True,
                escape=True)
    # fp.write("} \n")
    # fp.write("\end{document}\n")











####################################################################
####################################################################
#
#    PREPARING BACTERIAL TABLE ...
#
####################################################################
####################################################################
# arch['Halophilic'] = ~arch['assembly_accession'].isin(arch_nohalo['assembly_accession'])
# join taxon id-s to free some space ...
# arch['taxonid'] = arch['species_taxid'].apply(str)+','+arch['taxid'].apply(str)


bact_cols = ['GenomicID',
    'Organism name',
    'OGT',
    'CUS',
    'BioProjectID',
    'TaxonID']
#########################################
f_rename = {'OptimumTemperature':'OGT','Organism_des':'Organism name','BioProject':'BioProjectID'}
# formatters
f_truth = lambda _:'Y'if _ else 'N'
f_identity = lambda _:_
# f_ext = lambda _: str(_[:100])
###########################
pd.set_option('max_colwidth',2000)
##############################
with open('bact_table.htm','w') as fp:
    # fp.write("\documentclass[8pt]{extreport}\n")
    # fp.write("\usepackage{booktabs}\n")
    # fp.write("\usepackage{longtable}\n")
    # fp.write("\usepackage[top=1in, bottom=1.25in, left=0.3in, right=1.25in]{geometry}\n")
    # fp.write("\\begin{document}\n")
    # fp.write("{\small \n")
    # to_latex output ...
    bact.rename(columns=f_rename).to_html(fp, columns=bact_cols, col_space=20,
                header=True, index=False, na_rep='NaN',
                formatters=None, float_format=None, sparsify=None,
                index_names=True, bold_rows=True,
                escape=True)
    # fp.write("} \n")
    # fp.write("\end{document}\n")













# \documentclass{report}
# \usepackage{booktabs}
# \usepackage{longtable}
# \usepackage[top=1in, bottom=1.25in, left=0.5in, right=1.25in]{geometry}
# \begin{document}
# \end{document}












































