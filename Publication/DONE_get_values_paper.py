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
print "Thermophiles ..."
############################
# how many Thermophiles (OGT>50) ...
number_arch_topt_genome_nohalo_THERM = arch_nohalo[arch_nohalo[topt] > 50].shape[0]
print "# of archaeal Thermophiles with sufficient annotation and OGT, ",number_arch_topt_genome_nohalo_THERM
####################
number_bact_topt_compgenome_THERM = bact[bact[topt] > 50].shape[0]
print "# of bacterial Thermophiles with full annotations and OGT, ",number_bact_topt_compgenome_THERM
############################


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


num_TrOp_arch_nohalo    = sum( get_one_trop(idx,arch_cai_org) for idx in arch_nohalo['assembly_accession'] )
num_TrOp_arch           = sum( get_one_trop(idx,arch_cai_org) for idx in arch['assembly_accession'])
num_TrOp_bact           = sum( get_one_trop(idx,bact_cai_org) for idx in bact['GenomicID'])


print "# of Archaeal species with Translational Optimization, ",num_TrOp_arch
print "# of Archaeal species with Translational Optimization (excluding Halophiles), ",num_TrOp_arch_nohalo
print "# of Bacterial species with Translational Optimization, ",num_TrOp_bact








