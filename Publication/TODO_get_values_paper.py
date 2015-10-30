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




###############################################
# complete_CDS_CAI_DNA.dat same thing ...

arch_cai_fname = os.path.join(arch_path,"complete_arch_CDS_CAI_DNA.dat")
bact_cai_fname = os.path.join(bact_path,"complete_CDS_CAI_DNA.dat")

arch_cai = pd.read_csv(arch_cai_fname)
bact_cai = pd.read_csv(bact_cai_fname)

bact_cai_org = bact_cai.groupby('GenomicID')
arch_cai_org = arch_cai.groupby('assembly_accession')


# TODO ......................................................
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











