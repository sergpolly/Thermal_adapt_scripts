import pandas as pd
import numpy as np
import re
from fuzzywuzzy import fuzz
from fuzzywuzzy import process



# archaea genomes description from /ftp...ncbi/genomes/genbank/archaea/...
asm = pd.read_csv("assembly_summary_arch_ncbi.txt",sep='\t')
# geomes description and connection with organism names/assembly id/taxonid established earlier ...
phyl = pd.read_csv("phylog_temps.txt")
# KZ database of manually retrieved OptimumTemperatures of different archaea ...
therm = pd.read_csv("termoprop.txt",sep='\t')

# the problem is that the taxonid@therm is not fully within the taxonid@asm (~150 out of ~250).

# let's compare by Assembly accession index ...

# we need a transformation:
# GCA_000144915.1_ASM14491v1 -> GCA_000007185.1
# because left one is what we have in phyl, while right one is the assembly accession (from asm) ...
asm_tmp = re.compile("(.+\_.+\.\d+)\_.+")
asm_getid = lambda line: asm_tmp.match(line).groups()[0]
asm_if_check = lambda line: bool(asm_tmp.match(line))

# check if all asm index in phyl comply ...
print " Check if all comply ~GCA_xxxxxxxx.1_ASM14491v1: ", phyl.assembly.apply(asm_if_check).all()

# now get the accession form ...
asm_acc = phyl.assembly.apply(asm_getid)
phyl['assembly_accession'] = asm_acc

# and now check the intersection between asm_acc and accession from asm ...
asm_interection = np.intersect1d( asm_acc, asm.assembly_accession)

print "%d  out of %d from old list (phyl) are in the new GenBank list (asm)"%(asm_interection.size,asm_acc.shape[0])

dat = pd.merge(phyl, asm, how='inner', on='assembly_accession', suffixes=('_x', '_y'), copy=True)

# that's the end f verification by now.
# we would need to check organism names then and afterwards proceed with the analysis ...

# cols = ['taxonid', 'topt', 'name', 'fname', 'assembly', 'assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'isolate', 'version_status', 'assembly_level', 'release_type', 'genome_rep', 'asm_name', 'ftp_path']
cols = ['taxonid', 'topt', 'name', 'assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'version_status', 'assembly_level', 'genome_rep', 'assembly', 'asm_name']

print "Combined data output ..."
print dat[cols].head(3)

# asm has 2 kinds of taxonomy info: taxid and species_taxid
# In the phyl we have just one kind: taxonid
# let's what is the intersection here:

print dat[~((dat.taxonid == dat.taxid) | (dat.taxonid == dat.species_taxid))][cols]

# # organism_name               infraspecific_name
# asm['strain_opts'] = asm.infraspecific_name.apply(lambda line: tuple(strain.strip() for strain in str(line).strip().strip('strain=').split(';')) )
#
# # full_names_stack = []
# # for x,y in asm[['organism_name','strain_opts']].itertuples(index=False):
# #     full_names_stack.append(x)
# #     for strain in y:
# #         full_names_stack.append(' '.join((x,strain)))
#
# # for name in therm.name:
# #     check_simple = name in asm.organism_name


# EUREKA MOMENT, TRY EXTRACTING(matching)
# Organisms matched in assembly database either by name, taxonid, or species taxonid ...
asm_tax_matched = therm.name.isin(asm.organism_name)| therm.taxonid.isin(asm.taxid)| therm.taxonid.isin(asm.species_taxid)

# now we have to check if there are any collapsing matches ...
# asm['tid'] = asm[['organism_name','taxid','species_taxid']].itertuples(index=False)
tmp_list = []
variants = []
for name,tax,check_name,check_tax,check_sptax in zip(therm.name,therm.taxonid,therm.name.isin(asm.organism_name),therm.taxonid.isin(asm.taxid),therm.taxonid.isin(asm.species_taxid)):
    #
    if check_name and check_tax and check_sptax:
        query_res = asm[(asm.organism_name==name) & (asm.taxid==tax) & (asm.species_taxid==tax)]
        variants.append(1)
        tmp_list.append(query_res)
    elif check_name and check_tax:
        query_res = asm[(asm.organism_name==name) & (asm.taxid==tax)]
        variants.append(2)
        tmp_list.append(query_res)
    elif check_tax and check_sptax:
        query_res = asm[(asm.taxid==tax) & (asm.species_taxid==tax)]
        variants.append(3)
        tmp_list.append(query_res)
    elif check_name and check_sptax:
        query_res = asm[(asm.organism_name==name) & (asm.species_taxid==tax)]
        variants.append(4)
        tmp_list.append(query_res)
    elif check_name:
        query_res = asm[asm.organism_name==name]
        variants.append(5)
        tmp_list.append(query_res)
    elif check_tax:
        query_res = asm[asm.taxid==tax]
        variants.append(6)
        tmp_list.append(query_res)
    elif check_sptax:
        query_res = asm[asm.species_taxid==tax]
        variants.append(7)
        tmp_list.append(query_res)
    else:
        pass


# filter those with refseq_category=='representative genome' and  assembly_level=='Complete Genome' or 'chromosome' ...
both = []
rep = []
gen = []
either = []
justor = []
for asm_group in tmp_list:
    genome_stat = asm_group.assembly_level.isin(['Complete Genome','Chromosome',]).any()
    rep_stat = asm_group.refseq_category.isin(['representative genome',]).any()
    #
    #
    if genome_stat and rep_stat:
        both.append(1)
    if rep_stat:
        rep.append(1)
    if genome_stat:
        gen.append(1)
    if (genome_stat != rep_stat):
        either.append(1)
    if (genome_stat or rep_stat):
        justor.append(1)
    #
    #
    #
    #

print "both",len(both)
print "rep",len(rep)
print "gen",len(gen)
print "either",len(either)
print "justor",len(justor)


















