import pandas as pd
import numpy as np
import re
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import copy

# archaea genomes description from /ftp...ncbi/genomes/genbank/archaea/...
asm = pd.read_csv("assembly_summary_arch_ncbi.txt",sep='\t')
# # geomes description and connection with organism names/assembly id/taxonid established earlier ...
# phyl = pd.read_csv("phylog_temps.txt")
# KZ database of manually retrieved OptimumTemperatures of different archaea ...
therm = pd.read_csv("termoprop.txt",sep='\t')


# Organisms matched in assembly database either by name, taxonid, or species taxonid ...
asm_tax_matched = therm.name.isin(asm.organism_name)| therm.taxonid.isin(asm.taxid)| therm.taxonid.isin(asm.species_taxid)
print "Number of organisms that were match in the current GenBank release either by name, taxid, or species_taxid: %d out of %d"%(asm_tax_matched.sum(),therm.shape[0])

# now we have to check if there are any collapsing matches ...
# asm['tid'] = asm[['organism_name','taxid','species_taxid']].itertuples(index=False)
matched_asm = []
variants = []
topt_aligned = []
for name,tax,topt,check_name,check_tax,check_sptax in zip(therm.name,therm.taxonid,therm.topt,therm.name.isin(asm.organism_name),therm.taxonid.isin(asm.taxid),therm.taxonid.isin(asm.species_taxid)):
    #
    if check_name and check_tax and check_sptax:
        query_res = asm[(asm.organism_name==name) & (asm.taxid==tax) & (asm.species_taxid==tax)]
        variants.append(1)
        matched_asm.append(query_res)
        topt_aligned.append(topt)
    elif check_name and check_tax:
        query_res = asm[(asm.organism_name==name) & (asm.taxid==tax)]
        variants.append(2)
        matched_asm.append(query_res)
        topt_aligned.append(topt)
    elif check_tax and check_sptax:
        query_res = asm[(asm.taxid==tax) & (asm.species_taxid==tax)]
        variants.append(3)
        matched_asm.append(query_res)
        topt_aligned.append(topt)
    elif check_name and check_sptax:
        query_res = asm[(asm.organism_name==name) & (asm.species_taxid==tax)]
        variants.append(4)
        matched_asm.append(query_res)
        topt_aligned.append(topt)
    elif check_name:
        query_res = asm[asm.organism_name==name]
        variants.append(5)
        matched_asm.append(query_res)
        topt_aligned.append(topt)
    elif check_tax:
        query_res = asm[asm.taxid==tax]
        variants.append(6)
        matched_asm.append(query_res)
        topt_aligned.append(topt)
    elif check_sptax:
        query_res = asm[asm.species_taxid==tax]
        variants.append(7)
        matched_asm.append(query_res)
        topt_aligned.append(topt)
    else:
        pass


# now as we extracted corresponding asm-entries, let's filter them by refseq_category=='representative genome' and  assembly_level=='Complete Genome' or 'Chromosome' ...
filtered_matched_asm = []
#
for asm_group in matched_asm:
    if asm_group.shape[0] > 1:
        comp_genome = asm_group.assembly_level.isin(['Complete Genome',]).any()
        chromosome = asm_group.assembly_level.isin(['Chromosome',]).any()
        # contig = asm_group.assembly_level.isin(['Contig',]).any()
        rep_stat = asm_group.refseq_category.isin(['representative genome',]).any()
        version_status = asm_group.version_status.isin(['latest',]).any()
        # 
        if comp_genome and rep_stat:
            filtered_matched_asm.append(asm_group[(asm_group.assembly_level=='Complete Genome')&(asm_group.refseq_category=='representative genome')].iloc[0:1])
        elif chromosome and rep_stat:
            filtered_matched_asm.append(asm_group[(asm_group.assembly_level=='Chromosome')&(asm_group.refseq_category=='representative genome')].iloc[0:1])
        elif comp_genome:
            filtered_matched_asm.append(asm_group[(asm_group.assembly_level=='Complete Genome')].iloc[0:1])
        elif chromosome:
            filtered_matched_asm.append(asm_group[(asm_group.assembly_level=='Chromosome')].iloc[0:1])
        elif rep_stat:
            filtered_matched_asm.append(asm_group[(asm_group.refseq_category=='representative genome')].iloc[0:1])
        else:
            print "There are no complete or representative pieces of DNA ..."
            print "taking latest ..."
            latest_mask = asm_group.version_status=='latest'
            if latest_mask.sum()==1:
                filtered_matched_asm.append(asm_group[latest_mask].iloc[0:1])
            else:
                asm_tmp = copy.deepcopy(asm_group)
                asm_tmp['seq_date'] = pd.to_datetime(asm_tmp.seq_rel_date)
                filtered_matched_asm.append(asm_tmp.sort('seq_date',ascending=False,inplace=False).iloc[0:1])
    else:
        filtered_matched_asm.append(asm_group.iloc[0:1])




# after all the ugly manipulations we finally get what we wanted ...
dat_topt = pd.concat(filtered_matched_asm)
dat_topt['OptimumTemperature'] = topt_aligned
dat_topt.to_csv("assembly_topt.dat",index=False)

# write file with ftp requests ...
with open('ftp_request.dat','w') as fp:
    for ppp in ggg.ftp_path:
        fp.write(ppp+'\n')


# original assembly file is somewhere in ftp.ncbi.nih.gov/genomes/genbank/archaea/
#

 # wget -b -r --user=anonymous --password=sergey.venev@umassmed.edu -i ftp_request.dat # log outputted to wget-log



