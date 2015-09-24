import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez as ncbi
import xml.etree.ElementTree as ET
import pandas as pd
import itertools
import numpy as np
import copy

# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER/ftp.ncbi.nih.gov/refseq/release/bacteria/genbank')
path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')

env_fname = "env_catalog_bacter.dat"
descriptions_fname = "genbank.inventory.description"

env_df = pd.read_csv(os.path.join(path,env_fname))
descriptions = pd.read_csv(os.path.join(path,descriptions_fname))

# Expect some silly warnings, like this:
# Columns (11) have mixed types

bio_match = re.compile("BioProject:(.+)")
taxon_match = re.compile("taxon:(.+)")

def extract_bio_uid(dbxref):
    uids = str(dbxref).strip().replace(';',' ').split(' ')
    matched_uids = [bio_match.match(uid) for uid in uids if bio_match.match(uid)]
    # bio_match.match(uid) is executed twice here, but who cares ...
    if matched_uids:
        accession = matched_uids[0].groups()[0]
        uid = re.match("[A-Z]+(\d+)",accession).groups()[0]
        return int(uid)
        # return just a single BioProject uid, I believe that would be correct 100% of times
    else:
        return np.nan

def extract_taxon(dbxref):
    uids = str(dbxref).strip().replace(';',' ').split(' ')
    matched_txs = [taxon_match.match(uid) for uid in uids if taxon_match.match(uid)]
    # taxon_match.match(uid) is executed twice here, but who cares ...
    if matched_txs:
        uid = matched_txs[0].groups()[0]
        return int(uid)
        # return just a single taxon uid, I believe that would be correct 100% of times
    else:
        return np.nan




# extract solely BioProject uid first ...
descriptions['BioProject'] = descriptions['DbxRefs'].apply(extract_bio_uid)

# extract solely Taxonomy uid first ...
descriptions['TaxonID'] = descriptions['SourceDbxRefs'].apply(extract_taxon)


# make BioProject UID an index 
env_df = env_df.set_index('BioProject', drop=True, append=False)

# DB-style outter join of the tables
des = descriptions.join(env_df, on='BioProject', lsuffix='_des', rsuffix='_env')

# Entries with known temperature ...
des_temp = des[(des.OptimumTemperature.notnull())]

# get those that are not plasmids
noplasm_idx = des_temp.SourcePlasmid.isnull().nonzero()[0]
des_noplasmid = copy.deepcopy(des_temp.iloc[noplasm_idx])

# criteria for the Complete Genome ...
crit_NC = lambda name : True if (name[:3]=='NC_') else False
crit_comp_genome = lambda descr: True if re.search('complete genome',descr) else False
crit_features = lambda nf: True if (nf>=900) else False
crit_genlen = lambda genlen: True if (genlen>=300000) else False
crit_WGS = lambda keywords: False if re.search('WGS',keywords) else True
crit_plasmid = lambda descr: False if re.search('plasmid',descr) else True

# evaluate those criteria ...
des_noplasmid['crit_NC'] = des_noplasmid.GenomicID.apply(crit_NC)
des_noplasmid['crit_WGS'] = des_noplasmid.Keywords.apply(crit_WGS)
des_noplasmid['crit_genlen'] = des_noplasmid.GenomicLen.apply(crit_genlen)
des_noplasmid['crit_features'] = des_noplasmid.FeaturesNum.apply(crit_features)
des_noplasmid['crit_comp_genome'] = des_noplasmid.Description.apply(crit_comp_genome)
des_noplasmid['crit_plasmid'] = des_noplasmid.Description.apply(crit_plasmid)

# get those items that satisfy these criteria ...
criteria_to_satisfy = ['crit_NC', 'crit_WGS', 'crit_genlen', 'crit_features', 'crit_comp_genome', 'crit_plasmid']
indexes_satisfy = (des_noplasmid[criteria_to_satisfy].sum(axis=1)==len(criteria_to_satisfy)).nonzero()[0]
# des_valid = copy.deepcopy(des_noplasmid.iloc[indexes_satisfy].OptimumTemperature.notnull())
des_valid = copy.deepcopy(des_noplasmid.iloc[indexes_satisfy])


# now collapse non-unique things,
# there shouldn't be many of them, at this point
# [CHECK] - just 3 of them found for the RELEASE69 ...
cleaned_data = des_valid.drop_duplicates(subset='TaxonID')


check_bp = cleaned_data.duplicated(subset='BioProject').nonzero()[0].size
check_org1 = cleaned_data.duplicated(subset='Organism_des').nonzero()[0].size
check_org2 = cleaned_data.duplicated(subset='Organism_env').nonzero()[0].size

if (check_bp+check_org1+check_org2)==0:
    print "Data looks extremely clean: no duplicates, complete genomes only, opt temps available ..."
elif check_bp:
    print "WARNING: There are duplicated BioProject UIDs in the final data"
elif check_org1:
    print "WARNING: There are duplicated Organisms in the final data"
elif check_org2:
    print "WARNING: There are duplicated Organisms in the final data"


# some information summary ...
thermophiles_counts = (cleaned_data.OptimumTemperature>=50).nonzero()[0].size
total_counts = cleaned_data.shape[0]
print "Cleaned data summary: %d thermophiles (OGT>=50) out of %d species in total."

# GenomicID is the best candidate to be an index ...
cleaned_data.to_csv("env_catalog_compgenome.dat",index=False)








