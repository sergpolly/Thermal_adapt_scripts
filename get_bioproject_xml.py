import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez as ncbi
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import itertools

# always tell NCBI who ... are you!?
ncbi.email = "yourmail@server.com"

# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER/ftp.ncbi.nih.gov/refseq/release/bacteria/genbank')
path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')

descriptions_fname = "genbank.inventory.description"


descriptions = pd.read_csv(os.path.join(path,descriptions_fname))
# ignore low_memory=False warning ...


bio_match = re.compile("BioProject:(.+)")

def extract_bio_uid(dbxref):
	uids = str(dbxref).strip().replace(';',' ').split(' ')
	matched_uids = [bio_match.match(uid) for uid in uids if bio_match.match(uid)]
	# bio_match.match(uid) is executed twice here, but who cares ...
	if matched_uids:
		accession = matched_uids[0].groups()[0]
		uid = re.match("[A-Z]+(\d+)",accession).groups()[0]
		return uid
		# return just a single BioProject uid, I believe that would be correct 100% of times
	else:
		return np.nan


# extract solely BioProject uid first ...
descriptions['BioProject'] = descriptions['DbxRefs'].apply(extract_bio_uid)

# # write down the extracted BioProject uids back to the file ...
# descriptions.to_csv(os.path.join(path,descriptions_fname), index=False)


# find unique BioProject UIDs ...
bio_uniq_uids = descriptions['BioProject'].unique()
# that are not nan/null 
bio_uniq_uids = bio_uniq_uids[pd.notnull(bio_uniq_uids)]


####################### OLD STUFF ...
# def parse_bio_uid(uid):
# 	# assuming uid as PRJNA00000
# 	uid_loc = str(uid)
# 	res = re.match('[A-Z]+(\d+)',uid_loc)
# 	if res:
# 		return res.groups()[0]
# 	else:
# 		return ''
####################### OVER


print >> sys.stdout, "Sending query of %d BioProject uids to NCBI ..."%bio_uniq_uids.size


################ TODO:
# find unique , not null BioProject entries and fetch their corresponding xmls from NCBI ...
# for each file in infobacter_fnames do ...
#
#

# fetching from NCBI: the entire unique BioProject uids:
handle = ncbi.epost(db='bioproject',id=','.join(bio_uniq_uids))
results = ncbi.read(handle)
# efetch query parameters:
counts, web_env, query_key = bio_uniq_uids.size, results['WebEnv'], results['QueryKey']
#
print >> sys.stdout, "NCBI respose received."
print >> sys.stdout, ""
print >> sys.stdout, "Preparing to fetch %d results ..."%bio_uniq_uids.size
# fetch the whole batch ...
batch = 50
with open(os.path.join( path, 'results.xml' ) , "w") as out_handle:
	for i_batch in range(0,counts,batch):
	    end = min(counts, i_batch + batch)
	    print "Fetching records %i thru %i" % (i_batch + 1, end)
	    fetch_handle = ncbi.efetch(db="bioproject", retstart=i_batch, retmax=batch, webenv=web_env, query_key=query_key)
	    data = fetch_handle.read()
	    fetch_handle.close()
	    out_handle.write(data)

print >> sys.stdout, "Fetching completed."

#
#
#








