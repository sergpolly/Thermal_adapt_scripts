import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
import numpy as np
import pandas as pd
import itertools


path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER/ftp.ncbi.nih.gov/refseq/release/bacteria/genbank')


# this script inventories all GbDb entries, collecting their description, organism name,
# number of annotated features, genome length, etc.
# all this information is stored in a single csv file for convenience.

# (!)Many of the entries belong to the same organism, because some of the genomes
# are split in pieces etc.

inv_fname = sys.argv[1]

genome_inventory = os.path.join( path, inv_fname )
# read genome inventory ...
with open(genome_inventory,'r') as fp:
    genome_files = [line.strip() for line in fp.readlines()]




print "Flat processing of the inventoried files got started ..."
print

data = {'GenomicID':[],
        'GenomicName':[],
        'Organism':[],
        'Description':[],
        'FeaturesNum':[],
        'GenomicLen':[],
        'NucsPresent':[],
        'Keywords':[],
        'DbxRefs':[]
        }


# methods of extraction of SeqRec info
get_organism = lambda seqrec: seqrec.annotations['organism'].replace(',',' ')
get_descr = lambda seqrec: seqrec.description.replace(',',' ')
get_feat_num = lambda seqrec: len(seqrec.features)
get_seqlen = lambda seqrec: len(seqrec.seq)
get_if_nucs = lambda seq: sum([seq.find(let) for let in list('ATGC')]) > 0
get_kws = lambda seqrec: ';'.join(seqrec.annotations['keywords'])
get_dbx = lambda seqrec: ';'.join(seqrec.dbxrefs)

total_counter = 0
for gb_fname in genome_files:
    gb_genomes = SeqIO.parse( os.path.join(path,gb_fname) ,format='genbank')
    print >> sys.stdout, "Processing file %s ..."%gb_fname
    # collect annotations for all the genomes inside each genbank file ... 
    file_counter = 0
    for seqrec in gb_genomes:
        # data dictionary must be filled 
        data['GenomicID'].append(seqrec.id)
        data['GenomicName'].append(seqrec.name)
        data['Organism'].append(get_organism(seqrec))
        data['Description'].append(get_descr(seqrec))
        data['FeaturesNum'].append(get_feat_num(seqrec))
        data['GenomicLen'].append(get_seqlen(seqrec))
        data['NucsPresent'].append(get_if_nucs(seqrec.seq))
        data['Keywords'].append(get_kws(seqrec))
        data['DbxRefs'].append(get_dbx(seqrec))
        file_counter += 1
        total_counter += 1
    print >> sys.stdout, "%d entries are processed from the file."%file_counter
    print >> sys.stdout, ""


print >> sys.stdout, " All files are processed: %d entries inventoried in total"%total_counter


annot_df = pd.DataFrame(data)
# columns = ['GenomicID', 'GenomicName', 'Organism', 'Description', 'FeaturesNum', 'GenomicLen', 'NucsPresent', 'Keywords', 'DbxRefs']

annot_df.to_csv(os.path.join( path, inv_fname+'.description' ),index=False)































