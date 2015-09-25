import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
import numpy as np
import pandas as pd
import itertools


# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER/ftp.ncbi.nih.gov/refseq/release/bacteria/genbank')
path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')


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

data = {'GenomicID':[], # id extracted from SeqRecord
        'GenomicName':[], # name extracted from SeqRecord
        'Organism':[], # Organism name extracted from annotations 
        'Description':[], # Description from annotations
        'FeaturesNum':[], # Number of features, including CDS,RNA,genes, etc.
        'GenomicLen':[], # Length of the genomic piece
        'NucsPresent':[], # Is there nucleotide sequence in the genbank 
        'Keywords':[], # Keywords from annotations 
        'DbxRefs':[], # DataBase references extracted from annotations
        'Taxonomy':[], # Taxonomy information from annotations
        'SourceDbxRefs':[], # Taxonomy ID, extracted from the features[0], which describes SOURCE 
        'SourceOrganism':[], # Organism name extracted form features[0]
        'SourceStrain':[], # Strain extracted form features[0]
        'SourcePlasmid':[] # Check if it is a plasmid, from features[0]
        }


# one can also extract things like taxon id,organism name, and strain is this piece is a plasmid ...
# {'db_xref': ['taxon:1923'],
#  'mol_type': ['genomic DNA'],
#  'organism': ['Streptomyces phaeochromogenes'],
#  'plasmid': ['pJV1'],
#  'strain': ['NRRL-B3559']}

# that would be seqrec.features[0].qualifiers['db_xref']
#               seqrec.features[0].qualifiers['organism']
#               seqrec.features[0].qualifiers['strain']
#               seqrec.features[0].qualifiers['plasmid']
# 'plasmid' in seqrec.features[0].keys()
# also make sure seqrec.features[0] is of type 'source' prior of doing any extractions ...
# if (seqrec.features[0].type == 'source')



# methods of extraction of SeqRec info
get_organism = lambda seqrec: seqrec.annotations['organism'].replace(',',' ')
get_descr = lambda seqrec: seqrec.description.replace(',',' ')
get_feat_num = lambda seqrec: len(seqrec.features)
get_seqlen = lambda seqrec: len(seqrec.seq)
get_if_nucs = lambda seq: sum([seq.find(let) for let in list('ATGC')]) > 0
get_kws = lambda seqrec: ' '.join(seqrec.annotations['keywords'])
get_dbx = lambda seqrec: ' '.join(seqrec.dbxrefs)
get_taxon = lambda seqrec: ';'.join(seqrec.annotations['taxonomy']) if ('taxonomy' in seqrec.annotations) else np.nan
get_source_dbx = lambda feature: ' '.join(feature.qualifiers['db_xref']) if ('db_xref' in feature.qualifiers) else np.nan
get_source_org = lambda feature: ' '.join(feature.qualifiers['organism']) if ('organism' in feature.qualifiers) else np.nan
get_source_strain = lambda feature: ' '.join(feature.qualifiers['strain']) if ('strain' in feature.qualifiers) else np.nan
get_source_plasmid = lambda feature: ' '.join(feature.qualifiers['plasmid']) if ('plasmid' in feature.qualifiers) else np.nan

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
        data['Taxonomy'].append(get_taxon(seqrec))
        #
        first_feature =  seqrec.features[0]
        if first_feature.type == 'source':
            data['SourceDbxRefs'].append(get_source_dbx(first_feature))
            data['SourceOrganism'].append(get_source_org(first_feature))
            data['SourceStrain'].append(get_source_strain(first_feature))
            data['SourcePlasmid'].append(get_source_plasmid(first_feature))
        else:
            data['SourceDbxRefs'].append(np.nan)
            data['SourceOrganism'].append(np.nan)
            data['SourceStrain'].append(np.nan)
            data['SourcePlasmid'].append(np.nan)
        file_counter += 1
        total_counter += 1
    print >> sys.stdout, "%d entries are processed from the file."%file_counter
    print >> sys.stdout, ""


print >> sys.stdout, " All files are processed: %d entries inventoried in total"%total_counter


annot_df = pd.DataFrame(data)
# columns = ['GenomicID', 'GenomicName', 'Organism', 'Description', 'FeaturesNum', 'GenomicLen', 'NucsPresent', 'Keywords', 'DbxRefs']

annot_df.to_csv(os.path.join( path, inv_fname+'.description' ),index=False)































