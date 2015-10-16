from Bio import SeqIO
import pandas as pd
import os
from time import sleep
import subprocess as sub
import re
import sys
from Bio import SeqUtils


root_path = "."
# # archaea/genomes
path = os.path.join(root_path,"ftp.ncbi.nlm.nih.gov/genomes/all")
# # asm_path = os.path.join(root_path,"archaea_links")
# result_path = os.path.join("..","archaea250")

asm = pd.read_csv(os.path.join(root_path,'assembly_topt.dat'))


# approximate data extraction chart ...
data = {# there might be several GenBank per assembly, so ...
        'UniqueAsmId':[], # assembly id for later joining ...
        'GenBankItems':[], # number of SeqRecords per assembly ...
        #
        'GenomicIDx_list':[], # id-s extracted from SeqRecords ...
        'GenomicNames_list':[], # name per SeqRecord ..
        #
        'Organism_list':[], # just a single organism, though ... (from annotations)
        'Descriptions_list':[], # Descriptions from annotations
        #
        'FeaturesTotal':[], # Number of features, including CDS,RNA,genes, etc.
        'FeaturesByGenBank_list':[], # Number of features, per GenBank
        #
        'CDSTotal':[], # Number of CDS per assembly - total
        'CDSByGenBank_list':[], # Number of CDS, per GenBank in assembly 
        #
        'GenomicLenTotal':[], # Total Length
        'GenomicLenByGenBank_list':[], # Length per GenBank (SeqRecord) ...
        #
        'NucsPresent':[], # AND applied for all of the pieces ... 
        'Keywords_list':[], # Keywords from annotations (combined) ...
        # DataBase references are already in assembly_summary file ...
        # Taxonomy information is already in assembly_summary file ...
        'PlasmidCheck_list':[] # Check if it is a plasmid, from features[0] of each SeqRecord ...
        }
# Similar field is used in the Bacterial pipeline ...
# # approximate data extraction chart ...
# data = {'GenomicID':[], # id extracted from SeqRecord
#         'GenomicName':[], # name extracted from SeqRecord
#         'Organism':[], # Organism name extracted from annotations 
#         'Description':[], # Description from annotations
#         'FeaturesNum':[], # Number of features, including CDS,RNA,genes, etc.
#         'GenomicLen':[], # Length of the genomic piece
#         'NucsPresent':[], # Is there nucleotide sequence in the genbank 
#         'Keywords':[], # Keywords from annotations 
#         'DbxRefs':[], # DataBase references extracted from annotations
#         'Taxonomy':[], # Taxonomy information from annotations
#         'SourceDbxRefs':[], # Taxonomy ID, extracted from the features[0], which describes SOURCE 
#         'SourceOrganism':[], # Organism name extracted form features[0]
#         'SourceStrain':[], # Strain extracted form features[0]
#         'SourcePlasmid':[] # Check if it is a plasmid, from features[0]
#         }

# iterate over all folders with assemblages and extract genomic information from their contents ...
for asm_suffix,asm_access in asm[['asm_name','assembly_accession']].itertuples(index=False):
    #
    asm_suffix = asm_suffix.replace(' ','_')
    #
    asm_root_name = '_'.join([asm_access,asm_suffix])
    asm_dir = os.path.join(path, asm_root_name)
    gb_fname = os.path.join(asm_dir,'_'.join([asm_root_name,'genomic.gbff']))
    fa_fname = os.path.join(asm_dir,'_'.join([asm_root_name,'genomic.fna']))
    # check if both files exist ...
    gb_check = os.path.isfile(gb_fname)
    fa_check = os.path.isfile(fa_fname)
    if not(gb_check and fa_check):
        print "Genomic files do not exist in %s"%asm_dir
        sys.exit(1)
    # if files exist ...
    #
    SeqRecs = SeqIO.parse(gb_fname,format='genbank')
    # initialize description fields ...
    GenBankItems = 0
    GenomicIDx = []
    GenomicNames = []
    Organism = []
    Descriptions = []
    FeaturesByGenBank = []
    CDSByGenBank = []
    GenomicLenByGenBank = []
    NucsPresent = True
    Keywords = []
    PlasmidCheck = []
    # start extraction ...
    for SeqRec in SeqRecs:
        GenBankItems += 1
        GenomicIDx.append(SeqRec.id)
        GenomicNames.append(SeqRec.name)
        Organism.append(SeqRec.annotations['organism'].replace(',',' '))
        Descriptions.append(SeqRec.description.replace(',',' '))
        FeaturesByGenBank.append( len(SeqRec.features) )
        # count CDS features only ...
        CDSByGenBank.append( len( [feature for feature in SeqRec.features if feature.type=='CDS'] ) )
        GenomicLenByGenBank.append( len(SeqRec.seq) )
        NucsPresent = ( NucsPresent and (sum([SeqRec.seq.find(let) for let in list('ATGC')]) > 0) )
        Keywords.append(';'.join(SeqRec.annotations['keywords']))
        #
        # first feature is usually description of the SOURCE ...
        first_feature =  SeqRec.features[0]
        #
        if first_feature.type == 'source':
            PlasmidCheck.append(';'.join(first_feature.qualifiers['plasmid']) if ('plasmid' in first_feature.qualifiers) else '')
        else:
            PlasmidCheck.append(np.nan)
        #
    #
    # Now for all of these SeqRecords (that belong to same Assembly) output the summary and store it ...
    data['UniqueAsmId'].append(asm_access)
    data['GenBankItems'].append(GenBankItems)
    data['GenomicIDx_list'].append(';'.join(GenomicIDx))
    data['GenomicNames_list'].append(';'.join(GenomicNames))
    data['Organism_list'].append(';'.join(Organism))
    data['Descriptions_list'].append(';'.join(Descriptions))
    data['FeaturesTotal'].append( sum(FeaturesByGenBank) )
    data['FeaturesByGenBank_list'].append( ';'.join(str(_) for _ in FeaturesByGenBank) )
    # CDS business ...
    data['CDSTotal'].append( sum(CDSByGenBank) )
    data['CDSByGenBank_list'].append( ';'.join(str(_) for _ in CDSByGenBank) )
    #
    data['GenomicLenTotal'].append( sum(GenomicLenByGenBank) )
    data['GenomicLenByGenBank_list'].append( ';'.join(str(_) for _ in GenomicLenByGenBank) )
    data['NucsPresent'].append(NucsPresent)
    data['Keywords_list'].append(';'.join(Keywords))
    data['PlasmidCheck_list'].append(';'.join(PlasmidCheck))
    #



annot_df = pd.DataFrame(data)
#   TODO:
# it appears to be working ...
# just check the output and analyze the data ....
#
# The major criteria would be the availability of the translated proteins ...
# However we must avoid plasmid genes (!)
# GenBankItems number matches up exactly with the 'PlasmidCheck_list'.apply(len(_.split(';'))),
# so go ahead and use split(';') no problem.
# so then whenever an item/assembly has >= XXX translated proteins (excluding plasmid based ones),
# simply accept that item aboard! It might be a good idea to check that the number of translated protein
# per plasmid does not exceed some YYY as well.


PLASMID_LIMIT = 650
# there are really big archaeal plasmids out there !!! see GCA_000517625.1 and GCA_000025685.1 for instance ...
PROTEIN_LIMIT = 600
# iterate over lists of plasmid-checks and CDS numbers extracted ...
total_cds_excl_plasmid = []
genome_len_excl_plasmid = []
for asm_id,plasmid_list,cds_list,genlen_list in annot_df[['UniqueAsmId','PlasmidCheck_list','CDSByGenBank_list','GenomicLenByGenBank_list']].itertuples(index=False):
    # split back the ;-separated strings ...
    plasmid_list = plasmid_list.split(';')
    cds_list = [int(_) for _ in cds_list.split(';')]
    genlen_list = [int(_) for _ in genlen_list.split(';')]
    # and analyze it ...
    total_cds = 0
    genlen_no_plasmid = 0
    for p,cds,genlen in zip(plasmid_list,cds_list,genlen_list):
        # if it's a plasmid - see it is suspicious one ...
        if (p!="")and(cds>=PLASMID_LIMIT):
            print "There appears to be too many CDS (%d) in the plasmid %s from %s ..."% (cds,p,asm_id)
            sys.exit(1)
        # if it's not a plasmid increment CDS counter ...
        # and count genome length excluding plasmids  ...
        if p=="":
            total_cds += cds
            genlen_no_plasmid += genlen
    #
    total_cds_excl_plasmid.append(total_cds)
    genome_len_excl_plasmid.append(genlen_no_plasmid)
    #

# so these are the final CDS counts for the final filtering of relevant organisms ...
annot_df['CDSNoPlasmid'] = total_cds_excl_plasmid
annot_df['GenomicLenNoPlasmid'] = genome_len_excl_plasmid


# so finally the total number of the organisms we were able to recover from GenBank ...
number_recovered_organisms = (annot_df['CDSNoPlasmid'] > PROTEIN_LIMIT).nonzero()[0].size
print "The total number of archaea with complete genomes (CDS>600) and known environmental T is %d"%number_recovered_organisms
# Selecting interesting organisms only, ones with enough translated CDSs to work with ...
interesting_orgs_idx = (annot_df['CDSNoPlasmid'] > PROTEIN_LIMIT)
# items selected this way are perfectly aligned in terms of assembly accession numbers (unique identifiers in our case) ...
interest_annot = annot_df[interesting_orgs_idx]
interest_asm = asm[interesting_orgs_idx]


# Next thing is to check if all the nucleotide sequences are directly accessible from the GenBank files!
nucs_present_all = interest_annot['NucsPresent'].all()
if not nucs_present_all:
    print "Modify downstream scripts to allow nucleotide sequences from fasta files, as not every GenBank has nucleotides available ..."
    sys.exit(1)


# if everything runs smooth, we should output some kind of summary in the end ...
# make sure we are avoiding plasmid sequences for the downstream analysis ...


# now we'd need to figure out what columns do we wnat to keep both from 'interest_asm' and from 'interest_annot' arrays ...
# ...
# import columns from asm:
annot_cols = ['UniqueAsmId','CDSNoPlasmid','GenomicLenNoPlasmid','NucsPresent','GenBankItems','PlasmidCheck_list']
asm_cols   = ['asm_name', 'assembly_accession', 'assembly_level', 'bioproject', 'biosample', 'ftp_path', 'organism_name', 'species_taxid', 'taxid', 'OptimumTemperature']

# join or merge em by 'UniqueAsmId'(in annot) and 'assembly_accession'(in asm) - they correspond to each other perfectly and are perfectly aligned as well ...
#
# Web-Github is alos great because it can index all the code existing in the project, so it's easy to look for
# the way, let's say, a particular function was applied ...
# and, as we remember 'pandas.join' is a legacy thing, so we keep using 'pandas.merge' instead!

complete_summary = pd.merge(interest_asm[asm_cols],
                            interest_annot[annot_cols],
                            how='left',
                            left_on='assembly_accession',
                            right_on='UniqueAsmId')


# we can output the table finally ...
out_summary_fname = os.path.join(root_path,'summary_organisms_interest.dat')
print "everything is ready, saving final summary of interesting organisms to %s ..."%out_summary_fname
complete_summary.to_csv(out_summary_fname, index=False)











