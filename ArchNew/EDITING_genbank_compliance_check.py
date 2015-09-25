from Bio import SeqIO
import pandas as pd
import os
from time import sleep
import subprocess as sub
import re
import sys
from Bio import SeqUtils


aacids = sorted(list('CMFILVWYAGTSNQDEHRKP'))
# aacids = sorted(aacids)
#
# def get_cds_count(gb):
#     cds_features_counter = 0
#     for feature in gb.features:
#         if (feature.type == 'CDS')and('translation' in feature.qualifiers.keys()):
#             cds_features_counter += 1
#     #return the results ...
#     return cds_features_counter
#
def extract_proteins(gb):
    proteome = ''
    for feature in gb.features:
        if (feature.type == 'CDS')and('translation' in feature.qualifiers.keys()):
            proteome += feature.qualifiers['translation'][0]
    #return the results ...
    return proteome
#
def check_seq(seq):
    nucs = sum([int(nuc in seq) for nuc in list('ATGC')])
    return nucs == len('ATGC')
#
# let's simpy avoid plasmids ...
plasmid_pattern = re.compile('plasmid')
def check_chromosome(gb):
    return not plasmid_pattern.search(gb.description)






root_path = "."
# # archaea/genomes
path = os.path.join(root_path,"ftp.ncbi.nlm.nih.gov/genomes/all")
# # asm_path = os.path.join(root_path,"archaea_links")
# result_path = os.path.join("..","archaea250")

asm = pd.read_csv(os.path.join(root_path,'assembly_topt.dat'))

#
#
# # data we are interested in:
# # (1)number of features in gb, (2)are there protein fasta file,
# # (3)nucleotide sequence in gb, (4) total number of genbank entries!
#
# # cds_features = []
# # num_proteins = []
# # nucseq_gb = []
# # gb_entries = []
# # accession = []
#
# GC = []
# # prot_dat = []
# aafs = {}
# for aa in aacids:
#     aafs[aa] = []
#
# genome_length = []
# proteome_length = []
#
#

# approximate data extraction chart ...
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



# iterate over all folders with assemblages and extract genomic information from their contents ...
for asm_suffix,asm_access in asm[['asm_name','assembly_accession']].itertuples(index=False):
    #
    asm_root_name = '_'.join(asm_access,asm_suffix)
    asm_dir = os.path.join(path, asm_root_name)
    gb_fname = os.path.join(asm_dir,'_'.join(asm_root_name,'genomic.gbff'))
    fa_fname = os.path.join(asm_dir,'_'.join(asm_root_name,'genomic.fna'))
    # check if both files exist ...
    gb_check = os.path.isfile(gb_fname)
    fa_check = os.path.isfile(fa_fname)
    if not(gb_check and fa_check):
        print "Genomic files do not exist in %s"%asm_dir
        sys.exit(1)
    # if files exist ...
    #
    # EXTRACT DESCRIPTION FROM FILES ...
    # TO BE CONTINUED ...







#################################################
##              OLD - DEPRECATED STUFF ...
#################################################
for asm in dat.assembly:
    # look up genome folder corresponding to the assembly...
    asm_folder_content = sub.check_output("ls %s|grep %s"%(os.path.join(data_path,asm),asm),shell=True)
    # check if there is a genbank file 
    if re.search("genomic.gbff",asm_folder_content):
        gb_fname = "%s/%s_genomic.gbff"%(os.path.join(data_path,asm),asm)
        genbank = SeqIO.parse(gb_fname,format='genbank')
        proteome = ''
        genome = ''
        for entry in genbank:
            if check_chromosome(entry):
                proteome += extract_proteins(entry)
                if check_seq(entry.seq):
                    genome += entry.seq
        # now filling in the arrays ...
        if not genome:
            if re.search("genomic.fna",asm_folder_content):
                fna_fname = "%s/%s_genomic.fna"%(os.path.join(data_path,asm),asm)
                fna = SeqIO.parse(fna_fname,format='fasta')
                genome = ''.join([str(fna_entry.seq) for fna_entry in fna])
            else:
                print "There are no nucleotide entries for assembly %s! Terminate!"%asm
                sys.exit(1)
        #################################
        # now filling the dat array with the data ...
        GC.append(SeqUtils.GC(genome))
        prot_len = len(proteome)
        # lengths ...
        genome_length.append(len(genome))
        proteome_length.append(prot_len)
        if prot_len:
            for aa in aacids:
                aafs[aa].append(proteome.count(aa)/float(prot_len))
        else:
            for aa in aacids:
                aafs[aa].append(0.0)
        #################################
    else:
        print "No genbank for assembly %s, Terminate!!!"%asm
        sys.exit(1)




# # dat['cds_features'] = cds_features
# # dat['num_proteins'] = num_proteins
# # dat['nucseq_gb'] = nucseq_gb
# # dat['gb_entries'] = gb_entries
# # dat['accession'] = accession


# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ###########################
# # Length of values does not match length of index
# ###########################
# dat['GC'] = GC
# for aa in aacids:
#     dat[aa] = aafs[aa]

# dat['genome_length'] = genome_length
# dat['proteome_length'] = proteome_length


# dat.to_csv('archaea_repeat.dat',index=False)











