import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez as ncbi
from Bio import SeqUtils
from Bio.SeqUtils import CodonUsage
# import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
# import itertools
import math
# import matplotlib.pyplot as plt
from Bio import Data

CODON_LEN = 3
RIBO_LIMIT = 20

# always tell NCBI who ... are you!?
ncbi.email = "sergey.venev@umassmed.edu"




# imac_path = '/Users/venevs/Desktop/archaea2015'
# mbook_path = '/Users/sergpolly/Desktop/Dropbox (UMASS MED - BIB)/archaea2015'
bib2_data_path = '/home/venevs/NCBI_BACT_DATA'
bib2_scr_path = '/home/venevs/NCBI_BACT_PROCESS'

# genbank files ...
gbdb_path = '/home/venevs/NCBI_BACT_DATA'
gbdb_fname = os.path.join(gbdb_path,'genbank.idx')
gbdb = SeqIO.index_db(gbdb_fname)

# nucleotides fasta files ...
# res = SeqIO.index_db(db_file,filenames=fnames,format='fasta',key_function=(lambda name: name.split('|')[3]))
fndb_path = '/home/venevs/NCBI_BACT_DATA/nucfasta'
fndb_fname = os.path.join(fndb_path,'nucfasta.idx')
fndb = SeqIO.index_db(fndb_fname,key_function=(lambda name: name.split('|')[3]))

# we'll calculate protein level CAI and store em in individual files ...
path_CAI = os.path.join(bib2_scr_path,'genomes_with_CAI')


#
# gbdb and fndb keys must be the same ...


# file with the organisms of interest 
org_dat_fname = os.path.join(bib2_scr_path,'catalog_with_accesion.dat')
org_dat = pd.read_csv(org_dat_fname)




###################################################################################
# CAI related functions ...
ribosomal_check = re.compile('ribosomal protein')
def extract_ribo_genes(genbank,nuc_seq):
    ribo_genes = []
    rejected_count = 0
    for feature in genbank.features:
        qual_keys = feature.qualifiers.keys()
        if (feature.type == 'CDS')and('translation' in qual_keys)and('product' in qual_keys):
            if ribosomal_check.search(feature.qualifiers['product'][0]):
                # this is the SUPER STRICT check for suitable ribo-genes ... 
                # extracting coding DNA sequences ...
                extracted_nucs = feature.extract(nuc_seq)
                # try to translate it back using some genetic code, which one ?
                genetic_table = feature.qualifiers['transl_table'][0] if ('transl_table' in qual_keys) else 11
                # genetic_table = 11 is the Bacterial, Archaeal and Plant Plastid Code, so it'll be our default, just in case ...
                try:
                    translation = extracted_nucs.translate(table=genetic_table,cds=True)
                # catch the translational error exception ...
                except Data.CodonTable.TranslationError:
                    # if gene is not true-CDS and it cannot be translated we'll automatically reject it!
                    rejected_count += 1
                    # and move on to the next gene in the loop ...
                    continue
                # even if it is a true-CDS, we'll check if it translates back to what is in the feature qualifiers ..
                if str(translation)==feature.qualifiers['translation'][0]:
                    # OK! extracted gene translates to the the provided feature qualifier ...
                    ribo_genes.append(extracted_nucs)
                else:
                    # rejecting the gene if it does not match with the feature qualifier translation ...
                    rejected_count += 1
    return (ribo_genes,rejected_count)

# this returns an array of tuples with assorted info on the CDS-es from a given genbank ...
def extract_genes_features(genbank,nuc_seq):
    genes_n_features = []
    rejected_count = 0
    for feature in genbank.features:
        feat_quals = feature.qualifiers 
        qual_keys = feat_quals.keys()
        if (feature.type == 'CDS')and('translation' in qual_keys)and('product' in qual_keys):
            # we'll be extracting only nice CDS-es, ones that can be translatet back to the translation present in features ...
            # extract coding nucleotide sequence:
            gene = feature.extract(nuc_seq)
            # try to translate it back using some genetic code, which one ?
            genetic_table = feat_quals['transl_table'][0] if ('transl_table' in qual_keys) else 11
            # genetic_table = 11 is the Bacterial, Archaeal and Plant Plastid Code, so it'll be our default, just in case ...
            try:
                translation = gene.translate(table=genetic_table,cds=True)
            # we'll be catching TranslationError exceptions here ...
            except Data.CodonTable.TranslationError:
                # if gene is not true-CDS and it cannot be translated we'll automatically reject it!
                rejected_count += 1
                # and move on to the next gene in the loop ...
                continue
            # even if it is a true-CDS, we'll check if it translates back to what is in the feature qualifiers ..
            if str(translation)==feat_quals['translation'][0]:
                # OK! extracted gene translates to the the provided feature qualifier ...
                prot, = feat_quals['translation']
                prot_id, = feat_quals['protein_id']
                product, = feat_quals['product']
                product = product.replace(',',' ')
                genes_n_features.append((prot_id,product,gene,prot))
            else:
                # rejecting the gene if it does not match with the feature qualifier translation ...
                rejected_count += 1
    return (genes_n_features,rejected_count)

def count_codons(dna_seq_list,cod_dict):
    # cod_dict must be zeroed out beforehand, that's a potential issue ...
    # codons = [dna_seq[i:i+CODON_LEN] for i in range(0,len(dna_seq),CODON_LEN)]
    for gene_sequence in dna_seq_list:
        gene_sequence = str(gene_sequence)
        len_gene_sequence = len(gene_sequence)
        if len_gene_sequence%CODON_LEN:
            print >> sys.stderr,"(Counting codons in ribo-genes) gene lengths is not ~3! Skipping it!!!"
        else:
            for codon in [gene_sequence[i:i+CODON_LEN] for i in range(0,len_gene_sequence,CODON_LEN)]:
                if codon in cod_dict.keys(): 
                    cod_dict[codon] += 1 
                else: 
                    # raise TypeError("Illegal codon %s in genes" % codon)
                    print >> sys.stderr, "(Counting codons in ribo-genes) Illegal codon %s in genes! Skipping it!" % codon

# CodonUsage.SynonymousCodons
def generate_codon_index(cod_dict):
    index = {}
    # now to calculate the index we first need to sum the number of times 
    # synonymous codons were used all together. 
    # here is the place for potential improvement, because CodonUsage.SynonymousCodons implies the standard genetic code,
    # while table=4 gives up 1 STOP codon in favor of a W(triptophane) one! While, table=11 uses the standart code, with the 
    # expanded range of start codons ... 
    for aa in CodonUsage.SynonymousCodons: 
        SynCodons = CodonUsage.SynonymousCodons[aa]
        ###################################
        # total number of codons coding for a given amino acid ...
        total = sum([cod_dict[scodon] for scodon in SynCodons])
        if not total:
            print >> sys.stderr, "ACHTUNG!!! this is unprecedented!!!"
            print >> sys.stderr, "Codons coding for %s amino acid are not encountered in ribo-genes of this organism!" % aa
            print >> sys.stderr, "Codons corresponding to %s will appear as Illegal for this organism!" % aa            
            continue
        else:
            ###################################
            # calculate the RSCU value for each of the codons
            # RSCU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons) 
            # Relative Synonymous Codon Usage (RSCU)...
            # (sum of all synonymous codons/num of synonymous codons) is expected number of this-type-of-codon assuming no codon bias!
            no_bias_usage = float(total)/len(SynCodons)
            RSCU = [ cod_dict[scodon]/no_bias_usage for scodon in SynCodons ]
            ###################################
            # now generate the index W=RCSUi/RCSUmax: 
            RSCU_max = max(RSCU)
            for scodon,RSCUi in zip(SynCodons,RSCU):
                index[scodon] = RSCUi/RSCU_max
            ####################################
    return index

def cai_for_gene(dna_seq,index): 
    dna_seq = str(dna_seq) 
    cai_value, cai_length, illegal_codons_count = 0, 0, 0
    for codon in [dna_seq[i:i+CODON_LEN] for i in range(0,len(dna_seq),CODON_LEN)]:
        if codon in index: 
            # these two codons are always one, exclude them
            if codon not in ['ATG', 'TGG']:
                if index[codon]>0.0:
                    cai_value += math.log(index[codon]) 
                    cai_length += 1 
        else:
            illegal_codons_count += 1
            print >> sys.stderr, "skipping illegal codon in sequence: %s" % codon
            # our index DOES include STOP codons so, we'll just comment 2 lines out ...
            # elif codon not in ['TGA', 'TAA', 'TAG']:  # some indices may not include stop codons 
            # raise TypeError("illegal codon in sequence: %s.\n%s" % (codon, index))
    if illegal_codons_count:
        print >> sys.stderr, "%d illegal codons were skipped in the %d nucs gene! %s" % (illegal_codons_count,len(dna_seq),dna_seq)        
    return math.exp(cai_value / (cai_length - 1.0)) 
###################################################################################











# ALL organisms appear to be unique in the catalog, so that there are NO double-chromosome bacteria in our dataset...
# that is strange, because I remember that I did encounter them somehow and that there are abnormally low amount of ribo-genes there ...
# so anyways, even if I'm missing something, hopefully this two effects would cancel each other ...

# next, we'll simply follow the same approach that we pursued with Archaea to extract genes and their CAIs ...
# In the Archeal case, we used assembly code as a unique Organism identifier, whereas in the case of Bacteria, 
# we limited ourselfes to the analysis of the complete genomes only (omitting WholeGenomeShotguns, and multi-chromosomal Organisms), 
# so, in Bacterial case, it is the genome accession number that is the unique Orgnismal identifier. (accession number = IndexId in our case ...)
#########################################################
#########################################################
#(1) for each IndexId - find ribosomalome genes and get the codon usage of highly expressed genes ...
#(2) for each IndexId calculate CAI of every coding gene (CDS) in the genome
#(3) output detailed info on protein/prot_id/prot_seq/CAI/etc for each IndexId in a separate file ...
#(4) done!!!
#########################################################
#########################################################
#########################################################
# remark: could not find sufficient number of ribosomal genes -> no CAI!
# meaning that CAI is only 'defined', in case Ribosomalome is 'defined'!!!
detailed_protein_info = []
rejected_ribo_genes_count = []
accepted_ribo_genes_count = []
rejected_CDS_count = []
accepted_CDS_count = []
genome_GC = []
for IndexId in org_dat.IndexId:
    ###################
    gbk = gbdb[IndexId]
    fna = fndb[IndexId]
    ###################
    # in the case of Bacteria, we know for sure that a single accession number refers to related
    # genbank and nucleotide fasta, so there is no need to check if gbk.id matches fna.id, and beacuse 
    # all the accession numbers are associated with the complete genome entries (accessions starting with NC)
    # so, that there is no need to check if entry is a chromosome ...
    all_ribo_genes,rejected_ribo_genes = extract_ribo_genes(gbk,fna.seq)
    # we'll keep track of accepted and rejected ribosomal CDS-es, just to see the scale of the problem ...
    number_of_ribo_genes = len(all_ribo_genes)
    rejected_ribo_genes_count.append(rejected_ribo_genes)
    accepted_ribo_genes_count.append(number_of_ribo_genes)
    # regardless of the CAI business, we'd need to calculate genome-wide GC here ...
    GC = SeqUtils.GC(fna.seq)
    genome_GC.append(GC)
    # check number of ribo genes and avoid using fishy organisms with small amount of ribo genes, limit is set at 30 ...
    if number_of_ribo_genes < RIBO_LIMIT:
        # for EVERY IndexId report if there will be detailed info on the genes/proteins or NOT...
        detailed_protein_info.append(False)
        print >> sys.stdout, "Warning! There are %d<%d ribosomal proteins for %s IndexId; Skip it!"%(number_of_ribo_genes,RIBO_LIMIT,IndexId)
        all_ribo_genes = None
        # this assembly is skipeed, and no file is generated for it at all ...
        # EXPLICITLY skip this, and go to the next IndexId in the IndexId-loop (just for readability!)
        rejected_CDS_count.append(0)
        accepted_CDS_count.append(0)
        continue
    else:
        # for EVERY IndexId report if there WILL be detailed info on the genes/proteins or not...
        detailed_protein_info.append(True)
        # codons dictionary ...
        CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 
                      'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 
                      'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 
                      'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 
                      'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0, 
                      'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0, 
                      'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 
                      'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 
                      'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0, 
                      'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0, 
                      'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 
                      'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 
                      'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
        # count codon in the aal_ribo_genes ...
        # this will change CodonsDict
        count_codons(all_ribo_genes,CodonsDict)
        # generate codon index using the Codon count from CodonsDict ...
        CodonIndex = generate_codon_index(CodonsDict)
        # (2) 
        # this is alright, because there are ribosomal proteins at the very least!
        # it implies, all these loops wouldn't be empty and the files creation is justified!
        with open(os.path.join(path_CAI,'%s_genes.dat'%IndexId),'w') as fp:
            fp.write('prot_id,cai,gene_product,gene_seq,prot_seq\n')
            # in the case of Bacteria, we know for sure that a single accession number refers to related
            # genbank and nucleotide fasta, so there is no need to check if gbk.id matches fna.id, and beacuse 
            # all the accession numbers are associated with the complete genome entries (accessions starting with NC)
            # so, that there is no need to check if entry is a chromosome ...
            extracted_genes_features, rejected_genes_counter = extract_genes_features(gbk,fna.seq)
            accepted_genes_counter = len(extracted_genes_features)
            for prot_id,gene_product,gene_seq,prot_seq in extracted_genes_features:
                cai = cai_for_gene(gene_seq,CodonIndex)
                fp.write('%s,%.4f,%s,%s,%s\n'%(prot_id,cai,gene_product,gene_seq,prot_seq))
        rejected_CDS_count.append(rejected_genes_counter)
        accepted_CDS_count.append(accepted_genes_counter)
        #############################
        # print >> sys.stderr,'Overall genes stat: %d rejected out of %d that were detected'%(rejected_genes_counter, rejected_genes_counter+accepted_genes_counter)
        #############################






org_dat['protein_details'] = detailed_protein_info
org_dat['rejected_ribo_genes_count'] = rejected_ribo_genes_count
org_dat['accepted_ribo_genes_count'] = accepted_ribo_genes_count
org_dat['rejected_CDS_count'] = rejected_CDS_count
org_dat['accepted_CDS_count'] = accepted_CDS_count
org_dat['GC'] = genome_GC
org_dat.to_csv('catalog_with_accesion.dat',index=False)














