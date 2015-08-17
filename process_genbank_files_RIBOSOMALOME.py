
from Bio import Entrez as ncbi
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
data_path = os.path.join(root_path,"ftp.ncbi.nih.gov/genomes/all")
# asm_path = os.path.join(root_path,"archaea_links")
result_path = os.path.join("..","archaea250")


# # let ncbi know who you are ...
# ncbi.email = "sergey.venev@umassmed.edu"

# # read table with organisms and their taxonid ...
dat = pd.read_csv(os.path.join(root_path,"assembly_description.dat"))

aacids = sorted(list('CMFILVWYAGTSNQDEHRKP'))
# aacids = sorted(aacids)

# def get_cds_count(gb):
#     cds_features_counter = 0
#     for feature in gb.features:
#         if (feature.type == 'CDS')and('translation' in feature.qualifiers.keys()):
#             cds_features_counter += 1
#     #return the results ...
#     return cds_features_counter

def extract_proteins(gb):
    proteome = ''
    for feature in gb.features:
        if (feature.type == 'CDS')and('translation' in feature.qualifiers.keys()):
            proteome += feature.qualifiers['translation'][0]
    #return the results ...
    return proteome


def check_seq(seq):
    nucs = sum([int(nuc in seq) for nuc in list('ATGC')])
    return nucs == len('ATGC')


# let's simpy avoid plasmids ...
plasmid_pattern = re.compile('plasmid')
def check_chromosome(gb):
    return not plasmid_pattern.search(gb.description)


# get_ribosomalome ...
RIBO_LIMIT = 30
ribosomal_check = re.compile('ribosomal protein')
def get_ribosomalome(genbank):
    ribosomalome = []
    for feature in genbank.features:
        qual_keys = feature.qualifiers.keys()
        if (feature.type == 'CDS')and('translation' in qual_keys)and('product' in qual_keys):
            if ribosomal_check.search(feature.qualifiers['product'][0]):
                ribosomalome.append(feature.qualifiers['translation'][0])
    return ribosomalome

################################
# # ribo counter check ...
# if ribo_counter < RIBO_LIMIT:
#     print >> sys.stderr, "Warning! There are under %d ribosomal proteins in %s organism"%(RIBO_LIMIT,genbank.id)
#     return ''
# #return the results ...
# return ''.join(ribosomalome)
################################






# data we are interested in:
# (1)number of features in gb, (2)are there protein fasta file,
# (3)nucleotide sequence in gb, (4) total number of genbank entries!

# cds_features = []
# num_proteins = []
# nucseq_gb = []
# gb_entries = []
# accession = []

GC = []
# prot_dat = []
aafs = {}
for aa in aacids:
    aafs[aa] = []

genome_length = []
proteome_length = []


for asm in dat.assembly:
    # look up genome folder corresponding to the assembly...
    asm_folder_content = sub.check_output("ls %s|grep %s"%(os.path.join(data_path,asm),asm),shell=True)
    # check if there is a genbank file 
    if re.search("genomic.gbff",asm_folder_content):
        gb_fname = "%s/%s_genomic.gbff"%(os.path.join(data_path,asm),asm)
        genbank = SeqIO.parse(gb_fname,format='genbank')
        ribosomalome = []
        genome = ''
        for entry in genbank:
            if check_chromosome(entry):
                # adding up the list of ribosomal proteins ...
                ribosomalome += get_ribosomalome(entry)
                if check_seq(entry.seq):
                    genome += entry.seq
        ################################
        # finalizing ribosomalome stuff ...
        ################################
        if len(ribosomalome) < RIBO_LIMIT:
            print >> sys.stderr, "Warning! There are under %d ribosomal proteins in %s assembly"%(RIBO_LIMIT,asm)
            ribosomalome = ''
        else:
            ribosomalome = ''.join(ribosomalome)
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
        ribo_len = len(ribosomalome)
        # lengths ...
        genome_length.append(len(genome))
        proteome_length.append(ribo_len)
        if ribo_len:
            for aa in aacids:
                aafs[aa].append(ribosomalome.count(aa)/float(ribo_len))
        else:
            for aa in aacids:
                aafs[aa].append(0.0)
        #################################
    else:
        print "No genbank for assembly %s, Terminate!!!"%asm
        sys.exit(1)




# dat['cds_features'] = cds_features
# dat['num_proteins'] = num_proteins
# dat['nucseq_gb'] = nucseq_gb
# dat['gb_entries'] = gb_entries
# dat['accession'] = accession


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###########################
# Length of values does not match length of index
###########################
dat['GC'] = GC
for aa in aacids:
    dat[aa] = aafs[aa]

dat['genome_length'] = genome_length
dat['riboprot_length'] = proteome_length


dat.to_csv('archaea_ribosomalome.dat',index=False)


# Warning! There are under 30 ribosomal proteins in GCA_000513855.1_ASM51385v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000389735.1_Acd1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000247545.1_ASM24754v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000710615.1_HAPD43.1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000755225.1_ASM75522v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000739575.1_ASM73957v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000226975.2_HalLac_AJ5_1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000259215.1_ASM25921v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000755245.1_ASM75524v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000739555.1_ASM73955v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000379085.1_ASM37908v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000427685.1_ASM42768v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000739595.1_ASM73959v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000421805.1_ASM42180v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000710605.1_HTGA29.1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000690595.1_LPSSH13 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000447865.2_ASM44786v2 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000745485.1_ASM74548v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000746075.1_ASM74607v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000513315.1_ANOR1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000529525.1_Methanobrevibacter_oralis_JMR01 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000621965.1_ASM62196v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000320505.2_ASM32050v2 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000371805.1_ASM37180v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000376965.1_ASM37696v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000430905.1_ASM43090v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000711215.1_ASM71121v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000235685.3_ASM23568v3 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000744315.1_ASM74431v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000711905.1_ASM71190v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000308215.1_ASM30821v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000007185.1_ASM718v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000258515.1_TheZil1.0 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000585495.1_ASM58549v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000008085.1_ASM808v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000192595.1_ASM19259v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000192615.1_ASM19261v1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000685395.1_NKMY2 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000494145.1_Aigarchaeota_archaeon_OTU_1 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000494165.1_Aigarchaeota_archaeon_OTU_2 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000494185.1_Aigarchaeota_archaeon_OTU_3 assembly
# Warning! There are under 30 ribosomal proteins in GCA_000494125.1_Aigarchaeota_archaeon_OTU_4 assembly



# dat['prot_dat'] = prot_dat
# for aa in aacids:
#     dat[aa] = aafs[aa]


# dat_real = dat[dat.prot_dat=='yes']


# for aa in aacids:
#     plt.clf()
#     for i,cl in enumerate(dat_real.groupby(by='class')):
#         plt.plot(cl[1].topt,cl[1][aa],color=kelly_colors_hex[i],marker='o',label=cl[0],ms=9)
#     plt.title('amino acid=%s'%aa)
#     leg=plt.legend(loc='best',fontsize=9)
#     leg.get_frame().set_alpha(0)
#     leg.get_frame().set_edgecolor('white')
#     plt.savefig("aas_%s_class.pdf"%aa)





# for aa in aacids:
#     plt.clf()
#     i=1
#     # for i,cl in enumerate(dat_real.groupby(by='class')):
#     cl = ('Methanococci',dat_real[dat_real['class'] == 'Methanococci'])
#     plt.plot(cl[1].topt,cl[1][aa],color=kelly_colors_hex[i],marker='o',label=cl[0],ms=9)
#     plt.title('amino acid=%s'%aa)
#     leg=plt.legend(loc='best',fontsize=9)
#     leg.get_frame().set_alpha(0)
#     leg.get_frame().set_edgecolor('white')
#     plt.savefig("methanococci_aas_%s.pdf"%aa)





# for aa in aacids:
#     plt.clf()
#     # i=1
#     # for i,cl in enumerate(dat_real.groupby(by='class')):
#     ddd = dat_real[dat_real['class'] == 'Methanococci']
#     # cl = ('Methanococci',dat_real[dat_real['class'] == 'Methanococci'])
#     for i,cl in enumerate(ddd.groupby(by='genus')):
#         plt.plot(cl[1].topt,cl[1][aa],color=kelly_colors_hex[i],marker='o',label=cl[0],ms=9)
#     plt.title('amino acid=%s'%aa)
#     leg=plt.legend(loc='best',fontsize=9)
#     leg.get_frame().set_alpha(0)
#     leg.get_frame().set_edgecolor('white')
#     plt.savefig("methanococci_aas_%s.pdf"%aa)














