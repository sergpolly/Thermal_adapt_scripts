
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
dat['proteome_length'] = proteome_length


dat.to_csv('archaea_repeat.dat',index=False)



# intermediate conclusion about the data:
# (dat.cds_features ==dat.num_proteins).all() -> True

# (dat.nucseq_gb <= dat.gb_entries).all() -> True
# (dat.nucseq_gb == dat.gb_entries).all() -> False!!!




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














