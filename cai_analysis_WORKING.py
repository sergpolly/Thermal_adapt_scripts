import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez as ncbi
from Bio import SeqUtils
from Bio.SeqUtils import CodonUsage
import numpy as np
import pandas as pd
import math
from Bio import Data

import cairi
# extract_genes_features(genbank,nuc_seq)

from multiprocessing import Pool



# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
# # dbx = SeqIO.index_db(os.path.join(path,'genbank.idx'))
# dbx = SeqIO.index_db(os.path.join(path,"subset.idx"))
# # subset.idx is a previously indexed db for the single genbank file with complete genomes only (the one ~5GB).
# dat = pd.read_csv(os.path.join(path,"env_catalog_compgenome.dat"))
# #
# #
# def do_work(seqrecid):
#     seqrec = dbx[seqrecid]
#     return cairi.extract_genes_features(seqrec,seqrec.seq)
# #
# work = list(dat.GenomicID)
# #
# print "launching process, to do %d pieces of work ..."%len(work)
# results = list(do_work(piece) for piece in work[:10])
# # #
# # print "file outputting ..."
# # with open("proteome_all.dat","w") as fp:
# #     fp.write("GenomicID,GC,ProtLen,"+",".join(aacids)+"\n")
# #     for result in results:
# #         fp.write(','.join(str(item) for item in result)+'\n')








if __name__ == "__main__":
    path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
    # dbx = SeqIO.index_db(os.path.join(path,'genbank.idx'))
    dbx = SeqIO.index_db(os.path.join(path,"subset.idx"))
    # subset.idx is a previously indexed db for the single genbank file with complete genomes only (the one ~5GB).
    dat = pd.read_csv(os.path.join(path,"env_catalog_compgenome.dat"))
    # #
    # outpath = os.path.join(path,'parout')
    # out = lambda x: "out_%s.dat"%x
    # #
    def do_work(seqrecid):
        seqrec = dbx[seqrecid]
        return cairi.extract_genes_features(seqrec,seqrec.seq)
        # pd.DataFrame(dat).to_csv(os.path.join(outpath,out(seqrecid)),index=False)
    #
    work = list(dat.GenomicID)
    #
    processes = int(sys.argv[1])
    pool = Pool(processes=processes)
    print "launching %d processes, to do %d pieces of work ..."%(processes,len(work))
    results = pool.map(do_work, work)
    #
    print "parallel calculations are OVER!!!"
    # #
    i = 0 
    for res in results:
    	print work[i], len(res['fid'])
    	i += 1
    # print "file outputting ..."
    # with open("proteome_all.dat","w") as fp:
    #     fp.write("GenomicID,GC,ProtLen,"+",".join(aacids)+"\n")
    #     for result in results:
    #         fp.write(','.join(str(item) for item in result)+'\n')



# THIS ARE THE GENBANKS THAT GOT STUCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ['NC_006526.2', 'NC_010645.1', 'NC_015436.1']



# sr1 = dbx['NC_006526.2']
# sr2 = dbx['NC_010645.1']
# sr3 = dbx['NC_015436.1']

# fid = 11013

# trying to recover cDNA @ genbank NC_003888.3 feature 12136 ...

# # ALL organisms appear to be unique in the catalog, so that there are NO double-chromosome bacteria in our dataset...
# # that is strange, because I remember that I did encounter them somehow and that there are abnormally low amount of ribo-genes there ...
# # so anyways, even if I'm missing something, hopefully this two effects would cancel each other ...

# # next, we'll simply follow the same approach that we pursued with Archaea to extract genes and their CAIs ...
# # In the Archeal case, we used assembly code as a unique Organism identifier, whereas in the case of Bacteria, 
# # we limited ourselfes to the analysis of the complete genomes only (omitting WholeGenomeShotguns, and multi-chromosomal Organisms), 
# # so, in Bacterial case, it is the genome accession number that is the unique Orgnismal identifier. (accession number = IndexId in our case ...)
# #########################################################
# #########################################################
# #(1) for each IndexId - find ribosomalome genes and get the codon usage of highly expressed genes ...
# #(2) for each IndexId calculate CAI of every coding gene (CDS) in the genome
# #(3) output detailed info on protein/prot_id/prot_seq/CAI/etc for each IndexId in a separate file ...
# #(4) done!!!
# #########################################################
# #########################################################
# #########################################################
# # remark: could not find sufficient number of ribosomal genes -> no CAI!
# # meaning that CAI is only 'defined', in case Ribosomalome is 'defined'!!!
# detailed_protein_info = []
# rejected_ribo_genes_count = []
# accepted_ribo_genes_count = []
# rejected_CDS_count = []
# accepted_CDS_count = []
# genome_GC = []
# for IndexId in org_dat.IndexId:
#     ###################
#     gbk = gbdb[IndexId]
#     fna = fndb[IndexId]
#     ###################
#     # in the case of Bacteria, we know for sure that a single accession number refers to related
#     # genbank and nucleotide fasta, so there is no need to check if gbk.id matches fna.id, and beacuse 
#     # all the accession numbers are associated with the complete genome entries (accessions starting with NC)
#     # so, that there is no need to check if entry is a chromosome ...
#     all_ribo_genes,rejected_ribo_genes = extract_ribo_genes(gbk,fna.seq)
#     # we'll keep track of accepted and rejected ribosomal CDS-es, just to see the scale of the problem ...
#     number_of_ribo_genes = len(all_ribo_genes)
#     rejected_ribo_genes_count.append(rejected_ribo_genes)
#     accepted_ribo_genes_count.append(number_of_ribo_genes)
#     # regardless of the CAI business, we'd need to calculate genome-wide GC here ...
#     GC = SeqUtils.GC(fna.seq)
#     genome_GC.append(GC)
#     # check number of ribo genes and avoid using fishy organisms with small amount of ribo genes, limit is set at 30 ...
#     if number_of_ribo_genes < RIBO_LIMIT:
#         # for EVERY IndexId report if there will be detailed info on the genes/proteins or NOT...
#         detailed_protein_info.append(False)
#         print >> sys.stdout, "Warning! There are %d<%d ribosomal proteins for %s IndexId; Skip it!"%(number_of_ribo_genes,RIBO_LIMIT,IndexId)
#         all_ribo_genes = None
#         # this assembly is skipeed, and no file is generated for it at all ...
#         # EXPLICITLY skip this, and go to the next IndexId in the IndexId-loop (just for readability!)
#         rejected_CDS_count.append(0)
#         accepted_CDS_count.append(0)
#         continue
#     else:
#         # for EVERY IndexId report if there WILL be detailed info on the genes/proteins or not...
#         detailed_protein_info.append(True)
#         # codons dictionary ...
#         CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 
#                       'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 
#                       'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 
#                       'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 
#                       'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0, 
#                       'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0, 
#                       'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 
#                       'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 
#                       'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0, 
#                       'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0, 
#                       'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 
#                       'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 
#                       'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
#         # count codon in the aal_ribo_genes ...
#         # this will change CodonsDict
#         count_codons(all_ribo_genes,CodonsDict)
#         # generate codon index using the Codon count from CodonsDict ...
#         CodonIndex = generate_codon_index(CodonsDict)
#         # (2) 
#         # this is alright, because there are ribosomal proteins at the very least!
#         # it implies, all these loops wouldn't be empty and the files creation is justified!
#         with open(os.path.join(path_CAI,'%s_genes.dat'%IndexId),'w') as fp:
#             fp.write('prot_id,cai,gene_product,gene_seq,prot_seq\n')
#             # in the case of Bacteria, we know for sure that a single accession number refers to related
#             # genbank and nucleotide fasta, so there is no need to check if gbk.id matches fna.id, and beacuse 
#             # all the accession numbers are associated with the complete genome entries (accessions starting with NC)
#             # so, that there is no need to check if entry is a chromosome ...
#             extracted_genes_features, rejected_genes_counter = extract_genes_features(gbk,fna.seq)
#             accepted_genes_counter = len(extracted_genes_features)
#             for prot_id,gene_product,gene_seq,prot_seq in extracted_genes_features:
#                 cai = cai_for_gene(gene_seq,CodonIndex)
#                 fp.write('%s,%.4f,%s,%s,%s\n'%(prot_id,cai,gene_product,gene_seq,prot_seq))
#         rejected_CDS_count.append(rejected_genes_counter)
#         accepted_CDS_count.append(accepted_genes_counter)
#         #############################
#         # print >> sys.stderr,'Overall genes stat: %d rejected out of %d that were detected'%(rejected_genes_counter, rejected_genes_counter+accepted_genes_counter)
#         #############################






# org_dat['protein_details'] = detailed_protein_info
# org_dat['rejected_ribo_genes_count'] = rejected_ribo_genes_count
# org_dat['accepted_ribo_genes_count'] = accepted_ribo_genes_count
# org_dat['rejected_CDS_count'] = rejected_CDS_count
# org_dat['accepted_CDS_count'] = accepted_CDS_count
# org_dat['GC'] = genome_GC
# org_dat.to_csv('catalog_with_accesion.dat',index=False)



