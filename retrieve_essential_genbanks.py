import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import numpy as np
import pandas as pd

# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER/ftp.ncbi.nih.gov/refseq/release/bacteria/genbank')
path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
path_nucs = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/nucfasta')

print "Loading original databases and description file ..."

# summary file parsing ...
summary_fname = "env_catalog_compgenome.dat"
summary = pd.read_csv(os.path.join(path,summary_fname))

# genbank DB ...
gbdb_fname = os.path.join(path,'genbank.idx')
gbdb = SeqIO.index_db(gbdb_fname)

# nucleotides fasta files ...
fndb_fname = os.path.join(path_nucs,'nucfasta.idx')
fndb = SeqIO.index_db(fndb_fname,alphabet=IUPAC.IUPACAmbiguousDNA(),key_function=(lambda name: name.split('|')[3]))
# gbdb and fndb keys must be the same ...


print "Original databases are loaded ..."
print
print "Output the reduced genbank with the required entries only ..."

gb_small_fname = "tmp_gbank.gb"
with open(os.path.join(path,gb_small_fname),"w") as fp:
    for idx in summary['GenomicID']:
        fp.write(gbdb.get_raw(idx))


# index in-memory style the smaller subset of the genbank database ...
gbdb_small = SeqIO.index(os.path.join(path,gb_small_fname),"genbank")
# records = SeqIO.to_dict(SeqIO.parse("Quality/example.fastq", "fastq"))
# maybe to_dict would be faster, as we are parsing all of the records anyways ...
# maybe having another tmp file for nucfasta records would make things faster
# try it! sometime. So far, this takes ~40 mins
print "tmp file is written and indexed ..."
print
print "Matching nucleotide sequences from fasta with corresponding SeqRecords in genbank ..."

sequences = []
for idx in summary['GenomicID']:
    seqrec = gbdb_small[idx]
    nuc_seq = fndb[idx].seq
    # next thing should never happen in practice ...
    if ( len(seqrec.seq)!=len(nuc_seq) ):
        print "Genome length doesn't match between fasta and genbank for %s record"%idx
        print "Proceed anyways ..."
    # give the genbank, its actual nucleotide sequence ...
    seqrec.seq = nuc_seq
    sequences.append(seqrec)

print "All set"
print 
print "Writing final genbank ..."

# size_seqs = sys.getsizeof(sequences)
# print "All interesting sequences are in memory! %d byte of interesting sequences."%size_seqs

print "Writing all of that to the condensed.gb file"
with open("condensed.raw.gb","w") as fp:
    SeqIO.write(sequences,fp,"genbank")










