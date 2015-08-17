from Bio import SeqIO


inventory = 'genbank.inventory'
db_file = 'genbank.idx'


# get the files names to be indexed ...
with open(inventory,'r') as fp:
    fnames = [line.strip() for line in fp.readlines()]

# index em!
res = SeqIO.index_db(db_file,filenames=fnames,format='genbank')

#quit ...
