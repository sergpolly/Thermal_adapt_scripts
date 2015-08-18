import sys
from Bio import SeqIO

# see corresponding description in the project's wiki

inventory = sys.argv[1]
db_file = sys.argv[2]
seq_format = sys.argv[3]

# sanity check is ommited for such a short and simple script


# get the files names to be indexed
with open(inventory,'r') as fp:
    fnames = [line.strip() for line in fp.readlines()]

# index them depending on the format
if format == "genbank":
    res = SeqIO.index_db(db_file,filenames=fnames,format=seq_format)
elif format == "fasta":
    get_index = lambda name: name.split('|')[3]
    res = SeqIO.index_db(db_file,filenames=fnames,format=seq_format,key_function=get_index)
else:
    print "Only genbank and fasta formats are accepted!"
    sys.exit(1)


