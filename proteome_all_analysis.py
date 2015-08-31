import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
from functools import partial
import time

from multiprocessing import Pool


aacids = list('CMFILVWYAGTSNQDEHRKP')

def get_aausage_proteome(seqrec):
    # seqrec = db[seqrec_id]
    features = seqrec.features
    proteome = []
    for feature in features:
        qualifiers = feature.qualifiers
        if (feature.type == 'CDS')and('translation' in qualifiers):
            proteome.append(qualifiers['translation'][0])
    #return the results ...
    proteome = ''.join(proteome)
    prot_len = float(len(proteome))
    aa_freq = tuple(proteome.count(aa)/prot_len for aa in aacids)
    #
    return (int(prot_len),) + aa_freq

def analyse_genome(db,seqrec_id):
    seqrec = db[seqrec_id]
    pl_aa_freq = get_aausage_proteome(seqrec)
    gc = SeqUtils.GC(seqrec.seq)
    id = seqrec.id
    return (id,gc) + pl_aa_freq


if __name__ == "__main__":
    path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')
    # dbx = SeqIO.index_db(os.path.join(path,'genbank.idx'))
    dbx = SeqIO.index_db(os.path.join(path,"subset.idx"))
    # subset.idx is a previously indexed db for the single genbank file with complete genomes only (the one ~5GB).
    dat = pd.read_csv(os.path.join(path,"env_catalog_compgenome.dat"))
    #
    #
    def do_work(seqrecid):
        return analyse_genome(dbx,seqrecid)
    #
    work = list(dat.GenomicID)
    #
    print "launching processes, to do %d pieces of work ..."%len(work)
    pool = Pool(processes=16)
    results = pool.map(do_work, work)
    #
    print "file outputting ..."
    with open("proteome_all.dat","w") as fp:
        fp.write("GenomicID,GC,ProtLen,"+",".join(aacids)+"\n")
        for result in results:
            fp.write(','.join(str(item) for item in result)+'\n')



# DbxRefs,Description,FeaturesNum,GenomicID,
# GenomicLen,GenomicName,Keywords,NucsPresent,
# Organism_des,SourceDbxRefs,SourceOrganism,
# SourcePlasmid,SourceStrain,Taxonomy,BioProject,
# TaxonID,Organism_env,OptimumTemperature,
# TemperatureRange,OxygenReq,Habitat,Salinity,
# crit_NC,crit_WGS,crit_genlen,crit_features,
# crit_comp_genome,crit_plasmid



# takes 4 min 40 seconds to accomplish for 4 processes ...
# takes 1 min 50 seconds to accomplish for 12 processes ...
# takes       48 seconds to accomplish for 32 processes ...
# everything ob bib6


