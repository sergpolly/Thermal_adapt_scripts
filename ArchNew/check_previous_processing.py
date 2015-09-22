import pandas as pd
import numpy as np
import re
# from fuzzywuzzy import fuzz
# from fuzzywuzzy import process



# archaea genomes description from /ftp...ncbi/genomes/genbank/archaea/...
asm = pd.read_csv("assembly_summary_arch_ncbi.txt",sep='\t')
# geomes description and connection with organism names/assembly id/taxonid established earlier ...
phyl = pd.read_csv("phylog_temps.txt")
# KZ database of manually retrieved OptimumTemperatures of different archaea ...
therm = pd.read_csv("termoprop.txt",sep='\t')


# the problem is that the taxonid@therm is not fully within the taxonid@asm (~150 out of ~250).

# let's compare by Assembly accession index ...

# we need a transformation:
# GCA_000144915.1_ASM14491v1 -> GCA_000007185.1
# because left one is what we have in phyl, while right one is the assembly accession (from asm) ...
asm_tmp = re.compile("(.+\_.+\.\d+)\_.+")
asm_check = lambda line: asm_tmp.match(line).groups()[0]
asm_if_check = lambda line: bool(asm_tmp.match(line))


# check if all asm index in phyl comply ...
print " Check if all comply ~GCA_xxxxxxxx.1_ASM14491v1: ", phyl.assembly.apply(asm_if_check).all()

# now get the accession form ...
asm_acc = phyl.assembly.apply(asm_check)

# and now check the intersection between asm_acc and accession from asm ...
asm_interection = np.intersect1d( asm_acc, asm.assembly_accession)


print "%d  out of %d from old list (phyl) are in the new GenBank list (asm)"%(asm_interection.size,asm_acc.shape[0])




# that's the end f verification by now.
# we would need to check organism names then and afterwards proceed with the analysis ...























































