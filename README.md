# Thermal_adapt_scripts
Scripts to pull and analyze thermal adaptation data for prokaryotes.


some immediate information on the data in ftp.ncbi.nih.gov/genomes/all:

Positive features:
1) All assemblies I dealt with, do have a genbank file inside the assembly folder (‘asm_genomic.gbff’)
2) Some, but not all of the assemblies, carry asm_protein.faa files with all the proteins in fast format
3) Number of proteins in .faa file does equal EXACTLY to the number of CDS features of corresponding genbank files. All of these CDS features contains qualifier ‘translation’, i.e. there is a protein sequence available for them

Negative features:
1) Not all of the genbank files carry nucleotide sequences (regardless of the protein annotations available) -> unpredictability makes it a negative feature.
2) Some assemblies carry suspiciously small number of proteins! e.g. 59,124 -> that’s too small for the entire genome. Thus there is a need for the cutoff.

Just remarks:
1) It is very likely that the order of entries in genbank and in nucleotide fasta is exactly the same. However I’d rather use dictionary to cross-access these entries!