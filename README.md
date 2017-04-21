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


THE MOST IMMEDIATE TODO:
The next big thing to finish up is the Archaeal part of the code.
Goto ~/Dropbox (UMASS MED - BIB)/protein_design/NEW_DATASETS_genbank and work out all the neccesary scripts to get Archaea to the same stage as Bacteria in terms of research.
most likely, we'd need to download Archaeal genbanks again (folders) ftp.ncbi.nih.gov/genomes/all ...
Next important file is `phylog_temps.txt`, we should try matching orgnisms from that file with the ones available at NCBI GenBank - that would be our to go point. Check how it was done before ?


###############################################################
# POST MBE REVISION NOTES
###############################################################

1) MBE major revision request to explicitly consider protein abundance in the model: some modeling required again.
2) All modeling-related soft is on COSMOS 'correct_galeprot_design_data/Argentina_MJ99/design'
3) Simulated frequencies are stored in files 'L64_ffMJ99_w0.060_T0.70_RND787202457-allfreq.dat'
4) 'data_processing_corrected_MJ99_Argentina.py' relies on the raw simulated data 
5) 'NEW_data_extraction_corrected_MJ99_Argentina.py' also relies on raw simulated data, it also generated 'simulation.dat'
6) 'simulation.dat' is used by many plotting scripts and various analyses downstream

TODO:
1) relaunch GPU protein design on COSMOS, like a copy of 'correct_galeprot_design_data/Argentina_MJ99'







