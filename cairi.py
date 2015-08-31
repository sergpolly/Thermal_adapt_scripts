import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils import CodonUsage
import math
from Bio import Data
from Bio import pairwise2


# number of nucleotide per codon ...
CODON_LEN = 3
# Number of extracted CDS coding for ribosomal proteins, that we accept as sufficient... 
RIBO_LIMIT = 20
# matched amino acids threshold for cDNA recovery procedure ...
MATCHED_RATIO = 0.5

# Codon dictionary for accounting, just for convenience -> get a deepcopy of it whenever you need to reset
CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 
              'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 
              'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 
              'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 
              'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0, 
              'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0, 
              'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 
              'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 
              'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0, 
              'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0, 
              'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 
              'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 
              'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}




def six_frames_translation(seq, genetic_table=11):
    """
    six_frames_translation(seq, genetic_table=11)
        seq: presumably coding DNA sequence in Seq or string format.
        genetic_table: index of the genetic table in use (Bacterial/Archaeal/PlantPlasmid, 11 is the default)
    Function returns a dictionary with keys being all 6 possible in-frame cDNA sequences, and
    values being corresponding protein translations until first in-frame STOP. Beware cDNA are not
    truncated at corresponding STOP codons, so they can be much long than their "to_stop" translations.
    """
    anti = seq.reverse_complement()
    length = len(seq)
    frames = {}
    for fid in range(CODON_LEN):
        fragment_length = CODON_LEN * ( (length-fid)//CODON_LEN )
        frames[str(seq[fid:fid+fragment_length])] = str(seq[fid:fid+fragment_length].translate(table=genetic_table,to_stop=True))
        frames[str(anti[fid:fid+fragment_length])] = str(anti[fid:fid+fragment_length].translate(table=genetic_table,to_stop=True))
    return frames



def get_putative_cDNA(ref_prot,frames):
    """
    get_putative_cDNA(frames)
        ref_prot: reference protein sequence with unknown cDNA (cDNA fuzzy location, etc.)
        frames: dictionary of 6*cDNA:translation key:value pairs, returned by six_frames_translation
    Takes 6-frames dictionary returned by six_frames_translation and looks for the cDNA best matching
    the reference protein. Returns a pair (best_cDNA, percent_matched), where percent_matched is
    the number of amino acids that matched to the reference protein sequence.
    """
    # presumbaly_top_align = pairwise2.align.globalxx(seq1,seq2)[0]
    # s1_aln, s2_aln, score, begin, end = presumbaly_top_align
    frame_align = ( (cDNA, pairwise2.align.globalxx(str(ref_prot),translation)[0][2] ) for cDNA,translation in frames.iteritems() if translation)
    # frames' translations aligned to the reference protein: iter of (cDNA, SCORE_aln)-pairs is returned ...
    # returning cDNA with the best alignment score ...
    best_cDNA, prot_score = max(frame_align,key=lambda x: x[1])
    #
    print prot_score, len(ref_prot)
    #
    # prot_score is simply the number of matches, in case we use the "dumb" globalxx aligner
    # see BioPython pairwise2 description for details http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
    # let's return the fraction of amino acids that we were able to align,
    # so that we could do an advised decision on whether that recovered cDNA is usable or not ...
    return (best_cDNA, prot_score/len(ref_prot))




###################################################################################
# CAI related functions ...
ribosomal_check = re.compile('ribosomal protein')
def extract_ribo_genes(genbank,nuc_seq):
    """ 
    extract_ribo_genes(genbank,nuc_seq)
        genbank: SeqRecord object with the annotated genome
        nuc_seq: Seq object with the nucleotide sequence of the genome in 'genbank'
    It goes through the genbank's list of features and extracts genes coding for ribosomal proteins:
    their feature id (in a given genbank), nucleotide and protein sequences and finally their acceptance status. 
    """
    data = {"fid":[],"nucs":[],"protein":[],"table":[],"status":[]}
    for fid,feature in enumerate(genbank.features):
        feat_quals = feature.qualifiers
        qual_keys = feat_quals.keys()
        if (feature.type == 'CDS')and('translation' in qual_keys)and('product' in qual_keys):
            if ribosomal_check.search(feat_quals['product'][0]):
                # this is the SUPER STRICT check for suitable ribo-genes ... 
                # extracting coding DNA sequences ...
                extracted_nucs = feature.extract(nuc_seq)
                # extracting available protein sequence ...
                extracted_protein = feat_quals['translation'][0]
                # try to translate it back using some genetic code, which one ?
                genetic_table = feat_quals['transl_table'][0] if ('transl_table' in qual_keys) else 11
                # genetic_table = 11 is the Bacterial, Archaeal and Plant Plastid Code, so it'll be our default, just in case ...
                #
                # fill in the dict, no matter what:
                data['fid'].append(fid)
                data['nucs'].append(extracted_nucs)
                data['protein'].append(extracted_protein)
                data['table'].append(genetic_table)
                #
                try:
                    translation = extracted_nucs.translate(table=genetic_table,cds=True)
                # catch the translational error exception ...
                except Data.CodonTable.TranslationError:
                    # if gene is not true-CDS and it cannot be translated we'll automatically reject it!
                    data['status'].append('reject')
                    # and move on to the next gene in the loop ...
                    continue
                # even if it is a true-CDS, we'll check if it translates back to what is in the feature qualifiers ..
                if str(translation)==feat_quals['translation'][0]:
                    # OK! extracted gene translates to the the provided feature qualifier ...
                    data['status'].append('accept')
                else:
                    # rejecting the gene if it does not match with the feature qualifier translation ...
                    data['status'].append('reject')
    return data



# this returns a dictionary of with assorted info on the CDS-es from a given genbank ...
def extract_genes_features(genbank,nuc_seq):
    """
    extract_genes_features(genbank,nuc_seq)
        genbank: SeqRecord object with the annotated genome
        nuc_seq: Seq object with the nucleotide sequence of the genome in 'genbank'
    It goes through the genbank's list of features and extracts CDS-es with translation and product description.
    their feature id (in a given genbank), nucleotide and protein sequences and finally their acceptance status. 
    """
    data = {"fid":[],"pid":[],"nucs":[],"protein":[],"cdna":[],"product":[],"table":[],"status":[]}
    #
    for fid,feature in enumerate(genbank.features):
        feat_quals = feature.qualifiers 
        qual_keys = feat_quals.keys()
        if (feature.type == 'CDS')and('translation' in qual_keys)and('product' in qual_keys):
            # we'll be extracting only nice CDS-es, ones that can be translated back to the translation present in features ...
            # extract coding nucleotide sequence:
            extracted_nucs = feature.extract(nuc_seq)
            # extracting available protein sequence ...
            extracted_protein = feat_quals['translation'][0]
            # try to translate it back using some genetic code, which one ?
            genetic_table = feat_quals['transl_table'][0] if ('transl_table' in qual_keys) else 11
            # genetic_table = 11 is the Bacterial, Archaeal and Plant Plastid Code, so it'll be our default, just in case ...
            #
            # fill in the dict, no matter what:
            data['fid'].append(fid)
            data['pid'].append(feat_quals['protein_id'][0])
            data['nucs'].append(extracted_nucs)
            data['protein'].append(extracted_protein)
            data['product'].append(feat_quals['product'][0].replace(',',' '))
            data['table'].append(genetic_table)
            #
            # trying clean CDS translation first ...
            try:
                translation = extracted_nucs.translate(table=genetic_table,cds=True)
            # we'll be catching TranslationError exceptions here ...
            except Data.CodonTable.TranslationError:
                # if extracted_nucs is not true-CDS, try manual translation thing ...
                translation = extracted_nucs.translate(table=genetic_table).strip('*')
                # we are not catching anything here, it's kinda dangerous though ...
            finally:
                # best situation, translation matches exactly ...
                if (str(translation) == extracted_protein):
                    # OK! extracted extracted_nucs translates to the the provided feature qualifier ...
                    data['status'].append("accepted")
                    data['cdna'].append('')
                # translations are matched >= 80% ...
                # examples: TGA is a stop codon in 11's table, yet it can code for U (Selenocysteine)
                elif pairwise2.align.globalxx(str(translation),extracted_protein)[0][2] >= 0.8*len(extracted_protein): 
                    # OK! extracted extracted_nucs translates to the the provided feature qualifier
                    # ALMOST exactly ...
                    data['status'].append("accepted")
                    data['cdna'].append('')
                # try something else ...
                else:
                    #LAST RESORT ...
                    # before rejecting cDNA completely, let's try recovering cDNA manually ...
                    print "trying to recover cDNA @ genbank %s feature %d ..."%(genbank.id,int(fid))
                    extracted_cDNA_trans = six_frames_translation(extracted_nucs, genetic_table=genetic_table)
                    putative_cDNA, matched_aa = get_putative_cDNA(extracted_protein,extracted_cDNA_trans)
                    # compare matched_aa with MATCHED_RATIO:
                    if (matched_aa >= MATCHED_RATIO):
                        data['status'].append("recovered")
                        data['cdna'].append(putative_cDNA)
                    else:
                        # final rejection ...
                        data['status'].append('rejected')
                        data['cdna'].append('')
            # finally clause over ...
        # if CDS,translation,product clause over ...
    # features loop over ...
    return data


#
def count_codons(dna_seq_list,cod_dict):
    """
    count_codons(dna_seq_list,cod_dict)
        dna_seq_list: list(or other iterable) of DNA sequences(Seq object or strings)
        cod_dict: dictionary with the valid codons as the keys.
    It goes over all sequences and count appearances of all codons.
    Sequences whose length is not proportional to 3(nucleotides per codon) are skipped!
    Codons containing anything but 'ATGC' are skipped as well.
    All sequences must be in frame, so that first 3 letters must constitute an actual codon (user responsibility).
    """
    # cod_dict must be zeroed out beforehand, that's a potential issue ...
    # codons = [dna_seq[i:i+CODON_LEN] for i in range(0,len(dna_seq),CODON_LEN)]
    for gene_sequence in dna_seq_list:
        gene_sequence = str(gene_sequence)
        len_gene_sequence = len(gene_sequence)
        if len_gene_sequence%CODON_LEN:
            print >> sys.stderr,"(Counting codons in count_codons) gene lengths is not ~3! Skipping it!!!"
        else:
            for codon in [gene_sequence[i:i+CODON_LEN] for i in range(0,len_gene_sequence,CODON_LEN)]:
                if codon in cod_dict.keys(): 
                    cod_dict[codon] += 1 
                else: 
                    # raise TypeError("Illegal codon %s in genes" % codon)
                    print >> sys.stderr, "(Counting codons in count_codons) Illegal codon %s in genes! Skipping it!" % codon


# CodonUsage.SynonymousCodons
def generate_codon_index(cod_dict,genetic_table=11):
    """
    generate_codon_index(cod_dict,genetic_table=11)
        cod_dict: dictionary with codons as keys and number of corresponding codon appearance as the value
        genetic_table: index of the genetic table in use (Bacterial/Archaeal/PlantPlasmid, 11 is the default)
    Function generates codon index, based on the input codon usage data (presumably from ribosomal proteins,
    or other highly expressed genes). One can modify the genetic table in use.
    """
    index = {}
    # Let's create the index dictionary for the specific genetic table provided ...
    # unambiguous DNA deals with 61 codons made of 
    genetic_code = Data.CodonTable.unambiguous_dna_by_id[genetic_table]
    SynonymousCodons = dict([(aa,[]) for aa in genetic_code.protein_alphabet.letters])
    # SynonymousCodons['STOP'] = genetic_code.stop_codons # STOP codons are excluded from analysis ...
    for codon,aa in genetic_code.forward_table.iteritems():
        SynonymousCodons[aa].append(codon)
    # now exclude amino acids coded by a single codons: M,W for standard genetic code ...
    unambiguous_aacids = [aa for aa,codons in SynonymousCodons.iteritems() if len(codons)<2]
    for aa in unambiguous_aacids:
        # if there is no ambiguity, there is no optimality ...
        _ = SynonymousCodons.pop(aa)
    #
    # now to calculate the index we first need to sum the number of times 
    # synonymous codons were used all together. 
    # here is the place for potential improvement, because CodonUsage.SynonymousCodons implies the standard genetic code,
    # while table=4 gives up 1 STOP codon in favor of a W(triptophane) one! While, table=11 uses the standart code, with the 
    # expanded range of start codons ... 
    for aa in SynonymousCodons: 
        SynCodons = SynonymousCodons[aa]
        ###################################
        # total number of codons coding for a given amino acid in the cod_dict ...
        total = sum([ cod_dict[scodon] for scodon in SynCodons ])
        if not total:
            print >> sys.stderr, "ACHTUNG!!! this is unprecedented!!!"
            print >> sys.stderr, "Codons coding for %s amino acid are not encountered in ribo-genes of this organism!" % aa
            print >> sys.stderr, "Codons corresponding to %s will appear as Illegal for this organism!" % aa            
            continue
        else:
            ###################################
            # calculate the RSCU value for each of the codons
            # RSCU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons) 
            # Relative Synonymous Codon Usage (RSCU)...
            # (sum of all synonymous codons/num of synonymous codons) is expected number of this-type-of-codon assuming no codon bias!
            no_bias_usage = float(total)/len(SynCodons)
            RSCU = [ cod_dict[scodon]/no_bias_usage for scodon in SynCodons ]
            ###################################
            # now generate the index W=RCSUi/RCSUmax: 
            RSCU_max = max(RSCU)
            for scodon,RSCUi in zip(SynCodons,RSCU):
                index[scodon] = RSCUi/RSCU_max
            ####################################
            # We have excluded stop codons and those that are singletons, but mind that:
            # index is alway 1, when amino acid is coded by a single codon ...
            # for trans_table 11 that would be 'ATG', 'TGG' coding for Met and W(Trypt.) respectively 
    return index



def cai_for_gene(dna_seq,index):
    """
    cai_for_gene(dna_seq,index)
    dna_seq: Seq object or Python string with the in-frame coding sequence
    index: Codon usage bias calculated fro presumably highly expressed CDSs, excluding singletons and STOP codons
    This function calculate CAI for the provided coding sequence.
    """
    dna_seq = str(dna_seq) 
    cai_value, cai_length, illegal_codons_count = 0.0, 0, 0
    for codon in [dna_seq[i:i+CODON_LEN] for i in range(0,len(dna_seq),CODON_LEN)]:
        if codon in index: 
            # STOP codons and singletons are already excluded ...
            if index[codon]>0.0:
                cai_value += math.log(index[codon]) 
                cai_length += 1
        # illegal codons are those, that contain something but 'ATGC' ...
        elif sum( [ (nuc in 'ATGC') for nuc in codon ] )<3:
            illegal_codons_count += 1
            print >> sys.stderr, "skipping illegal codon in sequence: %s" % codon
            # this line is from BioPython:
            # raise TypeError("illegal codon in sequence: %s.\n%s" % (codon, index))
            # we'll be skipping illegal codons instead ...
    if illegal_codons_count:
        print >> sys.stderr, "%d illegal codons were skipped in the %d nucs gene! %s" % (illegal_codons_count,len(dna_seq),dna_seq)        
    return math.exp(cai_value / (cai_length - 1.0)) 
###################################################################################

















