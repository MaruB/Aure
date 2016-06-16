
from Bio import SeqIO # To parse a FASTA file
from sys import argv
from pandas import DataFrame
import pandas as pd

script, fasta_file_input = argv

# don't know how to put options like -f input fasta etc

# raise an error if the input file is not .fasta
# do some checks to see that the fasta is ok
# write an output error file if necessary

CodonsDict = {'TTT':0, 'TTC':0, 'TTA':0, 'TTG':0, 'CTT':0, 
'CTC':0, 'CTA':0, 'CTG':0, 'ATT':0, 'ATC':0, 
'ATA':0, 'ATG':0, 'GTT':0, 'GTC':0, 'GTA':0, 
'GTG':0, 'TAT':0, 'TAC':0, 'TAA':0, 'TAG':0, 
'CAT':0, 'CAC':0, 'CAA':0, 'CAG':0, 'AAT':0, 
'AAC':0, 'AAA':0, 'AAG':0, 'GAT':0, 'GAC':0, 
'GAA':0, 'GAG':0, 'TCT':0, 'TCC':0, 'TCA':0, 
'TCG':0, 'CCT':0, 'CCC':0, 'CCA':0, 'CCG':0, 
'ACT':0, 'ACC':0, 'ACA':0, 'ACG':0, 'GCT':0, 
'GCC':0, 'GCA':0, 'GCG':0, 'TGT':0, 'TGC':0, 
'TGA':0, 'TGG':0, 'CGT':0, 'CGC':0, 'CGA':0, 
'CGG':0, 'AGT':0, 'AGC':0, 'AGA':0, 'AGG':0, 
'GGT':0, 'GGC':0, 'GGA':0, 'GGG':0}


def count_codons(fasta_file):

    """ Counts codons from a cds fasta file.
    Needs to start with ATG and (bla bla requirements)"""
    
    handle = open(fasta_file, 'r')


    # initialize an empty list to store different dictionaries


    output = []

    # initialize empty list to store Gene IDs
    gene_ids = []

    # iterate over sequence and count all the codons in the FastaFile.
    for cur_record in SeqIO.parse(handle, "fasta"):
        # make the codon dictionary local
        codon_count = CodonsDict.copy()
        dna_sequence = str(cur_record.seq).upper() # do it with if or not necessary?

        # read from beginning to end, 3 by 3
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            if codon in codon_count:
                codon_count[codon] += 1

            else: # like this or should write into a file too?
                raise TypeError("Illegal codon %s in gene: %s" % (codon, cur_record.id))

        gene_ids.append(cur_record.id)
        output.append(codon_count)


    df_output = pd.DataFrame(output)
    df_output['gene_id'] = gene_ids
    cols = df_output.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_output = df_output[cols]
    df_output.to_csv('output.csv')

    # close handle
    handle.close()

count_codons(fasta_file_input)


