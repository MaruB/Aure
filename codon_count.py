
from Bio import SeqIO # To parse a FASTA file
from sys import argv
from pandas import DataFrame
import pandas as pd
import logging

script, fasta_file_input = argv

# the cds from ENSEMBL also has the genes from Mt and Chl
# keep them or take out those within script or ask for an
# input that only has genes from chromosomes?

# don't know how to put options like -f input fasta etc

# raise an error if the input file is not .fasta
# is that necessary? put options to parse other kinds of files?
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


# set up logging of errors

# logger = logging.getLogger(__name__)  # recommended format

logger = logging.getLogger('codon_count')
logger.setLevel(logging.INFO)

# will print to the console warnings and above
console = logging.StreamHandler()
console.setLevel(logging.WARNING)

# will print to the log.txt file info and above
# writes over any existing log.txt file erasing previous contents
log_handler = logging.FileHandler('log.txt', mode = 'w')
log_handler.setLevel(logging.INFO)


# setting format for the logger
formatter = logging.Formatter('%(asctime)s - %(name)s - '
                              '%(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
console.setFormatter(formatter)

logger.addHandler(log_handler)
logger.addHandler(console)


def count_codons(fasta_file):

    """ Counts codons from a cds fasta file.
    Needs to start with ATG and (bla bla requirements)"""

    handle = open(fasta_file, 'r')
    logger.info('Reading FASTA file')
    # issue warning if fasta does not pass controls

    # initialize an empty list to store different codons dictionaries
    # one for each gene ID
    output = []

    # initialize empty list to store Gene IDs
    gene_ids = []

    # initialize variable with number of illegal codons
    illegal_codons = 0


    # iterate over sequence and count all the codons in the FastaFile
    for cur_record in SeqIO.parse(handle, 'fasta'):

        # needs to check if the format and contents are appropriate
        # eg: starts with ATG, ends with stop codon, is a multiple of 3,
        # only contains ATCG(no)

        # make the codon dictionary local
        codon_count = CodonsDict.copy()
        # change DNA sequences to uppercase
        dna_sequence = str(cur_record.seq).upper()  # do it with if or not necessary?

        # read from beginning to end, 3 by 3
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]

            if codon in codon_count:

                codon_count[codon] += 1

            else:
                illegal_codons += 1
                logger.info('ID %s contains an illegal codon', cur_record.id)

        # add gene ids and codon counts to their respective list
        gene_ids.append(cur_record.id)
        output.append(codon_count)
    # issue summary warning if there are illegal codons in the fasta file
    if illegal_codons == 0:
        logger.info('No illegal codons in the fasta file')
    else:
        logger.warning('There are %d illegal codons in the fasta file',
                       illegal_codons)

    # how to print the gene name with different fasta formats
    # gff?

    # create data frame from lists
    df_output = pd.DataFrame(output)
    df_output['gene_id'] = gene_ids
    cols = df_output.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_output = df_output[cols]
    # df_output.to_csv('output.csv')  # make the file name related to the input
    # for testing
    df_output.to_csv('output_test.csv') # make the file name related to the input
    # maybe create a folder for output?
    logger.info('Creating output file: "output.csv"')

    logger.warning('Check log.txt for more information')  # only if warnings?

    # close handle
    handle.close()


count_codons(fasta_file_input)


