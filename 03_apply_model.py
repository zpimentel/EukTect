import os
import sys
import pickle

import numpy as np
import pandas as pd
from Bio import SeqIO


def is_valid_kmer(seq):
    ''' Used to remove any kmers with N's in them '''
    return set(seq).issubset({"A", "T", "C", "G"})


def is_valid_contig(seq):
    ''' Used to remove any kmers with N's in them '''
    return set(seq).issubset({"A", "T", "C", "G", "N"})


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)


if __name__ == "__main__":

    infile = sys.argv[1]
    if not os.path.exists(infile):
        raise Exception("Input file does not exist.")

    if is_fasta(infile) is False:
        raise Exception("Input file is not a valid fasta file.")

    # read in the model and list of canonical 5-mers
    clf = pickle.load(open("kmer_SVM.sav", 'rb'))
    kmer_list = pickle.load(open("kmer_list.pkl", 'rb'))

    # count canonical 5-mers from the user-defined fasta file
    kmer_dict = {}
    kmer_count = 0
    with open(infile, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if is_valid_contig(record.seq) is False:
                raise Exception("Fasta file contains non-A,T,G,C,N characters")
            if len(record.seq) < 2500:
                print(f"Warning: {record.id} is only {len(record.seq)} bp and thus the prediction may not be very accurate")
            if kmer_count < 1:
                kmer_count += 1
                for i in range(len(record.seq)):
                    kmer = record.seq[i:i+5].upper()
                    if str(kmer.reverse_complement()) in kmer_list:
                        kmer = str(kmer.reverse_complement())
                    if len(kmer) == 5:
                        if is_valid_kmer(str(kmer)) is True:
                            if str(kmer) in kmer_dict:
                                kmer_dict[str(kmer)] += 1
                            else:
                                kmer_dict[str(kmer)] = 1

    # make sure all kmers are accounted for in the dictionary
    for kmer in kmer_list:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = 0

    # convert to proportional by dividing each count by sum of kmer counts per contig
    kmer_df = pd.DataFrame.from_dict(kmer_dict, orient='index')
    kmer_df = kmer_df.reindex(sorted(kmer_list))
    kmer_df = kmer_df.div(kmer_df.sum(axis=0), axis=1)
    contig_kmers = kmer_df[0].values.tolist()
    contig_kmers = np.array(contig_kmers).reshape(1, -1)

    # perform prediction
    predictions = clf.predict(contig_kmers)
    if predictions == 1:
        print("This appears to be eukaryotic.")
    elif predictions == 0:
        print("This appears to be non-eukaryotic.")
    else:
        print("Sorry, something went wrong!")
