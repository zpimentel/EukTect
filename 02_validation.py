import os
import gzip
import pickle
import requests
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from sklearn.metrics import accuracy_score


def is_valid_sequence(seq):
    ''' Used to remove any kmers with N's in them '''
    return set(seq).issubset({"A", "T", "C", "G"})


if __name__ == "__main__":

    if not os.path.exists('./validation_genomes'):
        os.makedirs('./validation_genomes')

    # read in the refseq genome summary, SVM model, and list of kmers
    df = pd.read_csv("assembly_summary_refseq_labeled.csv")
    clf = pickle.load(open("kmer_SVM.sav", 'rb'))
    kmer_list = pickle.load(open("kmer_list.pkl", 'rb'))

    # Sample genomes to include in validation set
    sampled_euks = df[df['Classification'] == 'eukaryote'].sample(10, random_state=100)
    sampled_bact = df[df['Classification'] == 'bacteria'].sample(10, random_state=200)
    sampled_archae = df[df['Classification'] == 'archaea'].sample(10, random_state=200)
    sampled_virus = df[df['Classification'] == 'virus'].sample(10, random_state=200)

    # Prepare list of ftps addresses of genomes to be downloaded
    other_taxa_list = [sampled_bact, sampled_archae, sampled_virus]
    sampled_others = pd.concat(other_taxa_list)

    sampled_euks['# assembly_accession'] = sampled_euks['# assembly_accession'].str.split('.').str[0]
    sampled_others['# assembly_accession'] = sampled_others['# assembly_accession'].str.split('.').str[0]

    sampled_euks_list = sampled_euks['# assembly_accession'].to_list()
    sampled_other_list = sampled_others['# assembly_accession'].to_list()

    # Download the genomes
    for ftp in sampled_euks['ftp_path'].to_list():
        genome_id = ftp.split("/")[-1]
        gen_url = os.path.join(ftp, genome_id + "_genomic.fna.gz").replace(" ", "_")
        r = requests.get(gen_url, allow_redirects=True)
        open("validation_genomes/" + genome_id + ".fna.gzip", 'wb').write(r.content)

    for ftp in sampled_others['ftp_path'].to_list():
        genome_id = ftp.split("/")[-1]
        gen_url = os.path.join(ftp, genome_id + "_genomic.fna.gz").replace(" ", "_")
        r = requests.get(gen_url, allow_redirects=True)
        open("validation_genomes/" + genome_id + ".fna.gzip", 'wb').write(r.content)

    # count kmers in the validation genomes
    kmer_dict = defaultdict(lambda: 0)
    contigs = set(" ")
    directory = "./validation_genomes"
    euk_count = 0
    bact_count = 0
    for filename in os.listdir(directory):
        if filename.endswith(".fna.gzip"):
            genome_bp_count = 0

            genome_id = "_".join(filename.split("_")[0:2]).split('.')[0]

            if genome_id in sampled_euks_list: genome_type = 'euk'
            elif genome_id in sampled_other_list: genome_type = 'other'
            else: genome_type = 'NA'

            print(genome_id, genome_type)

            with gzip.open(os.path.join(directory, filename), "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if genome_bp_count < 50000:
                        if len(record.seq) >= 5000:
                            for i in range(0, len(record.seq), 5000):
                                genome_bp_count += 5000
                                if genome_type == 'euk':
                                    euk_count += 1
                                    contig_5kb_name = "euk_" + str(euk_count)
                                elif genome_type == 'other':
                                    bact_count += 1
                                    contig_5kb_name = "other_" + str(bact_count)
                                contigs.add(contig_5kb_name)

                                if genome_bp_count < 50000:
                                    contig_5kb = record.seq[i:i+5000]
                                    for j in range(len(contig_5kb)):
                                        kmer = contig_5kb[j:j+5].upper()
                                        if kmer.reverse_complement() in kmer_list:
                                            kmer = kmer.reverse_complement()
                                        if len(kmer) == 5:
                                            if is_valid_sequence(kmer) is True:
                                                kmer_dict[(contig_5kb_name, str(kmer))] += 1
                                else:
                                    break
                    else:
                        break

    # save the kmer counts for validation genomes
    f = open("kmer_matrix_validation.tsv", "w")
    f.write("Contig" + "\t")
    f.write("\t".join(contigs))
    f.write("\n")
    for kmer in kmer_list:
        f.write(kmer + "\t")
        for contig in contigs:
            f.write(str(kmer_dict[(contig, kmer)]) + "\t")
        f.write("\n")
    f.close()

    # drop the dictionary from memory just in case this causes issues - its large
    kmer_dict = []

    # read in counts of kmers per contig, remove contigs with not many kmer counts
    # convert to proportional by dividing each count by sum of kmer counts per contig
    kmer_df = pd.read_csv("kmer_matrix_validation.tsv", sep="\t", index_col=False)
    kmer_df = kmer_df.set_index("Contig")
    kmer_df = kmer_df.loc[:, kmer_df.sum(axis=0) > 1000]
    kmer_df = kmer_df.div(kmer_df.sum(axis=0), axis=1)
    kmer_df = kmer_df.reindex(sorted(kmer_list))

    # make list of lists containing kmer frequencies for training
    contig_kmers = []
    for column in kmer_df:
        kmer_list_contig = kmer_df[column].tolist()
        contig_kmers.append(kmer_list_contig)

    # make a list containing labels for training
    prelabels = kmer_df.columns.tolist()
    labels = []
    for i in prelabels:
        if i.split("_")[0] == "euk":
            labels.append(1)
        else:
            labels.append(0)

    # run preductions on validation genomes and check accuracy
    predictions = clf.predict(contig_kmers)
    print(f"Accuracy in the validation set is {accuracy_score(labels, predictions) * 100}%")
