import os
import gzip
import pickle
import requests
import subprocess
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from sklearn import svm
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split


def download_file(url, outfile):
    r = requests.get(url, allow_redirects=True)
    open(outfile, 'wb').write(r.content)


def is_valid_sequence(seq):
    ''' Used to remove any kmers with N's in them '''
    return set(seq).issubset({"A", "T", "C", "G"})


def unzip_db(filename, outdir):
    cmd = ['unzip', filename, '-d', outdir]
    unzip_call = subprocess.call(cmd, shell=False, stderr=subprocess.STDOUT)


def assign_genomes(taxid, node_dict):
    """
    Given a NCBI taxonomy ID, this function figures out the domain that ID belongs to.
    """
    eukaryote = 2759
    bacteria = 2
    archaea = 2157
    virus = 10239

    taxa_list = [eukaryote, bacteria, archaea, virus]
    classification = "not classified"

    while taxid not in taxa_list:
        try:
            taxid = int(node_dict[str(taxid)])
            return(assign_genomes(taxid, node_dict))  # get the parent taxid and try again
        except Exception:
            return(classification)

    if taxid == eukaryote: classification = 'eukaryote'
    if taxid == bacteria: classification = 'bacteria'
    if taxid == archaea: classification = 'archaea'
    if taxid == virus: classification = 'virus'

    return(classification)


if __name__ == "__main__":

    if not os.path.exists('./genomes'):
        os.makedirs('./genomes')

    # Download files and prepare data
    taxdmp = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip'
    download_file(taxdmp, 'taxdmp.zip')
    unzip_db('taxdmp.zip', 'taxonomy')
    df = pd.read_csv("assembly_summary_refseq.txt", skiprows=1, sep="\t")

    # Parse nodes.dmp from NCBI taxonomy (map each taxid to its parent taxid)
    nodes_file = open("./taxonomy/nodes.dmp")
    node_dict = {}
    for row in nodes_file:
        row = row.strip().split('\t|\t')
        child = row[0]
        parent = row[1]
        node_dict[child] = parent

    # Figure out whether each genome is euk, bacteria, archaea, or virus
    df['Classification'] = df['taxid'].apply(assign_genomes, args=(node_dict,))
    df.to_csv("assembly_summary_refseq_labeled.csv")

    # Sample genomes to include in training and testing set
    sampled_euks = df[df['Classification'] == 'eukaryote'].sample(50, random_state=234)
    sampled_bact = df[df['Classification'] == 'bacteria'].sample(15, random_state=345)
    sampled_archae = df[df['Classification'] == 'archaea'].sample(15, random_state=345)
    sampled_virus = df[df['Classification'] == 'virus'].sample(20, random_state=345)

    # Prepare list of ftps addresses of genomes to be downloaded
    other_taxa_list = [sampled_bact, sampled_archae, sampled_virus]
    sampled_others = pd.concat(other_taxa_list)

    sampled_euks[sampled_euks.columns[0]] = sampled_euks[sampled_euks.columns[0]].str.split('.').str[0]
    sampled_others[sampled_others.columns[0]] = sampled_others[sampled_others.columns[0]].str.split('.').str[0]

    sampled_euks_list = sampled_euks[sampled_euks.columns[0]].to_list()
    sampled_other_list = sampled_others[sampled_others.columns[0]].to_list()

    # Download the genomes
    for ftp in sampled_euks['ftp_path'].to_list():
        genome_id = ftp.split("/")[-1]
        gen_url = os.path.join(ftp, genome_id + "_genomic.fna.gz").replace(" ", "_")
        r = requests.get(gen_url, allow_redirects=True)
        open("genomes/" + genome_id + ".fna.gzip", 'wb').write(r.content)

    for ftp in sampled_others['ftp_path'].to_list():
        genome_id = ftp.split("/")[-1]
        gen_url = os.path.join(ftp, genome_id + "_genomic.fna.gz").replace(" ", "_")
        r = requests.get(gen_url, allow_redirects=True)
        open("genomes/" + genome_id + ".fna.gzip", 'wb').write(r.content)

    # Compute counts of canonical kmers (word size 5) from chunks of 5kb from a total of 50kb from each genome
    kmers_visited = set()
    kmer_dict = defaultdict(lambda: 0)
    contigs = set(" ")
    directory = "./genomes"
    euk_count = 0
    bact_count = 0
    for filename in os.listdir(directory):
        if filename.endswith(".fna.gzip"):
            genome_bp_count = 0

            genome_id = "_".join(filename.split("_")[0:2]).split('.')[0]

            if genome_id in sampled_euks_list: genome_type = 'euk'
            elif genome_id in sampled_other_list: genome_type = 'other'
            else: genome_type = 'NA'

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
                                        if kmer.reverse_complement() in kmers_visited:
                                            kmer = kmer.reverse_complement()
                                        if len(kmer) == 5:
                                            if is_valid_sequence(kmer) is True:
                                                kmers_visited.add(str(kmer))
                                                kmer_dict[(contig_5kb_name, str(kmer))] += 1
                                else:
                                    break
                    else:
                        break

    # write the counts of kmers per contig to file
    f = open("kmer_matrix2.tsv", "w")
    f.write("Contig" + "\t")
    f.write("\t".join(contigs))
    f.write("\n")
    for kmer in kmers_visited:
        f.write(kmer + "\t")
        for contig in contigs:
            f.write(str(kmer_dict[(contig, kmer)]) + "\t")
        f.write("\n")
    f.close()

    # drop the dictionary from memory just in case this causes issues - its large
    kmer_dict = []

    # read in counts of kmers per contig, remove contigs with not many kmer counts
    # convert to proportional by dividing each count by sum of kmer counts per contig
    kmer_df = pd.read_csv("kmer_matrix2.tsv", sep="\t", index_col=False)
    kmer_df = kmer_df.set_index("Contig")
    kmer_df = kmer_df.loc[:, kmer_df.sum(axis=0) > 1000]
    kmer_df = kmer_df.div(kmer_df.sum(axis=0), axis=1)
    kmer_df = kmer_df.reindex(sorted(kmers_visited))

    # make list of lists containing kmer frequencies for training
    contig_kmers = []
    for column in kmer_df:
        kmer_list = kmer_df[column].tolist()
        contig_kmers.append(kmer_list)

    # make a list containing labels for training
    prelabels = kmer_df.columns.tolist()
    labels = []
    for i in prelabels:
        if i.split("_")[0] == "euk":
            labels.append(1)
        else:
            labels.append(0)

    # save the list of canonical kmers for validation purposes
    with open('kmer_list.pkl', 'wb') as f:
        pickle.dump(kmers_visited, f)

    # split data for training and testing
    X_train, X_test, y_train, y_test = train_test_split(contig_kmers, labels, test_size=0.3, random_state=42)

    # fit the SVM
    clf = svm.SVC()
    clf.fit(X_train, y_train)

    # save the model
    model_outfile = 'kmer_SVM.sav'
    pickle.dump(clf, open(model_outfile, 'wb'))

    # run predictions on testing set and compute accuracy
    predictions = clf.predict(X_test)
    print(f"Accuracy in the testing set is {accuracy_score(y_test, predictions) * 100}%")
