import random

import os


def get_random_sequence_kmers(count):

    input_file_path = os.path.split(__file__)[0]
    input_file_path += '/data/reverseKmers.txt'

    f = open(input_file_path, "r")
    seqKmers = []
    for line in f.readlines():
        kmer = line.strip()
        seqKmers.append(kmer)
    random.shuffle(seqKmers)
    seqKmers = seqKmers[0:1000]
    return seqKmers


def get_sequence_kmers(input_file):

    input_file_path = os.path.split(__file__)[0]
    input_file_path += '/data/' + input_file
    f = open(input_file_path, "r")

    sequence_kmer_list = []
    for line in f.readlines():
        kmer = line.strip()
        sequence_kmer_list.append(kmer)

    return sequence_kmer_list


def get_design_kmers():

    input_file_path = os.path.split(__file__)[0]
    input_file_path += '/data/design_kmers.txt'

    f = open(input_file_path, "r")
    designKmers = []
    for line in f.readlines():
        if line[0] != ">":
            designKmers.append(line.strip())
    return designKmers
