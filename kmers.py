import random


def get_random_sequence_kmers(count):
    f = open("reverseKmers.txt", "r")
    seqKmers = []
    for line in f.readlines():
        kmer = line.strip()
        seqKmers.append(kmer)
    random.shuffle(seqKmers)
    seqKmers = seqKmers[0:1000]
    return seqKmers


def get_design_kmers():
    f = open("design_kmers.txt", "r")
    designKmers = []
    for line in f.readlines():
        if line[0] != ">":
            designKmers.append(line.strip())
    return designKmers
