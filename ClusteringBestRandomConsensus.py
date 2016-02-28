import kmers

import random


def findConsensusString(strings):
    nucleotides = "ACGT"
    profile = profileMatrix(strings,len(strings[0]))
    consensusString = ""
    for i in xrange(0,len(strings[0])):
        maxScore = 0
        maxIndex = 0

        for j in xrange(0,4):
            if profile[j][i] >= maxScore:
                maxScore = profile[j][i]
                maxIndex = j
        maxIndexArray = []
        for j in xrange(0,4):
            if profile[j][i] == maxScore:
                maxIndexArray.append(j)
        randomIndex = random.choice(maxIndexArray)
        consensusString += nucleotides[randomIndex]
    return consensusString


def findHierarchicalClusters(seqKmers,clusterSize):
    clusters = formInitialClusters(seqKmers)
    print len(clusters)
    print "initial clusters formed"
    while len(clusters) > clusterSize:
        closestClusters = findClosestClusters(clusters)
        clusterMembers = clusters[closestClusters[0]]
        clusterMembers.extend(clusters[closestClusters[1]])
        consensusString = findConsensusString(clusterMembers)
        clusters[consensusString] = clusterMembers
        del clusters[closestClusters[0]],clusters[closestClusters[1]]
        print len(clusters)
    return clusters


def get_cluster_dict(sequence_kmer_list, design_kmer_list):

    cluster_dict = findHierarchicalClusters(sequence_kmer_list, len(design_kmer_list))

    return cluster_dict


def main():

    sequence_kmer_list = kmers.get_random_sequence_kmers(1000)

    design_kmer_list = kmers.get_design_kmers()

    cluster_dict = get_cluster_dict(sequence_kmer_list, design_kmer_list)

    pickle.dump(cluster_dict, open("cluster_dict.p", "wb"))


if __name__ == "__main__":

    main()
