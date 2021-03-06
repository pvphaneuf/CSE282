import kmers

from collections import defaultdict

import sys

import copy

# import pickle


def formInitialClusters(seqKmers):
    initialClusters = defaultdict(list)
    for kmer in seqKmers:
        initialClusters[kmer].append(kmer)
    return initialClusters


def profileMatrix(motif,k):
    nucleotides = "ACGT"
    profile = [[0 for x in range(k)] for x in range(4)]
    for i in xrange(len(motif)):
        for j in xrange(k):
            profile[nucleotides.index(motif[i][j])][j] += 1
    for i in xrange(0,4):
        for j in xrange(0,k):
            #profile[i][j] += 1
            profile[i][j] = profile[i][j]/(1.0 * (len(motif)))
    return profile


def findConsensusString(strings):
    nucleotides = "ACGT"
    profile = profileMatrix(strings,len(strings[0]))
    consensusString = ""
    for i in xrange(0,len(strings[0])):
        maxScore = 0
        maxIndex = 0
        for j in xrange(0,4):
            if profile[j][i] > maxScore:
                maxScore = profile[j][i]
                maxIndex = j
        consensusString += nucleotides[maxIndex]
    return consensusString


def findHammingDistance(a,b):
    hammingDistance = 0
    for i in xrange(0,len(a)):
        if a[i] != b[i]:
            hammingDistance += 1
    return hammingDistance


def findClosestClusters(clusters, isClusterSizeHardStop):
    clusterList = clusters.keys()
    distanceMatrix = [[-1 for i in range(len(clusterList))] for j in range(len(clusterList))]
    minDistance = sys.maxint
    closestClusters = [0,0]
    for i in xrange(len(clusterList)):
        for j in xrange(len(clusterList)):
            if i != j and distanceMatrix[i][j] == -1:
                distanceMatrix[i][j] = findHammingDistance(clusterList[i],clusterList[j])
                distanceMatrix[j][i] = copy.copy(distanceMatrix[i][j])
                if distanceMatrix[i][j] < minDistance:
                    minDistance = distanceMatrix[i][j]
                    closestClusters[0] = clusterList[i]
                    closestClusters[1] = clusterList[j]
    # print "done finding closest clusters to merge"
    if not isClusterSizeHardStop and minDistance > 5:
        return None
    return closestClusters


def findHierarchicalClusters(seqKmers, clusterSize, isClusterSizeHardStop):

    clusters = formInitialClusters(seqKmers)

    while len(clusters) > clusterSize:

        closestClusters = findClosestClusters(clusters, isClusterSizeHardStop)

        if not isClusterSizeHardStop and closestClusters == None:
            break

        clusterMembers = clusters[closestClusters[0]]
        clusterMembers.extend(clusters[closestClusters[1]])

        consensusString = findConsensusString(clusterMembers)

        clusters[consensusString] = clusterMembers

        if closestClusters[0] != consensusString:
            del clusters[closestClusters[0]]

        if closestClusters[1] != consensusString:
            del clusters[closestClusters[1]]

    return clusters


def get_cluster_dict(sequence_kmer_list, design_kmer_list, isClusterSizeHardStop):

    cluster_dict = findHierarchicalClusters(sequence_kmer_list, len(design_kmer_list), isClusterSizeHardStop)

    return cluster_dict


def main():

    design_kmer_list = kmers.get_design_kmers()

    isClusterSizeHardStop = False

    input_file_name = "500_1"

    sequence_kmer_list = kmers.get_sequence_kmers(input_file_name)

    cluster_dict = get_cluster_dict(sequence_kmer_list,
                                                   design_kmer_list,
                                                   isClusterSizeHardStop)

    # pickle.dump(cluster_dict, open("cluster_dict.p", "wb"))


if __name__ == "__main__":

    main()
