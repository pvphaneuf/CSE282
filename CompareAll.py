import kmers

from BruteForceMapping import get_brute_force_mapping

from BipartiteMatching import get_bipartite_matching

import Clustering

import ClusteringBestRandomConsensus

import ClusteringBestRandomClusterMerge

import pickle


class MappingMethod:

    BRUTE_FORCE = 0

    BIPARTITE = 1


class ClusteringFeatures:

    CLUSTER_SIZE_HARD_STOP = 0

    RANDOM_BEST_CONSENSUS = 1

    RANDOM_BEST_CLUSTER = 2


def get_mapping(sequence_kmer_list,
                design_kmer_list,
                mapping_paramenter,
                clustering_parameter_list):

    use_cluster_size_hard_stop = False
    if ClusteringFeatures.CLUSTER_SIZE_HARD_STOP in clustering_parameter_list:
        use_cluster_size_hard_stop = True

    use_random_best_consensus = False
    if ClusteringFeatures.RANDOM_BEST_CONSENSUS in clustering_parameter_list:
        use_random_best_consensus = True

    use_random_best_cluster_merge = False
    if ClusteringFeatures.RANDOM_BEST_CLUSTER in clustering_parameter_list:
        use_random_best_cluster_merge = True

    if not use_random_best_consensus and not use_random_best_cluster_merge:
        cluster_dict = Clustering.get_cluster_dict(sequence_kmer_list, design_kmer_list, use_cluster_size_hard_stop)

    if use_random_best_consensus and not use_random_best_cluster_merge:
        cluster_dict = ClusteringBestRandomConsensus.get_cluster_dict(sequence_kmer_list, design_kmer_list, use_cluster_size_hard_stop)

    if not use_random_best_consensus and use_random_best_cluster_merge:
        cluster_dict = ClusteringBestRandomClusterMerge.get_cluster_dict(sequence_kmer_list, design_kmer_list, use_cluster_size_hard_stop)

    consensus_kmer_list = cluster_dict.keys()

    if mapping_paramenter == MappingMethod.BRUTE_FORCE:

        consensus_mapping_list = get_brute_force_mapping(consensus_kmer_list, design_kmer_list)

    else:

        consensus_mapping_list = get_bipartite_matching(consensus_kmer_list, design_kmer_list)

    sequence_kmer_mapping_list = []

    for consensus_mapping in consensus_mapping_list:

        consensus_kmer = consensus_mapping[0]

        consensus_mapping_score = consensus_mapping[1]

        design_kmer = consensus_mapping[2]

        sequence_kmer_list = cluster_dict[consensus_kmer]

        for sequence_kmer in sequence_kmer_list:

            sequence_kmer_mapping_list.append((sequence_kmer, consensus_mapping_score, design_kmer))

    # pickle.dump(sequence_kmer_mapping_list, open("sequence_kmer_mapping_list.p", "wb"))

    return sequence_kmer_mapping_list


def main():

    sequence_kmer_list = kmers.get_random_sequence_kmers(1000)

    design_kmer_list = kmers.get_design_kmers()

    clustering_parameter_list = [ClusteringFeatures.CLUSTER_SIZE_HARD_STOP]

    mapping_paramenter = MappingMethod.BIPARTITE

    mapping_list = get_mapping(sequence_kmer_list,
                               design_kmer_list,
                               mapping_paramenter,
                               clustering_parameter_list)

    len(mapping_list)

    # brute_force_mapping_list = get_brute_force_mapping(sequence_kmer_list, design_kmer_list)
    #
    # bipartite_matching_list = get_bipartite_matching(sequence_kmer_list, design_kmer_list)

    # print(len(brute_force_mapping_list))
    #
    # print(len(bipartite_matching_list))


if __name__ == "__main__":

    main()
