import kmers

from BruteForceMatching import get_brute_force_mapping

from BipartiteMatching import get_bipartite_matching

from Clustering import get_cluster_dict


def get_brute_force_with_clustering_mapping(sequence_kmer_list, design_kmer_list):

    cluster_dict = get_cluster_dict(sequence_kmer_list, design_kmer_list)

    consensus_kmer_list = cluster_dict.keys()

    brute_force_mapping_list = get_brute_force_mapping(consensus_kmer_list, design_kmer_list)

    return brute_force_mapping_list


def main():

    sequence_kmer_list = kmers.get_random_sequence_kmers(1000)

    design_kmer_list = kmers.get_design_kmers()

    brute_force_mapping_list = get_brute_force_mapping(sequence_kmer_list, design_kmer_list)

    bipartite_matching_list = get_bipartite_matching(sequence_kmer_list, design_kmer_list)

    clustering_matching_list = get_brute_force_with_clustering_mapping(sequence_kmer_list, design_kmer_list)

    print(len(brute_force_mapping_list))

    print(len(bipartite_matching_list))

    print(len(clustering_matching_list))


if __name__ == "__main__":

    main()
