import pickle

import kmers

import Clustering

from BipartiteMatching import get_bipartite_matching


def output_seq_kmer_mapping_list(input_name, consensus_mapping_list, cluster_dict):

    sequence_kmer_mapping_list = []

    matching_score_sum = 0

    mapping_score_sum = 0

    for consensus_mapping in consensus_mapping_list:

        consensus_kmer = consensus_mapping[0]

        consensus_matching_score = consensus_mapping[1]
        matching_score_sum += consensus_matching_score

        design_kmer = consensus_mapping[2]

        sequence_kmer_list = cluster_dict[consensus_kmer]

        for sequence_kmer in sequence_kmer_list:

            mapping_score_sum += consensus_matching_score

            sequence_kmer_mapping_list.append((sequence_kmer, consensus_matching_score, design_kmer))

    output_list_file_name = input_name + "_mapping.p"
    pickle.dump(sequence_kmer_mapping_list, open(output_list_file_name, "wb"))

    score_dict = {"matching_score_sum": matching_score_sum, "mapping_score_sum": mapping_score_sum}
    output_score_file_name = input_name + "_score.p"
    pickle.dump(score_dict, open(output_score_file_name, "wb"))



def main():

    prefix_500 = "500"
    prefix_1000 = "1000"
    prefix_2000 = "2000"

    suffix_start_count = 1
    suffix_end_count = 10

    design_kmer_list = kmers.get_design_kmers()

    use_cluster_size_hard_stop = True

    for suffix_count in range(suffix_start_count, suffix_end_count + 1):

        input_file_name = prefix_500 + '_' + str(suffix_count)
        sequence_kmer_list = kmers.get_sequence_kmers(input_file_name)
        cluster_dict = Clustering.get_cluster_dict(sequence_kmer_list,
                                                   design_kmer_list,
                                                   use_cluster_size_hard_stop)

        consensus_kmer_list = cluster_dict.keys()
        consensus_mapping_list = get_bipartite_matching(consensus_kmer_list, design_kmer_list)
        output_seq_kmer_mapping_list(input_file_name, consensus_mapping_list, cluster_dict)

    for suffix_count in range(suffix_start_count, suffix_end_count + 1):
        input_file_name = prefix_1000 + '_' + str(suffix_count)
        sequence_kmer_list = kmers.get_sequence_kmers(input_file_name)
        cluster_dict = Clustering.get_cluster_dict(sequence_kmer_list,
                                                   design_kmer_list,
                                                   use_cluster_size_hard_stop)

        consensus_kmer_list = cluster_dict.keys()
        consensus_mapping_list = get_bipartite_matching(consensus_kmer_list, design_kmer_list)
        output_seq_kmer_mapping_list(input_file_name, consensus_mapping_list, cluster_dict)

    for suffix_count in range(suffix_start_count, suffix_end_count + 1):
        input_file_name = prefix_2000 + '_' + str(suffix_count)
        sequence_kmer_list = kmers.get_sequence_kmers(input_file_name)
        cluster_dict = Clustering.get_cluster_dict(sequence_kmer_list,
                                                   design_kmer_list,
                                                   use_cluster_size_hard_stop)

        consensus_kmer_list = cluster_dict.keys()
        consensus_mapping_list = get_bipartite_matching(consensus_kmer_list, design_kmer_list)
        output_seq_kmer_mapping_list(input_file_name, consensus_mapping_list, cluster_dict)


if __name__ == "__main__":

    main()