import kmers

from skbio.alignment import StripedSmithWaterman

import networkx


def get_bipartite_matching(S_kmer_list, V_kmer_list):

    S_node_list = []
    V_node_list = []
    weighted_edge_list = []
    node_kmer_dict = {}

    s_node_str_prefix = 's'
    v_node_str_prefix = 'v'

    alignment_score_threshold = 15

    for s_kmer_idx in range(len(S_kmer_list)):

        s_node_str = s_node_str_prefix + str(s_kmer_idx)

        if s_node_str not in S_node_list:
            S_node_list.append(s_node_str)

            s_kmer = S_kmer_list[s_kmer_idx]
            node_kmer_dict[s_node_str] = s_kmer


        s_kmer = S_kmer_list[s_kmer_idx]

        sw_score_generator = StripedSmithWaterman(s_kmer,
                                                  match_score = 1,
                                                  mismatch_score = -1,
                                                  gap_open_penalty = 1,
                                                  gap_extend_penalty = 1,
                                                  mask_length = 0)

        for v_kmer_idx in range(len(V_kmer_list)):

            v_node_str = v_node_str_prefix + str(v_kmer_idx)

            if v_node_str not in V_node_list:
                V_node_list.append(v_node_str)

                v_kmer = V_kmer_list[v_kmer_idx]
                node_kmer_dict[v_node_str] = v_kmer

            v_kmer = V_kmer_list[v_kmer_idx]

            sw_score = sw_score_generator(v_kmer)['optimal_alignment_score']

            if sw_score >= alignment_score_threshold:

                edge = (s_node_str, v_node_str, sw_score)

                weighted_edge_list.append(edge)

    bipartite_graph = networkx.Graph()
    bipartite_graph.add_nodes_from(S_node_list, bipartite=0)
    bipartite_graph.add_nodes_from(V_node_list, bipartite=1)
    bipartite_graph.add_weighted_edges_from(weighted_edge_list)

    matching_node_dict = networkx.max_weight_matching(bipartite_graph, maxcardinality=True)

    match_list = []
    for node_str_1, node_str_2 in matching_node_dict.items():

        if node_str_1[0] == 's':

            s_node_str = node_str_1
            v_node_str = node_str_2

            s_kmer = node_kmer_dict[s_node_str]
            v_kmer = node_kmer_dict[v_node_str]

            sw_score_generator = StripedSmithWaterman(s_kmer,
                                                  match_score = 1,
                                                  mismatch_score = -1,
                                                  gap_open_penalty = 1,
                                                  gap_extend_penalty = 1,
                                                  mask_length = 0)

            sw_score = sw_score_generator(v_kmer)['optimal_alignment_score']

            match_list.append((s_kmer, sw_score, v_kmer))

    return match_list


def main():

    sequence_kmer_list = kmers.get_random_sequence_kmers(1000)

    design_kmer_list = kmers.get_design_kmers()

    mapping_list = get_bipartite_matching(sequence_kmer_list, design_kmer_list)

    print(len(mapping_list))


if __name__ == "__main__":

    main()
