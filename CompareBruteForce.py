import BruteForceMapping

import kmers

import pickle


def output_seq_kmer_mapping_list(input_name, mapping_list):

    mapping_score_sum = sum(mapping[1] for mapping in mapping_list)

    output_list_file_name = "brute_force_" + input_name + "_mapping.p"
    pickle.dump(mapping_list, open(output_list_file_name, "wb"))

    score_dict = {"mapping_score_sum": mapping_score_sum}
    output_score_file_name = "brute_force_" + input_name + "_score.p"
    pickle.dump(score_dict, open(output_score_file_name, "wb"))


def main():

    print("poop")

    prefix_500 = "500"
    prefix_1000 = "1000"
    prefix_2000 = "2000"

    suffix_start_count = 1
    suffix_end_count = 10

    design_kmer_list = kmers.get_design_kmers()

    for suffix_count in range(suffix_start_count, suffix_end_count + 1):

        input_file_name = prefix_500 + '_' + str(suffix_count)
        sequence_kmer_list = kmers.get_sequence_kmers(input_file_name)
        mapping_list = BruteForceMapping.get_brute_force_mapping(sequence_kmer_list, design_kmer_list)
        output_seq_kmer_mapping_list(input_file_name, mapping_list)

    for suffix_count in range(suffix_start_count, suffix_end_count + 1):
        input_file_name = prefix_1000 + '_' + str(suffix_count)
        sequence_kmer_list = kmers.get_sequence_kmers(input_file_name)
        mapping_list = BruteForceMapping.get_brute_force_mapping(sequence_kmer_list, design_kmer_list)
        output_seq_kmer_mapping_list(input_file_name, mapping_list)

    for suffix_count in range(suffix_start_count, suffix_end_count + 1):
        input_file_name = prefix_2000 + '_' + str(suffix_count)
        sequence_kmer_list = kmers.get_sequence_kmers(input_file_name)
        mapping_list = BruteForceMapping.get_brute_force_mapping(sequence_kmer_list, design_kmer_list)
        output_seq_kmer_mapping_list(input_file_name, mapping_list)


if __name__ == "__main__":

    main()