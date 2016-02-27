from skbio.alignment import StripedSmithWaterman

import kmers


SW_SCORE_THRESHOLD = 15


def get_brute_force_mapping(A_kmer_list, B_kmer_list):

    mapping_list = []

    for a_kmer in A_kmer_list:

        sw_score_generator = StripedSmithWaterman(a_kmer,
                                                  match_score=1,
                                                  mismatch_score=-1,
                                                  gap_open_penalty=1,
                                                  gap_extend_penalty=1,
                                                  mask_length=0)

        max_score = 0

        max_score_b_kmer = ""

        for b_kmer in B_kmer_list:

            sw_score = sw_score_generator(b_kmer)['optimal_alignment_score']

            if sw_score >= SW_SCORE_THRESHOLD and sw_score > max_score:

                max_score = sw_score

                max_score_b_kmer = b_kmer

        if max_score != 0:

            mapping_list.append((a_kmer, max_score, max_score_b_kmer))

    return mapping_list


def main():

    sequence_kmer_list = kmers.get_random_sequence_kmers(1000)

    design_kmer_list = kmers.get_design_kmers()

    mapping_list = get_brute_force_mapping(sequence_kmer_list, design_kmer_list)

    print(len(mapping_list))


if __name__ == "__main__":

    main()
