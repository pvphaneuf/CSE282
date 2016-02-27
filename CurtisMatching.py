def spacerSW_trunc(sequencing_files, seqFile, barcodes, nprocs, spacer2seq):
    def runSW_trunc(record_set, sequencing_files, seqFile, spacer2seq, barcodes, procNum, out_queue, verbose=False):
        # import swalign
        # from ssw_wrap import Aligner

        # sw_handle = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(1, -1))
        # Find all record ids in the correct order
        # import pydevd
        # # pydevd.settrace('137.110.95.186', port=8890, stdoutToServer=True, stderrToServer=True)
        # pydevd.settrace('137.110.17.249', port=8890, stdoutToServer=True, stderrToServer=True)
        if sequencing_files[seqFile][-5:] != 'fastq':
            record_dict = SeqIO.index(sequencing_files[seqFile][:-8] + 'fastq', 'fastq')
        else:
            record_dict = SeqIO.index(sequencing_files[seqFile], 'fastq')


        maxes = []
        for i, r in enumerate(record_set):

            percentDone = float(i)/len(record_set)*100.

            if percentDone % 1 < .5 and float(i - 1)/len(record_set)*100. % 1 > .5:
                newLine = "SW: from process %d, %d percent done" % (procNum, percentDone)
                sys.stdout.write("\r\x1b[K"+newLine.__str__())
                sys.stdout.flush()

            record = record_dict[r]

            maxK = ''
            maxAln = 0
            altSpacers = []
            aligned_query = ''

            if seqFile == 1:
                # If it is the second of a mate-paired file, then align to the reverse complement of the sequence
                sw_barcodes = StripedSmithWaterman(str(record.seq.reverse_complement()), match_score=1, mismatch_score=-1,
                                     gap_open_penalty=1, gap_extend_penalty=1, mask_length=0)


                barcode2_aln = sw_barcodes(barcodes[2])
                if barcode2_aln.optimal_alignment_score< 90:
                    barcode1_aln = sw_barcodes(barcodes[0])
                    if barcode1_aln.optimal_alignment_score<24:
                        maxes += [(r, [], 0, [], [])]
                        continue
                    else:
                        spacer_begin = barcode1_aln.query_end+1
                        spacer_end = barcode1_aln.query_end+21
                else:
                    spacer_begin = max(0, barcode2_aln.query_begin -20)
                    spacer_end = barcode2_aln.query_begin
                if spacer_begin> spacer_end:
                    maxes += [(r, [], 0, [], [])]
                    continue
                sw = StripedSmithWaterman(str(record.seq.reverse_complement()[spacer_begin : spacer_end]), match_score=1, mismatch_score=-1,
                                     gap_open_penalty=1, gap_extend_penalty=1)
                # sw = StripedSmithWaterman(str(record.seq.reverse_complement()[spacer_begin : spacer_end]), match_score=1, mismatch_score=-1,
                #                      gap_open_penalty=1, gap_extend_penalty=1, mask_length=0)




                # sw = Aligner(record.seq.reverse_complement(), match=1, mismatch=1, gap_open=1,
                #              gap_extend=1, report_cigar=False, report_secondary=False)

            else:
                sw_barcodes = StripedSmithWaterman(str(record.seq), match_score=1, mismatch_score=-1,
                                     gap_open_penalty=1, gap_extend_penalty=1, mask_length=0)


                barcode2_aln = sw_barcodes(barcodes[2][:82])
                if barcode2_aln.optimal_alignment_score< 60:
                    barcode1_aln = sw_barcodes(barcodes[0])
                    if barcode1_aln.optimal_alignment_score<24:
                        maxes += [(r, [], 0, [], [])]
                        continue
                    else:
                        spacer_begin = barcode1_aln.query_end+1
                        spacer_end = barcode1_aln.query_end+21
                else:
                    spacer_begin = max(0, barcode2_aln.query_begin -20)
                    spacer_end = barcode2_aln.query_begin
                if spacer_begin> spacer_end:
                    maxes += [(r, [], 0, [], [])]
                    continue



                #todo implement QC based on barcode alignment score

                sw = StripedSmithWaterman(str(record.seq[spacer_begin : spacer_end]), match_score=1, mismatch_score=-1,
                                     gap_open_penalty=1, gap_extend_penalty=1)
                # sw = StripedSmithWaterman(str(record.seq[spacer_begin : spacer_end]), match_score=1, mismatch_score=-1,
                #                      gap_open_penalty=1, gap_extend_penalty=1, mask_length=0)
                # sw = Aligner(record.seq, match=1, mismatch=1, gap_open=1,
                #              gap_extend=1, report_cigar=False, report_secondary=False)

            for k in spacer2seq:

                # if seqFile == 0:
                #     aln = sw(spacer2seq[k] + barcodes[1][:3]) # score_filter and distance_filter can be set
                # if seqFile == 1:
                #     aln = sw(spacer2seq[k] + barcodes[2][:3]) # score_filter and distance_filter can be set

                # ipdb.set_trace()
                # todo for reverse reads, and check if alignment extends into the 2nd barcode
                aln = sw(spacer2seq[k])
                if aln.optimal_alignment_score >= maxAln:
                    if aln.optimal_alignment_score == maxAln:
                        altSpacers += [k]
                    else:
                        maxAln = aln.optimal_alignment_score
                        maxK = k
                        altSpacers = []
                        aligned_query = aln.aligned_query_sequence

            # Employ tie-breaker, add the 20 preceding bp and find highest score
            if altSpacers != [] and maxAln > 0:
                altSpacers = [maxK] + altSpacers

                if seqFile == 0:
                    altScores = [sw(barcodes[0][-20:] + spacer2seq[g]).optimal_alignment_score for g in altSpacers]
                else:
                    altScores = [sw(spacer2seq[g] + barcodes[-1][:20]).optimal_alignment_score for g in altSpacers]
                maxAltScore = max(altScores)
                altSpacers = [altSpacers[j] for j, g in enumerate(altScores) if g == maxAltScore]

                if len(altSpacers) == 1:
                    maxK = altSpacers[0]
                    altSpacers = []
                    aligned_query = sw(spacer2seq[maxK]).aligned_query_sequence

            if verbose:

                print altSpacers
                print maxK
                print i

                nameSet = set()
                for f in altSpacers:
                    nameSet = nameSet.union([f.split('_')[0]])

                if maxAln > 18 and len(altSpacers) > 1:
                    aln = sw(spacer2seq[maxK])
                    print aln.aligned_target_sequence
                    print aln.aligned_query_sequence
                    # print sw_handle.align(record.seq, spacer2seq[maxK]).dump()
            maxes += [(r, maxK, maxAln, altSpacers, aligned_query)]

        out_queue.put(maxes)

    # spacer2seq = spacer2sequence(dat)

    # Find all record ids in the correct order
    if sequencing_files[seqFile][-5:] != 'fastq':
        handle = open(sequencing_files[seqFile][:-8] + 'fastq', "rU")
    else:
        handle = open(sequencing_files[seqFile], "rU")
    records = [record.id for record in SeqIO.parse(handle, "fastq")]
    handle.close()

    # Set up the parallelization
    out_q = RetryQueue()
    chunksize = int(len(records) / float(nprocs))
    procs = []

    allMaxes = []
    start = time.time()

    # ipdb.set_trace()
    # runSW(records[chunksize * 3: chunksize * (3 + 1)], sequencing_files, seqFile, spacer2seq, barcodes2, 3, out_q, verbose=True)
    # ipdb.set_trace()


    # runSW(records[:10], sequencing_files, seqFile, spacer2seq, barcodes, 0, out_q)

    try:
        for i in range(nprocs):
            if i != nprocs - 1:
                p = multiprocessing.Process(
                    target=runSW_trunc,
                    args=(records[chunksize * i: chunksize * (i + 1)],
                          sequencing_files,
                          seqFile,
                          spacer2seq,
                          barcodes,
                          i,
                          out_q))
            else:
                p = multiprocessing.Process(
                    target=runSW_trunc,
                    args=(records[chunksize * i: len(records)],
                          sequencing_files,
                          seqFile,
                          spacer2seq,
                          barcodes,
                          i,
                          out_q))
            procs.append(p)


        for p in procs:
            p.start()

        for i in range(nprocs):
            allMaxes += out_q.get()

        # Wait for all worker processes to finish
        for p in procs:
            p.join()

        end = time.time()

        print "\nSW done for %d reads in %d seconds" % (len(records), end - start)

    except KeyboardInterrupt:

        print "Keyboard interrupt sent!"
        return None


    # allScores = [r[2] for r in allMaxes]
    # plt.figure()
    # bins = np.linspace(min(allScores), max(allScores), int(max(allScores) - min(allScores) + 1))
    # plt.hist(allScores, bins, alpha=1.0)
    # plt.show()

    return allMaxes