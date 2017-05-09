#!/usr/bin/python

import logging


def load_design_file(bed_file):
    design_by_start = dict()
    design_by_end = dict()

    logging.info('Reading in design file')

    fh = open(bed_file, 'r')
    for line in fh:
        if line[0:7] == 'browser' or line[0:5] == 'track':
            continue
        line = line.rstrip()
        chrom, start, end, name, score, strand = line.split("\t")
        # new chromosome
        if not chrom in design_by_start:
            design_by_start[chrom] = dict()
            design_by_end[chrom] = dict()
        # new start
        if not start in design_by_start[chrom]:
            design_by_start[chrom][start] = dict()
        # new end
        if not end in design_by_end[chrom]:
            design_by_end[chrom][end] = dict()
        # store amplicon
        design_by_start[chrom][start][end] = name
        # design_counts[name] = 0 # counts
        design_by_end[chrom][end][start] = name
    # print design_by_end
    fh.close

    logging.debug('Created design dictionaries')

    return design_by_start, design_by_end


def check_design_conflicts(fuzziness, design_by_start):
    logging.info('Checking design conflicts (fuzziness : %s)', fuzziness)

    for chrom in sorted(design_by_start):
        print chrom
        last_start = 0
        # print type(last_start)
        # print type(fuzziness)

        for start in sorted(design_by_start[chrom]):
            if last_start + fuzziness >= int(start) - fuzziness:
                print "Collision by start for", str(chrom), "start", str(last_start), "and", str(start), "."
            last_start = int(start)
            last_end = 0
            for end in sorted(design_by_start[chrom][start]):
                if last_end >= int(end) - fuzziness:
                    print "Collision by end for identical start", str(chrom) + ":" + str(start), " end ", str(
                        last_end), "and", str(end), "."
                last_end = int(end)

    logging.debug('Finished checking design conflicts')
