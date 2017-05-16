#!/usr/bin/python

import pysam
import re


def get_pad_alleles(alleles):
    padded_alleles = [[0 for x in range(2)] for x in range(len(alleles) - 1)]
    for i in range(1, len(alleles)):
        if len(alleles[0]) == len(alleles[i]):
            padded_alleles[i - 1][0] = alleles[0]
            padded_alleles[i - 1][1] = str(alleles[i])
        elif len(alleles[0]) > len(alleles[i]):
            padded_alleles[i - 1][0] = alleles[0]
            padded_alleles[i - 1][1] = str(alleles[i][0]) + ''.join('-' * (len(alleles[0]) - len(alleles[i])))
            if len(alleles[i]) > 1:
                padded_alleles[i - 1][1] += str(alleles[i][1:])
        elif len(alleles[0]) < len(alleles[i]):
            padded_alleles[i - 1][1] = alleles[i]
            padded_alleles[i - 1][0] = alleles[0][0] + ''.join('-' * (len(alleles[i]) - len(alleles[0])))
            if len(alleles[0]) > 1:
                padded_alleles[i - 1][0] += alleles[0][1:]
    return (padded_alleles)

def get_variant_type(filter_modus, variant, padded_alleles, var_type_count):
    for multi_allel_var in padded_alleles:
        a1 = str(multi_allel_var[0])
        a2 = str(multi_allel_var[1])
        var_type = 'total'
        var_type_count[var_type] = var_type_count[var_type] + 1
        if len(a1) == 1:
            if len(a2) == 1:
                var_type = 'snv'
                var_type_count[var_type] =var_type_count[var_type] + 1
            elif len(a2) > 1:
                var_type = 'ins'
                var_type_count[var_type] = var_type_count[var_type] + 1
            elif len(a2) == 0:
                print "Invalid deletion notation: " + a1 + "/" + a2
                var_type = 'del'
                var_type_count[var_type] = var_type_count[var_type] + 1
        elif len(a2) == 1:
            if len(a1) > 1:
                var_type = 'del'
                var_type_count[var_type] = var_type_count[var_type] + 1
            elif len(a1) == 0:
                print "Invalid insertion notation: " + a1 + "/" + a2
                var_type = 'ins'
                var_type_count[var_type] = var_type_count[var_type] + 1
        else:
            var_type = 'sub'
            var_type_count[var_type] = var_type_count[var_type] + 1

    if filter_modus == 'somatic':
        var_state = variant.INFO['SS']
        var_type_count[filter_modus][var_state] = var_type_count[filter_modus][var_state] + 1

    if len(padded_alleles) > 1:
        var_type = 'multi_allel'
        var_type_count[var_type] = var_type_count[var_type] + 1

    return var_type_count

def get_nr_amplicons(var_chrom, var_pos, design_by_start):
    acount = 0
    for ele_start in [x for x in design_by_start[var_chrom].keys() if int(x) <= int(var_pos)]:
        acount += len([x for x in design_by_start[var_chrom][ele_start].keys() if int(x) >= int(var_pos)])
    return acount

def check_var_position_read(read, position, var_read_pos_cutoff):
    try:
        fail_check = 1
        var_position_read_var = read.reference_end - position
        var_position_read_length = read.query_length
        var_position_read = var_position_read_length - var_position_read_var
        for index_var, ref_position in read.get_aligned_pairs(with_seq=False):
            if ref_position == position:
                # print "YES", index_var, ref_position, position, read.query_length, var_read_pos_cutoff
                if index_var >= (read. query_length - var_read_pos_cutoff):
                    fail_check = 0
                elif index_var <= (var_read_pos_cutoff):
                    fail_check = 0
                else:
                    fail_check = 1
                break
            else:
                pass
    except:
        fail_check = 1

    return fail_check

def check_read_for_variant(read ,variant ,ref_seqs ,seq_length ,padded_alleles, ref_fasta):
    # 0. edit distance == 0, reference read.
    ########################################
    # if read.get_tag('NM') == 0 and :
    #	print "nm == 0"
    #	print read.query_sequence.upper()
    #	return(0)
    # 	=> disabled to check for non-informative indels at end of read in stretch;

    # print "NM == 0, skip?"
    # 1. mismatches: process cigar (credit goes here : http://tinyurl.com/j3wa3us).
    ################################################################################
    # define variables
    read_idx = 0
    read_out_pos = 0
    ref_out_idx = 0
    read_out_idx = 0
    ref_out_pos = 0
    ref_idx = read.reference_start + 1  # read is zerobased, variant is one-based.
    ref_local_idx = 0
    read_length = read.query_length
    read_sequence = read.query_sequence.upper()
    has_variant = 0
    ref_out = ''
    read_out = ''
    # get reference_sequence
    try:
        region = '%s:%i-%i' % (read.reference_name, read.reference_start + 1, read.reference_end)
    except:
        # no/invalid cigar => no end_position => skip read.
        return (-1)

    if region in ref_seqs:
        ref_sequence = ref_seqs[region]
    else:
        # remove headerline and fuse multiple lines to a single line.
        ref_sequence = (''.join(pysam.faidx(ref_fasta, region))).split('\n', 1)[-1].replace('\n', '').upper()
        ref_seqs[region] = ref_sequence
    # go over cigar to build read/reference alignment.
    for operation, length in read.cigartuples:

        # hard clip and padding are not really in the read, skip
        if operation in (5, 6):
            continue
        # Insert or softclip : query_sequence that is not in reference_sequence
        elif operation in (4, 1):
            for i in range(0, length):
                if ref_idx == variant.POS:
                    read_out_pos = read_out_idx
                    ref_out_pos = ref_out_idx

                # move along the read.
                read_out += read_sequence[read_idx]
                read_out_idx += 1
                ref_out += '-'
                ref_out_idx += 1
                read_idx += 1

            continue
        ## deletion (or intron) : reference_sequence that is not in query_sequence
        elif operation in (2, 3):
            for i in range(0, length):
                if ref_idx == variant.POS:
                    read_out_pos = read_out_idx
                    ref_out_pos = ref_out_idx
                read_out += '-'
                read_out_idx += 1
                ref_out += ref_sequence[ref_local_idx]
                ref_out_idx += 1
                # move along reference
                ref_idx += 1
                ref_local_idx += 1
            continue
        ## alignment match. Can be "M" (general), "=" (match) or "X" (mismatch)
        elif operation in (0, 7, 8):
            for i in range(0, length):
                if ref_idx == variant.POS:
                    read_out_pos = read_out_idx
                    ref_out_pos = ref_out_idx
                ref_out += ref_sequence[ref_local_idx]
                ref_out_idx += 1
                read_out += read_sequence[read_idx]
                read_out_idx += 1

                ref_local_idx += 1
                ref_idx += 1
                read_idx += 1

    # 3. determine which allele the read holds.
    ###########################################
    matches = [0 for x in range(len(padded_alleles) + 1)]
    for i in range(0, len(padded_alleles)):
        alen = len(padded_alleles[i][0])
        # special indel  handling needed....
        if (str(padded_alleles[i][0]) + str(padded_alleles[i][1])).find('-') >= 0:
            if str(padded_alleles[i][0]).find('-') >= 0:
                pseq = str(padded_alleles[i][1])[str(padded_alleles[i][0]).find('-'):]
            elif str(padded_alleles[i][1]).find('-') >= 0:
                pseq = str(padded_alleles[i][0])[str(padded_alleles[i][1]).find('-'):]

            # we now have the sequence that is deleted/inserted
            # if stretch of pseq from pos to end of read => non informative
            # example:
            #  pseq : AC
            #  padded: GAC / G--
            #  ref : GCGATCGATGACACACACAC
            #  alt : GCGATCGATGACACACACAC
            #   => cannot be known if ref/alt due to no 3' fragment to align. is discarded by genotyper
            if re.search("^(" + pseq + ")+[" + pseq + "]{0," + str(len(pseq)) + "}$",
                         ref_out[(ref_out_pos + alen - len(pseq)):]):
                #	log.write("non informative stretch in "+ref_out+"\n")
                return (-1)

            # if stretch of pseq (e.g. homopolymer) in middle of read && length of stretch is equal in ref+read ==> reference read.
            # example:
            # pseq:  A
            # padded:  [['C-', 'CA']]
            # ref :  CCAAAAAAAAAAAATCCCC
            # read:  CCAAAAAAAAAAAATCCCC
            #  => would match alternate using logic below.
            target = ref_out[(ref_out_pos + alen - len(pseq)):]

            try:
                m = re.search("^((?:" + pseq + "){2,}[" + pseq + "]{0," + str(len(pseq)) + "}?)", target)
                stretch_length = len(m.group(0))
                if ref_out[ref_out_pos:(ref_out_pos + stretch_length)] == read_out[read_out_pos:(read_out_pos + stretch_length)]:
                    matches[0] += 1
                    continue
            except:
                pass

            # insertion without a stretch : reference if read == ref for span of alleles
            if ref_out[ref_out_pos:(ref_out_pos + max(5, alen))] == read_out[read_out_pos:(read_out_pos + max(5, alen))]:
                matches[0] += 1
                continue

        # reference for snv or del
        if read_out[read_out_pos:(read_out_pos + alen)] == padded_alleles[i][0]:
            matches[0] += 1
        # alt for snv or del or ins
        elif read_out[read_out_pos:(read_out_pos + alen)] == padded_alleles[i][1]:
            matches[i + 1] += 1

    if max(matches[1:]) > 0:
        has_variant = 1
    elif max(matches) == 0:
        has_variant = -1
    return (has_variant)

