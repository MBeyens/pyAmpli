#!/usr/bin/python

import os, sys, logging, pysam
from pyampli import read_methods, variant_methods


def bam_index_file(bam_file):
    logging.info('Indexing BAM file (%s)', bam_file)
    pysam.index(bam_file)


def time_bam_bai(bam_file):
    try:
        print "BAMMM ::", os.stat(bam_file + '.bai')
        if os.stat(bam_file).st_mtime >= os.stat(bam_file + '.bai').st_mtime:
            logging.error('BAI file timestamp is older than BAM. Please re-index your BAM file (%s)', bam_file)
            sys.exit(0)
    except OSError:
        logging.error('No index of BAM file found. pyAmpli will re-index your BAM file now (%s)', bam_file)
        bam_index_file(bam_file)


def bam_loader(input_arguments):
    seq_lengths = {'bam_normal': dict(), 'bam_tumor': dict()}
    bam_file_match = ['bam_normal', 'bam_tumor']
    logging.info('Checking BAM/BAI files and calculating sequence lengths')
    for bam_file_type in input_arguments.keys():
        if bam_file_type in bam_file_match:
            bam_file_path = input_arguments[bam_file_type]
            time_bam_bai(bam_file_path)
            bam_file = pysam.AlignmentFile(bam_file_path, "rb")
            for ele in bam_file.header['SQ']:
                seq_lengths[bam_file_type][ele['SN']] = ele['LN']
            bam_file.close()
    return seq_lengths


def bam_read_var_checker(input_arguments, variant, mapping_quality, base_quality, design_by_start,
                         design_by_end, seq_lengths, padded_alleles, ref_fasta, has_variant,
                         reads_all, amplicons_all, pass_var_position_read_list, var_read_pos_cutoff):
    bam_file_match = ['bam_normal', 'bam_tumor']
    logging.debug('Checking reads: 1) reference and alternative nucleotides; 2) amplicon range')
    for bam_file_type in input_arguments.keys():
        if bam_file_type in bam_file_match:
            bam_file_path = input_arguments[bam_file_type]
            bam_file = pysam.AlignmentFile(bam_file_path, "rb")
            ref_seqs = dict()
            for read in bam_file.fetch(variant.CHROM, variant.POS - 1, variant.POS):
                reads_all['total'][bam_file_type] = reads_all['total'][bam_file_type] + 1
                pass_quality_read = read_methods.check_quality_read(read, variant.POS, mapping_quality, base_quality)
                if pass_quality_read == 1:
                    reads_all['total_quality_reads'][bam_file_type] = reads_all['total_quality_reads'][bam_file_type] + 1
                    amplicon_name = read_methods.get_exact_amplicon(read, design_by_start, design_by_end)
                    has_variant[bam_file_type] = variant_methods.check_read_for_variant(read, variant, ref_seqs, seq_lengths[bam_file_type][variant.CHROM], padded_alleles, ref_fasta)

                    if has_variant[bam_file_type] == 0:
                        if not amplicon_name in amplicons_all['ref_amps'][bam_file_type]:
                            amplicons_all['ref_amps'][bam_file_type][amplicon_name] = 1
                        else:
                            amplicons_all['ref_amps'][bam_file_type][amplicon_name] += 1
                        reads_all['refs'][bam_file_type] = reads_all['refs'][bam_file_type] + 1

                    elif has_variant[bam_file_type] == 1:
                        if not amplicon_name in amplicons_all['alt_amps'][bam_file_type]:
                            amplicons_all['alt_amps'][bam_file_type][amplicon_name] = 1
                        else:
                            amplicons_all['alt_amps'][bam_file_type][amplicon_name] += 1
                        reads_all['alts'][bam_file_type] = reads_all['alts'][bam_file_type] + 1

                    elif has_variant[bam_file_type] == -1:
                        reads_all['other_allel'][bam_file_type] = reads_all['other_allel'][bam_file_type] + 1

                    if has_variant[bam_file_type] == 1 or has_variant[bam_file_type] == -1:
                        reads_all['total_read_var'][bam_file_type] = reads_all['total_read_var'][bam_file_type] + 1
                        pass_var_position_read = variant_methods.check_var_position_read(read, variant.POS, var_read_pos_cutoff)
                        if pass_var_position_read == 1:
                            pass_var_position_read_list['passed'][bam_file_type] = pass_var_position_read_list['passed'][bam_file_type] + 1
                        else:
                            pass_var_position_read_list['failed'][bam_file_type] = pass_var_position_read_list['failed'][bam_file_type] + 1
                else:
                    reads_all['total_failed_reads'][bam_file_type] = reads_all['total_failed_reads'][bam_file_type] + 1


    return reads_all, amplicons_all, pass_var_position_read_list
