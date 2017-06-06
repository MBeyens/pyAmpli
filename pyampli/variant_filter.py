#!/usr/bin/python

import os, logging, multiprocessing, subprocess
from pyampli import variant_methods, process_bam, germline_filter, somatic_filter, make_vcf, progress


def multi_processing_variants(filter_modus, input_arguments, config_parameters, unfiltered_vcf_file, design_by_start,
                              design_by_end, seq_lengths):
    logging.info('Multiprocessing variant filtering on %s cores', config_parameters['multi_cores'])

    cmd1 = 'egrep -v "^#" %s | wc -l' % (input_arguments['vcf'])
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True)
    amount_total_variants, error1 = p1.communicate()

    process_count = 0
    variant_number = 0
    #progress.init_progressbar()
    for variant in unfiltered_vcf_file:
        worker_arguments = variant, filter_modus, input_arguments, config_parameters, design_by_start, design_by_end, seq_lengths, variant_number
        pool = multiprocessing.Process(target=extract_variant, args=worker_arguments)
        pool.start()
        process_count += 1
        variant_number += 1
        if process_count == int(config_parameters['multi_cores']):
            pool.join()
            process_count = 0
        progress.update_progress(variant_number, amount_total_variants)
    # Wait for all processes to be closed
    pool.join()

    return variant_number


def extract_variant(variant, filter_modus, input_arguments, config_parameters, design_by_start, design_by_end,
                    seq_lengths, variant_number):
    var_type_count = {'snv': 0, 'ins': 0, 'del': 0, 'sub': 0, 'total': 0, 'multi_allel': 0,
                      'somatic': {'0': 0, '1': 0, '2': 0, '3': 0, '5': 0},
                      'ref': {'A': 0, 'G': 0, 'T': 0, 'C': 0},
                      'alt': {'A': 0, 'G': 0, 'T': 0, 'C': 0}}
    amplicons_all = {'ref_amps': {'bam_normal': dict(), 'bam_tumor': dict()},
                     'alt_amps': {'bam_normal': dict(), 'bam_tumor': dict()}}
    reads_all = {'total': {'bam_normal': 0, 'bam_tumor': 0},
                 'total_read_var': {'bam_normal': 0, 'bam_tumor': 0},
                 'total_quality_reads': {'bam_normal': 0, 'bam_tumor': 0},
                 'total_failed_reads': {'bam_normal': 0, 'bam_tumor': 0},
                 'refs': {'bam_normal': 0, 'bam_tumor': 0},
                 'alts': {'bam_normal': 0, 'bam_tumor': 0},
                 'other_allel': {'bam_normal': 0, 'bam_tumor': 0}}

    pass_var_position_read_list = {'passed': {'bam_normal': 0, 'bam_tumor': 0},
                                   'failed': {'bam_normal': 0, 'bam_tumor': 0}}

    has_variant = {'bam_normal': 0, 'bam_tumor': 0}

    ref_seqs = dict()
    var_allel_list = list()

    for e in variant.alleles:
        var_allel_list.append(str(e))
    padded_alleles = variant_methods.get_pad_alleles(var_allel_list)
    gt = variant.samples[0]['GT'].replace('|', '/').split('/')

    var_type_count = variant_methods.get_variant_type(filter_modus, variant, padded_alleles, var_type_count)

    nr_amplicons_design = variant_methods.get_nr_amplicons(variant.CHROM, variant.POS, design_by_start)

    after_reads_all, after_amplicons_all, after_pass_var_position_read_list = process_bam.bam_read_var_checker(
        input_arguments,
        variant,
        config_parameters['general_settings']['mq'],
        config_parameters['general_settings']['bq'],
        design_by_start,
        design_by_end,
        seq_lengths,
        padded_alleles,
        config_parameters['reference'],
        has_variant,
        reads_all,
        amplicons_all,
        pass_var_position_read_list,
        config_parameters['general_settings']['min_read_pos'])

    number_after_amplicons_all = {'ref_amps': {'bam_normal': 0, 'bam_tumor': 0},
                                  'alt_amps': {'bam_normal': 0, 'bam_tumor': 0},
                                  'total_amps': {'bam_normal': 0, 'bam_tumor': 0}}
    new_variant_info_field = {'ampF_R': {'bam_normal': 0, 'bam_tumor': 0},
                              'ampF_A': {'bam_normal': 0, 'bam_tumor': 0}}

    bam_file_match = ['bam_normal', 'bam_tumor']
    for bam_file_type in input_arguments.keys():
        if bam_file_type in bam_file_match:
            nr_ref_amp = len(after_amplicons_all['ref_amps'][bam_file_type].keys())
            if 'NA' in after_amplicons_all['ref_amps'][bam_file_type].keys():
                nr_ref_amp -= 1
            number_after_amplicons_all['ref_amps'][bam_file_type] = nr_ref_amp

            nr_alt_amp = len(after_amplicons_all['alt_amps'][bam_file_type].keys())
            if 'NA' in after_amplicons_all['alt_amps'][bam_file_type].keys():
                nr_alt_amp -= 1
            number_after_amplicons_all['alt_amps'][bam_file_type] = nr_alt_amp

            l = set(after_amplicons_all['ref_amps'][bam_file_type].keys() + after_amplicons_all['alt_amps'][
                bam_file_type].keys())
            nr_amps_seen = len(l)
            if 'NA' in l:
                nr_amps_seen -= 1

            if nr_amps_seen > 0:
                new_variant_info_field['ampF_R'][bam_file_type] = nr_ref_amp / float(nr_amps_seen)
                new_variant_info_field['ampF_A'][bam_file_type] = nr_alt_amp / float(nr_amps_seen)
            else:
                new_variant_info_field['ampF_R'][bam_file_type] = new_variant_info_field['ampF_A'][bam_file_type] = 0

	    number_after_amplicons_all['total_amps'][bam_file_type] = nr_amps_seen

    if filter_modus == 'germline':
        variant = germline_filter.filter_variants_germline(config_parameters,
                                                           number_after_amplicons_all,
                                                           new_variant_info_field,
                                                           variant,
                                                           nr_amplicons_design,
                                                           reads_all,
                                                           pass_var_position_read_list)
        variant = germline_filter.add_variant_info_fields(number_after_amplicons_all,
                                                          new_variant_info_field,
                                                          variant)
    if filter_modus == 'somatic':
        variant = somatic_filter.filter_variants_somatic(config_parameters,
                                                         number_after_amplicons_all,
                                                         new_variant_info_field,
                                                         variant, nr_amplicons_design,
                                                         reads_all,
                                                         pass_var_position_read_list)
        print 'this is a variant :: ', variant
        variant = somatic_filter.add_variant_info_fields(number_after_amplicons_all,
                                                         new_variant_info_field, variant)

    new_vcf_name = make_vcf.create_tmp_vcf(filter_modus,
                                           input_arguments['filename'],
                                           input_arguments['outdir'],
                                           input_arguments['vcf'],
                                           variant,
                                           variant_number)


def merge_tmp_vcfs(input_arguments, variant_number):
    tmp_results_dir_vcfs = os.path.abspath(input_arguments['outdir']) + '/tmp/*.vcf'
    tmp_results_dir = os.path.abspath(input_arguments['outdir']) + '/tmp/'
    cmd_tmp_vcfs = 'grep -v '"'#'"' --no-filename ' + tmp_results_dir_vcfs + ' | sort -k1,1V -k2,2n >> ' + \
                   input_arguments['filter_vcf'] + ' && rm -Rf ' + tmp_results_dir
    os.system(cmd_tmp_vcfs)
    logging.info('Merged %s temporary vcf files to one vcf (%s)', variant_number, input_arguments['filter_vcf'])
