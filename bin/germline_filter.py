#!/usr/bin/python

import logging

def filter_variants_germline(config_parameters, number_after_amplicons_all, new_variant_info_field, variant, nr_amplicons_design):
    if number_after_amplicons_all['alt_amps']['bam_normal'] < int(config_parameters['general_settings']['min_amp']) and new_variant_info_field['ampF_A']['bam_normal'] < float(config_parameters['general_settings']['min_frac']):
        variant.add_filter('LowAmp')
    else:
        variant.add_filter('PASS')

    return variant

def add_variant_info_fields(number_after_amplicons_all, new_variant_info_field, variant):
    variant.add_info('AmpFA', new_variant_info_field['ampF_A']['bam_normal'])
    variant.add_info('AmpFR', new_variant_info_field['ampF_R']['bam_normal'])
    variant.add_info('AmpCA', number_after_amplicons_all['alt_amps']['bam_normal'])
    variant.add_info('AmpCR', number_after_amplicons_all['ref_amps']['bam_normal'])
    variant.add_info('AmpC', number_after_amplicons_all['total_amps']['bam_normal'])

    if 'AD' in variant.FORMAT and len(variant.alleles) == 2:
        # some stats.
        offset[var_type]['ref'].append(abs(refs - variant.samples[0]['AD'][0]))
        dp[var_type]['ref'].append(variant.samples[0]['AD'][0])
        offset[var_type]['alt'].append(abs(variant.samples[0]['AD'][1]))
        dp[var_type]['alt'].append(variant.samples[0]['AD'][1])
        # add to variant record.
        AmpF_O_R = refs - variant.samples[0]['AD'][0]
        AmpF_O_A = alts - variant.samples[0]['AD'][1]
        variant.add_info('AmpF_OR', AmpF_O_R)
        variant.add_info('AmpF_OA', AmpF_O_A)

    return variant