#!/usr/bin/python


def filter_variants_germline(config_parameters, number_after_amplicons_all, new_variant_info_field, variant, nr_amplicons_design,reads_all,pass_var_position_read_list):

    try:
        total_read_pos_ratio = (100 * (pass_var_position_read_list['passed']['bam_tumor']) / float(reads_all['total']['bam_tumor']))
    except ZeroDivisionError:
        total_read_pos_ratio = 100

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

    return variant