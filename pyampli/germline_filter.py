#!/usr/bin/python


def filter_variants_germline(config_parameters, number_after_amplicons_all, new_variant_info_field, variant, nr_amplicons_design,reads_all,pass_var_position_read_list):
    try:
        total_read_pos_ratio = (100 * (pass_var_position_read_list['passed']['bam_normal']) / float(reads_all['total']['bam_normal']))
    except ZeroDivisionError:
        total_read_pos_ratio = 100
    print total_read_pos_ratio

    if variant.samples[0]['DP'] < int(config_parameters['germline_settings']['min_depth_normal']):
        variant.add_filter('DepthFail')

    elif total_read_pos_ratio <= 10:
        variant.add_filter('PositionFail')

    elif nr_amplicons_design < 2:
        variant.add_filter('OneAmpPass')

    elif number_after_amplicons_all['alt_amps']['bam_normal'] < int(config_parameters['general_settings']['min_amp']) and new_variant_info_field['ampF_A']['bam_normal'] < float(config_parameters['general_settings']['min_frac']):
        variant.add_filter('LowAmpFail')

    else:
        if number_after_amplicons_all['total_amps']['bam_normal'] <= 2:
            if number_after_amplicons_all['ref_amps']['bam_normal'] == 2 and number_after_amplicons_all['alt_amps']['bam_normal'] == 1 and number_after_amplicons_all['total_amps']['bam_normal'] >= 1:
                variant.add_filter('MatchAmpPass')

            elif number_after_amplicons_all['ref_amps']['bam_normal'] == 1 and number_after_amplicons_all['alt_amps']['bam_normal'] == 2 and number_after_amplicons_all['total_amps']['bam_normal'] >= 1:
                variant.add_filter('MatchAmpPass')

            elif number_after_amplicons_all['ref_amps']['bam_normal'] == 1 and number_after_amplicons_all['alt_amps']['bam_normal'] == 1 and number_after_amplicons_all['total_amps']['bam_normal'] == 1:
                variant.add_filter('MatchAmpPass')

            elif number_after_amplicons_all['ref_amps']['bam_normal'] == 2 and number_after_amplicons_all['alt_amps']['bam_normal'] == 2 and number_after_amplicons_all['total_amps']['bam_normal'] == 2:
                variant.add_filter('MatchAmpPass')

            else:
                variant.add_filter('LowAmpFail')

        elif number_after_amplicons_all['total_amps']['bam_normal'] > 2:
            variant.add_filter('AmpPass')

        else:
            variant.add_filter('LowAmpFail')

    return variant

def add_variant_info_fields(number_after_amplicons_all, new_variant_info_field, variant):
    variant.add_info('AmpFA', new_variant_info_field['ampF_A']['bam_normal'])
    variant.add_info('AmpFR', new_variant_info_field['ampF_R']['bam_normal'])
    variant.add_info('AmpCA', number_after_amplicons_all['alt_amps']['bam_normal'])
    variant.add_info('AmpCR', number_after_amplicons_all['ref_amps']['bam_normal'])
    variant.add_info('AmpC', number_after_amplicons_all['total_amps']['bam_normal'])

    return variant