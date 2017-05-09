#!/usr/bin/python

import os, logging, pysam, vcf


def open_new_vcf(new_vcf_name):
    vcf_out_reader = vcf.Reader(open(new_vcf_name, 'r'))
    vcf_out_writer = vcf.Writer(open(new_vcf_name, 'w'), vcf_out_reader)
    return vcf_out_writer


def create_tmp_vcf(filter_modus, filename, outdir, vcf_path, variant, variant_number):
    tmp_results_dir = os.path.abspath(outdir) + '/tmp/'
    vcf_file = pysam.VariantFile(vcf_path)
    new_header = vcf_file.header

    if filter_modus == 'germline':
        new_header = create_germline_vcf(new_header)
    if filter_modus == 'somatic':
        new_header = create_somatic_vcf(new_header)

    new_vcf_name = tmp_results_dir + os.path.basename(vcf_path)[:-3] + filename + '_' + str(variant_number) + '.vcf'
    new_vcf_file = pysam.VariantFile(new_vcf_name, 'w', header=new_header)
    new_vcf_file.close()
    logging.debug('Created tmp vcf file (%s)', new_vcf_name)

    vcf_out_reader = vcf.Reader(open(new_vcf_name, 'r'))
    vcf_out_writer = vcf.Writer(open(new_vcf_name, 'w'), vcf_out_reader)
    vcf_out_writer.write_record(variant)
    logging.debug('Tmp vcf file ready for writing')

    return new_vcf_name


def create_germline_vcf(new_header):
    logging.debug('Created germline file header')
    new_header.info.add('AmpFA', 1, "Float",
                        "Amplicon Fraction for Alternate allele. Combined for multiallelic variants")
    new_header.info.add('AmpFR', 1, "Float", "Amplicon Fraction for Reference allele.")
    new_header.info.add('AmpCR', 1, "Integer", "Amplicon Count for Reference allele.")
    new_header.info.add('AmpCA', 1, "Integer", "Amplicon Count for Alternate allele.")
    new_header.info.add('AmpC', 1, "Integer", "Total Amplicon Count")
    new_header.info.add('AmpF_OA', 1, "Integer",
                        "Amplicon Count offset compared to AD field, for alternate allele. Combined for multiallelic variants")
    new_header.info.add('AmpF_OR', 1, "Integer",
                        "Amplicon Count offset compared to AD field, for reference allele. Combined for multiallelic variants")
    return new_header


def create_somatic_vcf(new_header):
    logging.debug('Created somatic vcf file header')
    new_header.info.add('AmpFA_n', 1, "Float",
                        "Amplicon Fraction for Alternate allele of normal. Combined for multiallelic variants")
    new_header.info.add('AmpFR_n', 1, "Float", "Amplicon Fraction for Reference allele of normal.")
    new_header.info.add('AmpCR_n', 1, "Integer", "Amplicon Count for Reference allele of normal.")
    new_header.info.add('AmpCA_n', 1, "Integer", "Amplicon Count for Alternate allele of normal.")
    new_header.info.add('AmpC_n', 1, "Integer", "Total Amplicon Count")
    new_header.info.add('AmpFA_t', 1, "Float",
                        "Amplicon Fraction for Alternate allele of tumor. Combined for multiallelic variants")
    new_header.info.add('AmpFR_t', 1, "Float", "Amplicon Fraction for Reference allele of tumor.")
    new_header.info.add('AmpCR_t', 1, "Integer", "Amplicon Count for Reference allele of tumor.")
    new_header.info.add('AmpCA_t', 1, "Integer", "Amplicon Count for Alternate allele of tumor.")
    new_header.info.add('AmpC_t', 1, "Integer", "Total Amplicon Count of tumor")
    new_header.info.add('AmpF_OA', 1, "Integer",
                        "Amplicon Count offset compared to AD field, for alternate allele. Combined for multiallelic variants")
    new_header.info.add('AmpF_OR', 1, "Integer",
                        "Amplicon Count offset compared to AD field, for reference allele. Combined for multiallelic variants")
    return new_header


def create_vcf(filter_modus, filename, outdir, vcf_path, input_arguments):
    vcf_file = pysam.VariantFile(vcf_path)
    new_header = vcf_file.header

    if filter_modus == 'germline':
        new_header = create_germline_vcf(new_header)
    if filter_modus == 'somatic':
        new_header = create_somatic_vcf(new_header)

    new_vcf_name = outdir + os.path.basename(vcf_path)[:-3] + filename + '.vcf'
    new_vcf_file = pysam.VariantFile(new_vcf_name, 'w', header=new_header)
    new_vcf_file.close()
    input_arguments['filter_vcf'] = new_vcf_name
    logging.info('Created new %s vcf for filtered variants (%s)', filter_modus, new_vcf_name)

    # vcf_out_writer = open_new_vcf(new_vcf_name)
    logging.debug('New vcf ready for processing')

    return input_arguments


def open_unfiltered_vcf(vcf_path):
    unfiltered_vcf_file = vcf.Reader(open(vcf_path, 'r'))
    logging.debug('Opened unfiltered VCF')
    return unfiltered_vcf_file
