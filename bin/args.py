#!/usr/bin/python

import logging, argparse


def argument_checkup():
    # type: 
    """
    Argument parsing
    """

    parser = argparse.ArgumentParser(prog='pyAMPLI')
    subparsers = parser.add_subparsers(help='commands')
    germline_args = subparsers.add_parser('germline',
                                          help='Input arguments for germline amplicon filter')
    germline_args.add_argument('-b',
                               '--bam',
                               type=str,
                               help='BAM file. Index located in the same directory',
                               required=True)
    germline_args.add_argument('-v',
                               '--vcf',
                               type=str,
                               help='VCF file',
                               required=True)
    germline_args.add_argument('-d',
                               '--design',
                               type=str,
                               help='Probe/amplicon design file. Contains field, is this order: Contact your manufacturer for details',
                               required=True)
    germline_args.add_argument('-f',
                               '--filename',
                               type=str,
                               help='Output file name',
                               default='filtered')
    germline_args.add_argument('-od',
                               '--outdir',
                               type=str,
                               help='Output directory',
                               default='filter')

    somatic_args = subparsers.add_parser('somatic',
                                         help='Input arguments for somatic amplicon filter')
    somatic_args.add_argument('-bn',
                              '--bam_normal',
                              type=str,
                              help='BAM file of the normal sample. Index located in the same directory',
                              required=True)
    somatic_args.add_argument('-bt',
                              '--bam_tumor',
                              type=str,
                              help='BAM file of the tumor sample. Index located in the same directory',
                              required=True)
    somatic_args.add_argument('-v',
                              '--vcf',
                              type=str,
                              help='VCF file',
                              required=True)
    somatic_args.add_argument('-d',
                              '--design',
                              type=str,
                              help='Probe/amplicon design file. Contains field, is this order: Contact your manufacturer for details',
                              required=True)
    somatic_args.add_argument('-f',
                              '--filename',
                              type=str,
                              help='Output file name',
                              default='filtered')
    somatic_args.add_argument('-od',
                              '--outdir',
                              type=str,
                              help='Output directory',
                              default=None)

    input_arguments = parser.parse_args()
    logging.debug('Input_parameters: %s', input_arguments)

    return input_arguments
