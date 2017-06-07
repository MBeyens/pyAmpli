#!/usr/bin/python

import logging, os, subprocess


def check_extension(file_type, infile):
    decomp_file_name = False
    if infile.endswith(str(file_type)):
        logging.debug('File extension as expected %s %s', file_type, infile)
    elif infile.endswith('.gz'):
        decomp_file_name = decompress_file(infile)
        logging.info('Decompressed file %s', decomp_file_name)
    else:
        logging.error('File extension (%s) is not correct or file compression is not supported (%s)', file_type, infile)
    return decomp_file_name


def check_file(file_type, infile):
    if os.path.isfile(infile):
        logging.info('Found %s file %s', file_type, infile)
        return 0
    else:
        if os.path.isdir(os.path.dirname(infile)):
            #logging.info('Cannot find %s file %s', file_type, infile)
            file_not_found = 'Cannot find %s file on location %s -> although directory exists!' % (file_type, infile)
        else:
            #logging.info('Cannot find %s file %s because directory does not exist!', file_type, infile)
            file_not_found = 'Cannot find %s file %s -> directory does not exist!' % (file_type, infile)
        return file_not_found

def decompress_file(infile):
    # Only works for .gz files atm
    outfile = infile[:-3]
    p = subprocess.Popen(['gzip', '-d', '-k', '-q', infile, '-N', outfile], stdout=subprocess.PIPE)
    p.communicate()
    return outfile


def checker(input_arguments, filter_modus, genome_file):
    logging.debug('Checking files')
    files_to_check = {}
    if filter_modus == 'germline':
        files_to_check = {'bam_normal': input_arguments.bam, 'vcf': input_arguments.vcf, 'bed': input_arguments.design, 'genome': genome_file}
    if filter_modus == 'somatic':
        files_to_check = {'bam_normal': input_arguments.bam_normal, 'bam_tumor': input_arguments.bam_tumor,
                          'vcf': input_arguments.vcf, 'bed': input_arguments.design, 'genome': genome_file}

    update_arguments = {}
    for file_type, infile in files_to_check.iteritems():
        file_not_found = check_file(file_type.split('_')[0], os.path.abspath(infile))
        if not file_not_found:
            decomp_file_name = check_extension(file_type.split('_')[0], os.path.abspath(infile))
            if decomp_file_name:
                update_arguments[file_type] = decomp_file_name
            else:
                update_arguments[file_type] = os.path.abspath(infile)
        else:
            break
    update_arguments['filename'] = input_arguments.filename
    update_arguments['outdir'] = os.path.abspath(input_arguments.outdir) + '/'

    return update_arguments, file_not_found
