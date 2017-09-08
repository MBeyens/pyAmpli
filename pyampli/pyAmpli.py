#!/usr/bin/python

###
# Load python dependencies
import os, sys, time
###
# Load script dependencies
from pyampli import log_system, args, config, check_input_files, get_design, make_vcf, process_bam, variant_filter


def main():
    t0 = time.time()
    ###
    # Initiate logging systems
    debug_mode = 0
    logging = log_system.initiate_log(debug_mode)
    logging.debug('Log initiated')

    ###
    # Argument parsing
    input_arguments = args.argument_checkup()
    print input_arguments.config
    
    ###
    # Load config
    config_parameters = config.config_load(input_arguments.config)
    logging.debug('Config loaded')

    ###
    # Create directory structure
    if input_arguments.outdir is None:
        results_dir = os.getcwd() + '/amplicon_filter/'
        tmp_results_dir = os.getcwd() + '/amplicon_filter/tmp/'
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        if not os.path.exists(tmp_results_dir):
            os.makedirs(tmp_results_dir)
    else:
        results_dir = os.path.abspath(input_arguments.outdir) + '/'
        tmp_results_dir = os.path.abspath(input_arguments.outdir) + '/tmp/'
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        if not os.path.exists(tmp_results_dir):
            os.makedirs(tmp_results_dir)
    logging.info('Results directory created : %s', results_dir)

    ###
    # Check filter modus
    try:
        input_arguments.bam_tumor
        filter_modus = 'somatic'
        logging.info('Selected somatic modus')
    except AttributeError:
        filter_modus = 'germline'
        logging.info('Selected germline modus')

    ###
    # Check input argument: exist and decompress
    input_arguments, file_not_found = check_input_files.checker(input_arguments, filter_modus, config_parameters['reference'])
    print "FILEEEE::", file_not_found

    if not file_not_found:
        ###
        # Load design file
        design_by_start, design_by_end = get_design.load_design_file(input_arguments['bed'])

        ###
        # Check design conflicts
        # TO DO
        # fuzziness = 0
        # get_design.check_design_conflicts(config_parameters['general_settings']['fuzziness'], design_by_start)

        ###
        # Create new vcf template
        input_arguments = make_vcf.create_vcf(filter_modus, input_arguments['filename'],
                                              input_arguments['outdir'], input_arguments['vcf'],
                                              input_arguments)

        ###
        # Get sequencing lengths from BAM(s)
        seq_lengths = process_bam.bam_loader(input_arguments)

        ###
        # Open unfiltered vcf file
        unfiltered_vcf_file = make_vcf.open_unfiltered_vcf(input_arguments['vcf'])

        ###
        # Filter unfiltered vcf file
        variant_number = variant_filter.multi_processing_variants(filter_modus, input_arguments, config_parameters,
                                                                  unfiltered_vcf_file, design_by_start, design_by_end,
                                                                  seq_lengths)

        ###
        # Merge multiprocessed files
        variant_filter.merge_tmp_vcfs(input_arguments, variant_number)

        ###
        # Summary statistics
        # get_summary_statistics

        ###
        # End analysis
        logging.info('pyAMPLI finshed!')
        t1 = time.time()
        print ("Total time running: %s seconds" % (str(t1-t0)))

    else:
        logging.info(file_not_found)
        logging.info('Fix file location or file name, and restart pyAmpli!')


if __name__ == "__main__":
    main()
