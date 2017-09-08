#!/usr/bin/python

import sys, logging, time


def update_progress(variant_number, amount_total_variants):
    progress = (float(variant_number) / float(amount_total_variants))
    barLength = 100  # Modify this to change the length of the progress bar
    status = "- running ({0}/{1})- ".format(variant_number, amount_total_variants)
    if isinstance(progress, int):
        progress = float(progress)
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    elif progress >= 1:

        time.sleep(3)
        block = int(round(barLength * 1))
        text = "\r --> Status : [{0}] {1}% {2} \n".format("#" * block + "-" * (barLength - block),
                                                                               round((1 * 100), 2), "finished")
        sys.stdout.write(text)
        sys.stdout.flush()
        logging.info('VCF processing done')
    else:
        block = int(round(barLength * progress))
        text = "\r --> Status : [{0}] {1}% {2} ".format("#" * block + "-" * (barLength - block),
                                                                               round((progress * 100), 2), status)
        sys.stdout.write(text)
        sys.stdout.flush()