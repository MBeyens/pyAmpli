#!/usr/bin/python

import logging


def check_quality_read(read, position, mapping_quality, base_quality):
    try:
        base_qual_index = read.reference_end - position
        if read.mapping_quality >= mapping_quality and read.query_qualities[base_qual_index] >= base_quality:
            pass_quality_read = 1
        else:
            pass_quality_read = 0
    except IndexError: ## PROBLEM OUT OF INDEX OF READ.QUERY_QUALITIES!
        #logging.warning('Quality index error (check_quality_read)')
        pass_quality_read = 0
    except TypeError: ## PROBLEM OUT OF INDEX OF READ.QUERY_QUALITIES!
        #logging.warning('Quality type error (check_quality_read)')
        pass_quality_read = 0

    return pass_quality_read


def get_exact_amplicon(read, design_by_start, design_by_end):
    if read.is_unmapped:
        return ("NA")
    if not str(read.reference_name) in design_by_start:
        return ('NA')
    if str(read.reference_start) in design_by_start[str(read.reference_name)]:
        # check which matches (take the last one if none matches.
        last_end = 0
        for end in design_by_start[str(read.reference_name)][str(read.reference_start)]:
            last_end = end
            if int(end) == int(read.reference_end):
                break
        if last_end != 0:
            return (design_by_start[str(read.reference_name)][str(read.reference_start)][str(last_end)])

    if str(read.reference_end) in design_by_end[str(read.reference_name)]:
        # check which matches (take the last one if none matches.
        last_start = 0
        for start in design_by_end[str(read.reference_name)][str(read.reference_end)]:
            last_start = start
            if int(start) == int(read.reference_start):
                break
        if last_start != 0:
            return (design_by_end[str(read.reference_name)][str(read.reference_end)][str(last_start)])
    return ("NA")
