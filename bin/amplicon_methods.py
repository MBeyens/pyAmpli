#!/usr/bin/python

import os, sys, logging

def Get_Exact_Amplicon(read):
    # print read.reference_name , ":" , read.reference_start, "-",read.reference_end
    if read.is_unmapped: return("NA")
    if not str(read.reference_name) in design_by_start:
        return('NA')
    if str(read.reference_start) in design_by_start[str(read.reference_name)]:
        # check which matches (take the last one if none matches.
        last_end = 0
        for end in design_by_start[str(read.reference_name)][str(read.reference_start)]:
            last_end = end
            if int(end) == int(read.reference_end):
                break
        if last_end != 0:
            return( design_by_start[str(read.reference_name)][str(read.reference_start)][str(last_end)])

    if str(read.reference_end) in design_by_end[str(read.reference_name)]:
        # check which matches (take the last one if none matches.
        last_start = 0
        for start in design_by_end[str(read.reference_name)][str(read.reference_end)]:
            last_start = start
            if int(start) == int(read.reference_start):
                break
        if last_start != 0:
            return( design_by_end[str(read.reference_name)][str(read.reference_end)][str(last_start)])
    return("NA")