#!/usr/bin/python

import os, sys, logging

# Count total amount for progression
cmd1 = 'egrep -v "^#" %s | wc -l' % (vcf_path)
p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True)
amount_total_variants, error1 = p1.communicate()

## statistics
offset = dict()
offset['snv'] = {'ref': [], 'alt': []}
offset['ins'] = {'ref': [], 'alt': []}
offset['del'] = {'ref': [], 'alt': []}
offset['sub'] = {'ref': [], 'alt': []}
dp = dict()
dp['snv'] = {'ref': [], 'alt': []}
dp['ins'] = {'ref': [], 'alt': []}
dp['del'] = {'ref': [], 'alt': []}
dp['sub'] = {'ref': [], 'alt': []}

error_get = 0
nr_filtered = 0
nr_reads = 0
nr_filtered_quality = 0
nr_no_match_amp = 0
nr_filtered_one_amp = 0
nr_filtered_normal = 0
nr_offset = 0
nr_variants = 0
nr_fail_nrm = 0
nr_str_offset = 0
ss_ref = 0
ss_germ = 0
ss_som = 0
ss_loh = 0
ss_unkn = 0
nr_fail_pos = 0

# print summary
print "\n"
print "OFFSET FROM GT:"
for vtype in ['snv', 'ins', 'del', 'sub']:
    if len(offset[vtype]['ref']) > 0:
        print "Type: ", vtype, " : "
        print "   Mean REF depth  :", float(sum(dp[vtype]['ref']) / len(dp[vtype]['ref']))
        print "   Mean REF offset :", float(sum(offset[vtype]['ref']) / len(offset[vtype]['ref']))
        print "   Mean ALT depth  :", float(sum(dp[vtype]['alt']) / len(dp[vtype]['alt']))
        print "   Mean ALT offset :", float(sum(offset[vtype]['alt']) / len(offset[vtype]['alt']))

print ""
print "Nr.Variants : ", nr_variants
if somatic_modus == 1:
    print "\tSomatic Status"
    print "\tReference : %s " % (ss_ref)
    print "\tGermline : %s " % (ss_germ)
    print "\tSomatic : %s " % (ss_som)
    print "\tLOH : %s " % (ss_loh)
    print "\tUnknown : %s " % (ss_unkn)
print "Nr.Filtered : ", nr_filtered, "(", (100 * nr_filtered / float(nr_variants)), "% )"
print "\twith one amplicon in design : ", nr_filtered_one_amp, "(", (
100 * nr_filtered_one_amp / float(nr_filtered)), "% )"
print "Nr. one amplicon with matching ref and alt amplicon : ", nr_no_match_amp, "(", (
100 * nr_no_match_amp / float(nr_variants)), "% )"
print "Nr. failed because high read fraction in normal (ONLY APPLIED SOMATIC) : ", nr_fail_nrm, "(", (
100 * nr_fail_nrm / float(nr_variants)), "% )"
print "Nr. failed because read position failed (ONLY APPLIED SOMATIC) : ", nr_fail_pos, "(", (
100 * nr_fail_pos / float(nr_variants)), "% )"
print "\tNr.Reads: ", nr_reads
print "\tNr.Reads below quality: ", nr_filtered_quality, "(", (100 * nr_filtered_quality / float(nr_reads)), "% )"
# print "Nr.BadCounts: ",nr_offset, "(",(100*nr_offset/float(nr_variants)),"% )"
# print "Fraction of badcounts located in repeat:", (100*nr_str_offset/float(nr_offset)),"% )"
print "\n"
