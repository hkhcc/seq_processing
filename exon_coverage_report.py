# -*- coding: utf-8 -*-
"""
This script reports the per-position coverage of a gene.

Input -
1. Path to indexed BAM file
2. Transcript name (e.g. ATP7B-201)

Output -
1. Per-exon coverage report.
"""
import os
import re
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt

from autoprimer import Gene, get_flanking_regions

# the maximum coverage in the BAM file (too small a value can result in overflow)
MAX_DEPTH = 1000000
# other custom settings
MIN_ROI_COVERAGE = 20
CDS_FLANK = 20

def cal_coverage(chromosome, start, end, bam_file, max_depth=MAX_DEPTH,min_phred=20, min_mapq=20):
    """Return the output from samtools depth"""
    run_samtools = subprocess.run(['samtools', 'depth',
                                   '-a', 
                                   '-r', 'chr' + str(chromosome) + ':' + str(start) + '-' + str(end),
                                   '-d', str(max_depth),
                                   '-q', str(min_phred),
                                   '-Q', str(min_mapq),
                                    bam_file],
                                    capture_output=True,
                                    text=True)
    return run_samtools.stdout

def parse_coverage(samtools_depth_output, min_depth=20):
    """Return two lists, all positions and failed positions"""
    per_position_result = samtools_depth_output.split('\n')
    all_pos = list()
    failed_pos = list()
    for line in per_position_result:
        line = line.rstrip()
        if line == '':
            continue # skip empty lines
        Chr, base, coverage = re.split(r'\s+', line)
        all_pos.append([Chr, int(base), int(coverage)])
        if int(coverage) < min_depth:
            failed_pos.append([Chr, int(base), int(coverage)])
    return (all_pos, failed_pos)

def main():
    if len(sys.argv) != 3:
        print('Usage: exon_coverage_report.py [BAM_PATH] [TRANSCRIPT_NAME]')
        sys.exit()
    # check that the BAM path exists
    assert os.path.isfile(sys.argv[1]), "BAM_PATH not valid!"
    # assume transcript name is something meaningful
    assert len(sys.argv[2].split('-')) == 2, "Valid transcript name must be in the fomrat [GENE]-[digits]"
    # now retrieve gene info
    gene_name, transcript_number = sys.argv[2].split('-')
    g = Gene(gene_name, version='GRCh37') # we need GRCh37 for hg19 coordinates
    # set the transcript
    g.set_transcript(sys.argv[2])
    # generate the list of coordinates
    exons = g.list_exon_regions()
    cdss = [g.exon_to_translated(exon) for exon in exons]
    rois = list()
    for cds in cdss:
        if cds[1] is None or cds[2] is None:
            rois.append(None)
            continue
        head, tail = get_flanking_regions(cds, flank=CDS_FLANK)
        rois.append([head, cds, tail])

    # determine the coverage for each region
    coverage_report = dict()
    coverage_report['BAM_PATH'] = sys.argv[1]
    coverage_report['CDS_FLANK'] = CDS_FLANK
    x_ticks_pos, x_ticks_label = list(), list()
    for i, roi in enumerate(rois, 1):
        coverage_report['Exon ' + str(i)] = list()
        if roi is None:
            print('Exon', i, 'has nothing to report!', file=sys.stderr)
            continue
        for label, segment in zip(['head', 'body', 'tail'], roi): # each roi has three segments: head, body, tail
            Chr = segment[0]
            start = segment[1]
            end = segment[2]
            coverage = cal_coverage(Chr, start, end, sys.argv[1])
            all_pos, failed_pos = parse_coverage(coverage, min_depth=MIN_ROI_COVERAGE)
            coverage_report['Exon ' + str(i)].append([start, end, all_pos, failed_pos])
            print('Exon', i, label, 'pass rate=', 1 - len(failed_pos)/len(all_pos), file=sys.stderr)
            print('Min. coverage:', min([p[2] for p in all_pos]), file=sys.stderr)
            if len(failed_pos) > 0:
                color = 'red'
            else:
                color = 'cyan'
            if label == 'head':
                plt.bar(i-0.25, np.max([p[2] for p in all_pos]), width=0.25, color='white', edgecolor='black', label="5' flank")
                plt.bar(i-0.25, np.average([p[2] for p in all_pos]), width=0.25, color=color, edgecolor='black', label="5' flank")
                plt.bar(i-0.25, np.min([p[2] for p in all_pos]), width=0.25, color='black', alpha=0.8, edgecolor='black', label="5' flank")
            elif label == 'body':
                plt.bar(i, np.max([p[2] for p in all_pos]), width=0.25, color='white', edgecolor='black', label='CDS')
                plt.bar(i, np.average([p[2] for p in all_pos]), width=0.25, color=color, edgecolor='black', label='CDS')
                plt.bar(i, np.min([p[2] for p in all_pos]), width=0.25, color='black', alpha=0.8, edgecolor='black', label='CDS')
            elif label == 'tail':
                plt.bar(i+0.25, np.max([p[2] for p in all_pos]), width=0.25, color='white', edgecolor='black', label="3' flank")
                plt.bar(i+0.25, np.average([p[2] for p in all_pos]), width=0.25, color=color, edgecolor='black', label="3' flank")
                plt.bar(i+0.25, np.min([p[2] for p in all_pos]), width=0.25, color='black', alpha=0.8, edgecolor='black', label="3' flank")
            x_ticks_pos.append(i)
            x_ticks_label.append('E' + str(i))
    plt.ylabel('Sequencing coverage (folds)')
    plt.yscale('log')
    y_ticks = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    y_ticks.append(MIN_ROI_COVERAGE)
    plt.yticks(ticks=y_ticks, labels=y_ticks)
    plt.xlabel('Exon number')
    plt.xticks(ticks=x_ticks_pos, labels=x_ticks_label)
    plt.axhline(y=MIN_ROI_COVERAGE, color='green', linestyle='--')
    plt.title(sys.argv[2] + ' coverage report\n Sample: ' + os.path.basename(sys.argv[1]) +
              ' Min. depth threshold:' + str(MIN_ROI_COVERAGE) + 'x')
    plt.show()
if __name__ == '__main__':
    main()


