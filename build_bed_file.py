# -*- coding: utf-8 -*-
"""
This script takes a list of transcript names and generates a sorted BED file.

Usage:
    build_bed_file.py [TRANSCRIPT-NAME-1] [TRANSCRIPT-NAME-2]
"""
import sys

from autoprimer import Gene

def generate_bed(transcripts):
    unsorted_output = list()
    for transcript in transcripts:
        assert len(transcript.split('-')) == 2, transcript + ' is not supported!'
        gene_name, transcript_number = transcript.split('-')
        g = Gene(gene_name, version='GRCh37')
        g.set_transcript(transcript)
        exons = g.list_exon_regions()
        for i, exon in enumerate(exons, 1):
            cds = g.exon_to_translated(exon)
            if cds[1] is not None and cds[2] is not None:
                unsorted_output.append(['chr' + str(cds[0]), cds[1], cds[2], gene_name.upper() + '.EX{exon:03}'.format(exon=i)])
    return unsorted_output

def sort_regions(region_list):
    sorted_regions = sorted(region_list, key=lambda x: x[0].replace('chr', '').zfill(2) + '-' + str(x[1]).zfill(10))
    return sorted_regions

def main():
    if len(sys.argv) < 2:
        print('Usage: build_bed_file.py [TRANSCRIPT-NAME-1]...')
        sys.exit()
    # pass the command line arguments to the generate_bed() function
    unsorted_output = generate_bed(sys.argv[1:])
    sorted_output = sort_regions(unsorted_output)
    for line in sorted_output:
        print('\t'.join([str(x) for x in line]))


if __name__ == '__main__':
    main()
