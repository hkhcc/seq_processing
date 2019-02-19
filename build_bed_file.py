# -*- coding: utf-8 -*-
"""
This script takes a list of transcript names and generates a sorted BED file.

Usage:
    build_bed_file.py [TRANSCRIPT-NAME-1] [TRANSCRIPT-NAME-2]
"""
import argparse

from autoprimer import Gene

def generate_bed(transcripts, flanking=10):
    unsorted_output = list()
    for transcript in transcripts:
        if len(transcript.split('-')) == 2:
            # user-specified transcript
            gene_name, transcript_number = transcript.split('-')
            g = Gene(gene_name, version='GRCh37')
            g.set_transcript(transcript)
        else:
            # use canonical transcript from Ensembl
            gene_name = transcript
            g = Gene(gene_name, version='GRCh37')
            g.set_transcript()
        exons = g.list_exon_regions()
        for i, exon in enumerate(exons, 1):
            cds = g.exon_to_translated(exon)
            if cds[1] is not None and cds[2] is not None:
                assert cds[1] < cds[2], 'cds[1] is not less than cds[2]!'
                flanked_cds = [cds[1] - flanking, cds[2] + flanking]
                unsorted_output.append(['chr' + str(cds[0]),
                                        flanked_cds[0], flanked_cds[1],
                                        gene_name.upper() + '.EX{exon:03}'.format(exon=i) +
                                        '.' + str(cds[1]) + '.' + str(cds[2])
                                        ]
                                        )
    return unsorted_output

def sort_regions(region_list):
    sorted_regions = sorted(region_list, key=lambda x: x[0].replace('chr', '').zfill(2) + '-' + str(x[1]).zfill(10))
    return sorted_regions

def main():
    parser = argparse.ArgumentParser(description='Builds a BED file from a list of provided genes/ transcripts.')
    parser.add_argument('transcripts', nargs='+', help='gene name or trancript, e.g. ATP7B or ATP7B-201')
    parser.add_argument('-f', dest='cds_flank', type=int, default=10, help='flanking region to add to each target interval (default: 10)')
    args = parser.parse_args()
    # pass the command line arguments to the generate_bed() function
    unsorted_output = generate_bed(args.transcripts, flanking=args.cds_flank)
    sorted_output = sort_regions(unsorted_output)
    for line in sorted_output:
        print('\t'.join([str(x) for x in line]))


if __name__ == '__main__':
    main()
