#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses exon_coverage_report.py to produce a gene-panel coverage 
report.

Input -
1. Path to indexed BAM file
2. Gene list (tab-separated file, gene name in first column).

Output -
1. Multi-page PDF report in the same directory as the BAM file.
"""
import re
import subprocess
import sys

def main():
    if len(sys.argv) != 3:
        print('Usage: gene_panel_coverage_report.py [BAM_PATH] [GENE_LIST_TXT]')
        sys.exit()
        
    # load the gene list
    genes = []
    with open(sys.argv[2]) as f:
        for line in f:
            if line.startswith('#'):
                continue
            m = re.search(r'([A-Z0-9\-]+)', line)
            if m:
                genes.append(m.group())
    print('# Genes:', genes, file=sys.stderr)
    
    # execute exon_coverage_report
    subprocess.run(['python3', 'exon_coverage_report.py', sys.argv[1]] + genes)

if __name__ == '__main__':
    main()
