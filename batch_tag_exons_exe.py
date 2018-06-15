#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 14:19:37 2017

@author: HCC604
"""
import glob
import os
import sys

folder = sys.argv[1]
ABI_files = glob.glob(folder+'/**/*.ab1', recursive=True)

for file in ABI_files:
    print('# Now processing:', file, file=sys.stderr)
    os.system('tag_exons.exe "' + file + '" ' +
              ' '.join(sys.argv[2:]))

