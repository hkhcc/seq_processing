#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python 3 script for interactive batch run.

Input
1. Sample folder path
2. Coverage target
3. CDS flank

Output
1. (Automatic running of bwa-vardict-germline.sh on all samples)
 
"""
import os
import sys

PIPELINES = ['bwa-vardict-germline.sh']
PARAMETERS = [['Sample folder path', 'Coverage target', 'CDS flank (bp)']]

script_dir = os.path.dirname(os.path.realpath(__file__))

# check that scripts exists
for script_file in PIPELINES:
    script_file_path = os.path.join(script_dir, script_file)
    print('# Now checking', script_file_path, file=sys.stderr)
    assert os.path.exists(script_file_path), script_file_path + ' does not exist!'

print(script_dir, file=sys.stderr)

print('# Welcome to Multiple Run Manager #', file=sys.stderr)

pipeline_number = 1
for p_and_p in zip(PIPELINES, PARAMETERS):
    print(pipeline_number, ')', p_and_p[0], file=sys.stderr)

choice = -1
while choice <= 0 or choice > len(PIPELINES):
    try:
        choice = int(input('Please choose a pipeline to execute:'))
    except:
        print('Input must be an integer!', file=sys.stderr)

print('Pipeline chosen:', choice, file=sys.stderr)

# now convert the entered value to index of the PIPELINES array
pipeline_index = choice - 1
chosen_pipeline = PIPELINES[pipeline_index]
required_parameters = PARAMETERS[pipeline_index]

para_list = []
for p in required_parameters:
    print('Please enter ', p, ':', sep='', file=sys.stderr)
    input_value = input()
    para_list.append(input_value)

print('# Now generating commands', file=sys.stderr)

# the bwa-vardict-germline.sh command
if pipeline_index == 0:
    commands = []
    run_folder = para_list[0].replace("'", '').replace('"', '').rstrip()
    assert os.path.isdir(run_folder), 'Run folder ' + run_folder + ' does not exist!'
    samples = os.listdir(run_folder)
    for sample in samples:
        command = os.path.join(script_dir, PIPELINES[0]) + ' ' # the pipeline script path
        command += sample + ' ' # the sample name
        # now we need to determine the path of the fastq files
        files_in_dir = os.listdir(os.path.join(run_folder, sample))
        fastq1, fastq2, gene_list = '', '', ''
        for file in files_in_dir:
            if file.endswith('1.fq.gz'):
                fastq1 = os.path.join(run_folder, sample, file)
            if file.endswith('2.fq.gz'):
                fastq2 = os.path.join(run_folder, sample, file)
            if file.endswith('.txt'):
                gene_list = os.path.join(run_folder, sample, file)
        print('Finished loading directory content.', file=sys.stderr)
        if fastq1 == '' or fastq2 == '' or gene_list == '':
            print('Could not retrieve the required files!', file=sys.stderr)
            print('Please check', os.path.join(run_folder, sample), file=sys.stderr)
        command += "'" + fastq1 + "' "
        command += "'" + fastq2 + "' "
        command += para_list[1] + " "
        command += para_list[2] + " "
        command += "'" + gene_list + "'"
        commands.append(command)

for command in commands:
    print(command, file=sys.stderr)
    os.system(command)




