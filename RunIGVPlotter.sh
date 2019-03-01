#!/bin/bash

if [[ $# -ne 7 ]] ; then
    echo 'igv_plotter wrapper script by Tom C.C. Ho (c) 2017-2019'
    echo '-----'
    echo "Usage: RunIGVPlotter.sh [SAMPLE_NAME] [SAMPLE_BAM_PATH] [SAMPLE_VCF_PATH] [CONTROL_BAM_PATH] [UNZIPPED_HG19_PATH] [IGV_PATH] [IGV_PLOTTER_PATH]"
    echo ''
    exit 0
fi

# essential file paths
PATH_TO_IGVPLOTTER=$7
PATH_TO_HG19=$5
PATH_TO_IGV=$6
PATH_TO_TEMP_FILE="./temp.loci"
PATH_TO_CONTROL=$4

# essential parameters

# generate target regions for plotting
grep "^[^#]" $3 | cut -f1,2 --output-delimiter ':' > $PATH_TO_TEMP_FILE

# perform automated IGV plotting
echo "Now plotting sample ${1}..."
$PATH_TO_IGVPLOTTER --igv-jar-path $PATH_TO_IGV -p SAM.SHOW_CENTER_LINE=TRUE -p SAM.ALLELE_THRESHOLD=0.005 -p SAM.ALLELE_USE_QUALITY=FALSE -p SAM.SHOW_SOFT_CLIPPED=TRUE -p SAM.COLOR_BY=READ_STRAND -m 8G --squish -v -L $PATH_TO_TEMP_FILE -o $1 $2 $PATH_TO_CONTROL
echo "Now removing temp file..."
rm $PATH_TO_TEMP_FILE
