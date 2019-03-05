#!/usr/bin/env bash

echo '# Starting PYNEH NGS analysis pipeline version 20190228...'

# Show the help message if the required number of arguments is not found
if [[ $# -ne 6 ]]; then
	echo "    Usage: $0 [SAMPLE_NAME] [FASTQ1] [FASTQ2] [COVERAGE] [FLANK_BP] [GENE_LIST_TXT]"
	exit
fi

# Check the required programs
echo '====='
echo "[Timestamp: `date`]"
echo '# Initializing...'
SCRIPT_DIR=`dirname "$(readlink -f "$0")"`
echo "    Script directory: $SCRIPT_DIR"
echo '== Programs =='
echo "    FastQC path: `which fastqc`"
echo "    BWA path: `which bwa`"
echo "    SAMTOOLS path: `which samtools`"
echo "    VARDICT (core) path: `which VarDict`"
echo "    VARDICT (strand bias test) path: `which teststrandbias.R`"
echo "    R path: `which R`"
echo "    VARDICT (var2vcf) path: `which var2vcf_valid.pl`"
echo "    ANNOVAR path: `which table_annovar.pl`"
echo "    PYTHON3 path: `which python3`"
BUILD_BED_FILE=`readlink -e $SCRIPT_DIR/build_bed_file.py`
echo "    BUILD_BED_FILE path: $BUILD_BED_FILE"
IGV_PATH=`readlink -e ~/Programs/IGV/IGV_2.4.19/lib/igv.jar`
echo "    IGV path: $IGV_PATH"
IGV_PLOTTER_PATH=`which igv_plotter`
echo "    IGV_PLOTTER path: $IGV_PLOTTER_PATH"
RUN_IGV_PLOTTER=`readlink -e $SCRIPT_DIR/RunIGVPlotter.sh`
echo "    RUN_IGV_PLOTTER path: $RUN_IGV_PLOTTER"
PDF_REPORT_TOOL=`readlink -e $SCRIPT_DIR/ngs_variant_report.py`
echo "    NGS_VARIANT_REPORT path: $PDF_REPORT_TOOL"
PDF_COVREPORT_TOOL=`readlink -e $SCRIPT_DIR/gene_panel_coverage_report.py`
echo "    NGS_COVERAGE_REPORT path: $PDF_COVREPORT_TOOL"
echo '== Databases =='
HG19_PATH=`readlink -e ~/bundle/hg19/ucsc.hg19.fasta.gz`
if [[ -r $HG19_PATH ]]; then
	echo '    hg19... OK'
else
	echo '    hg19... Failure!'
	exit
fi

UNZIP_HG19_PATH=`readlink -e ~/bundle/hg19/ucsc.hg19.fasta`
if [[ -r $UNZIP_HG19_PATH ]]; then
	echo '    hg19 (unzip)... OK'
else
	echo '    hg19 (unzip)... Failure!'
	exit
fi

HUMANDB_PATH=`readlink -e ~/humandb`
if [[ -r $HUMANDB_PATH ]]; then
	echo '    humandb... OK'
else
	echo '    humandb... Failure!'
	exit
fi

CONTROL_BAM_PATH=`readlink -e ~/localdb/control.bam`
if [[ -r $CONTROL_BAM_PATH ]]; then
	echo '    control BAM file... OK'
else
	echo '    control BAM file... Failure!'
	exit
fi

# Check the input files
echo '# Validating inputs...'
echo "    Sample name: $1"
echo "    FASTQ 1: $2"
if [[ -r $2 ]]; then
	echo "    ..File check... OK"
	FASTQ1=$2
else
	echo "    ..File check... Failure!"
	exit
fi

echo "    FASTQ 2: $3"
if [[ -r $3 ]]; then
	echo "    ..File check... OK"
	FASTQ2=$3
else
	echo "    ..File check... Failure!"
	exit
fi

DIR_A=`dirname $2`
DIR_B=`dirname $3`
if [[ "$DIR_A" != "$DIR_B" ]]; then
	echo "    FASTQ files in different directories!"
	exit
else
	echo "# Setting output directory to: $DIR_A"
	OUTPUT_DIR=$DIR_A
fi

COVERAGE=$4
echo "    Minimum coverage set to: $COVERAGE x"

FLANK_BP=$5
echo "    CDS flanking set to: $FLANK_BP bp"

echo "    Gene list text file: $6"
if [[ -r $6 ]]; then
	echo "    ..File check... OK"
	GENE_LIST_TXT=$6
else
	echo "    ..File check... Failure!"
	exit
fi

# Perform mapping
BAM_PATH="$OUTPUT_DIR/$1.bam"
BAI_PATH="$OUTPUT_DIR/$1.bam.bai"
if [ -r $BAM_PATH ] && [ -r $BAI_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! Files found. Skipping Step 1. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 1: Mapping...'
	echo "# Writing output to $BAM_PATH"
	`bwa mem -R "@RG\tID:$1\tSM:$1\tLB:$1\tPL:ILMN_LIKE" -t 8 $HG19_PATH $FASTQ1 $FASTQ2 | samtools sort -@8 -o $BAM_PATH`
	echo "# Indexing $BAM_PATH"
	`samtools index -@8 $BAM_PATH $BAI_PATH`
fi

# Generate BED file
BED_PATH="$OUTPUT_DIR/$1.bed"
if [ -r $BED_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! File found. Skipping Step 2. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 2: Generate BED file for variant calling...'
	echo "# Writing output to $BED_PATH"
	GENE_STRING=`grep -v '^#'  $GENE_LIST_TXT | grep -v '^$' | tr '\r\n' '\n' | tr '\n' ' '`
	`python3 $BUILD_BED_FILE -f $FLANK_BP $GENE_STRING > $BED_PATH`
fi

# Perform variant calling
RAW_VCF_PATH="$OUTPUT_DIR/$1.raw.vcf"
if [ -r $RAW_VCF_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! File found. Skipping Step 3. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 3: Variant calling...'
	echo "# Writing output to $RAW_VCF_PATH"
	`VarDict -G $UNZIP_HG19_PATH -f 0.05 -I 1000 -L 1001 -k 1 -N $1 -th 8 --dedup -b $BAM_PATH -z -c 1 -S 2 -E 3 -g 4 $BED_PATH | teststrandbias.R | var2vcf_valid.pl -N $1 -E -f 0.05 > $RAW_VCF_PATH`
fi

# Perform variant annotation
ANNO_VCF_PREFIX="$OUTPUT_DIR/$1.annovar"
ANNO_VCF_PATH="$OUTPUT_DIR/$1.annovar.hg19_multianno.vcf"
if [ -r $ANNO_VCF_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! File found. Skipping Step 4. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 4: Variant annotation...'
	echo "# Writing output to $ANNO_VCF_PATH"
	`table_annovar.pl $RAW_VCF_PATH $HUMANDB_PATH -buildver hg19 -out $ANNO_VCF_PREFIX -remove -protocol refGene,ensGene,1000g2015aug_all,1000g2015aug_eas,gnomad_genome,exac03,gnomad_exome,dbnsfp35c,dbscsnv11,intervar_20180118,clinvar_20180603,avsnp150 -operation g,g,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput`
fi

# Perform IGV plotting

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 5: Automatic IGV plotting...'
echo "# Writing output to $OUTPUT_DIR/snapshots"

cd $OUTPUT_DIR

if [ -r "$OUTPUT_DIR/snapshots" ]; then
	echo "IGV snapshot directory found."
else
	echo "Creating IGV snapshot directory..."
	mkdir $OUTPUT_DIR/snapshots
fi

cd "$OUTPUT_DIR/snapshots"

echo "Current working directory: `pwd`"
$RUN_IGV_PLOTTER $1 $BAM_PATH $ANNO_VCF_PATH $CONTROL_BAM_PATH $UNZIP_HG19_PATH $IGV_PATH $IGV_PLOTTER_PATH
cd $SCRIPT_DIR
echo "Current working directory: `pwd`"

# Compile PDF variant report

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 6: Compile PDF variant report'
echo "# Output will be written to $OUTPUT_DIR/$1.pdf"
python3 $PDF_REPORT_TOOL $ANNO_VCF_PATH
if [ -r "$ANNO_VCF_PATH.pdf" ]; then
	echo 'Renaming output...'
	mv "$ANNO_VCF_PATH.pdf" "$OUTPUT_DIR/$1.pdf"
else
	echo 'Output file not found!'
	exit
fi

# Compile PDF coverage report

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 7: Compile PDF coverage report'
echo "# Output will be written to $OUTPUT_DIR/$1.coverage.pdf"
python $PDF_COVREPORT_TOOL $BAM_PATH $COVERAGE $FLANK_BP $GENE_LIST_TXT
if [ -r "$BAM_PATH.coverage.pdf" ]; then
	echo 'Renaming output...'
	mv "$BAM_PATH.coverage.pdf" "$OUTPUT_DIR/$1.coverage.pdf"
else
	echo 'Output file not found!'
	exit
fi

# Generate quality statistics

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 8: Generate FastQC reports'
echo "# Output will be written to $OUTPUT_DIR/fastqc"
if [ -d "$OUTPUT_DIR/fastqc" ]; then
	echo 'FastQC output directory found. Step 8 will be skipped.'
	echo "Note: Remove $OUTPUT_DIR/fastqc in order to trigger the FastQC step."
else
	cd $OUTPUT_DIR
	echo "Current working directory: `pwd`"
	fastqc $FASTQ1 $FASTQ2
	mkdir fastqc
	mv *.html *.zip ./fastqc
	cd $SCRIPT_DIR
	echo "Current working directory: `pwd`"
fi


