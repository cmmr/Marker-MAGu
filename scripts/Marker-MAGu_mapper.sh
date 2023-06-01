#!/bin/bash

### description:
### Takes metagenomic short reads, e.g. from Illumina machine
### optional QC and removal of host/spike-in reads
### alignment to marker gene database of phage, bacteria, archaea, and micro-eukaryotes (Metaphlan4 + Trove of Gut Virus Genomes)
### Makes taxonomical profile


### arguments
READS=$1
SAMPLE=$2
CPUS=$3
OUT_DIR=$4
QUAL=$5
FILTER_SEQS=$6
TEMP_DIR=$7
KEEP=$8
MM_DB=$9
MARKERMAGU_DIR=${10}

MDYT=$( date +"%m-%d-%y---%T" )
echo "Time Update: Starting main bash mapper script for Marker-MAGu @ $MDYT"

#arguments check
if [ $# -ne 9 ] ; then 
    echo "expected 9 arguments passed on the command line:"
    echo "read file(s), sample ID, CPUs, output directory, trim by quality?, filter seqs?, temp directory path, keep temp?, Marker-MAGu script directory"
    echo "exiting"
    exit
fi


## temporary directory variable setup
if [ "$TEMP_DIR" == "default" ] ; then
    TEMP_DIR="${OUT_DIR}/${SAMPLE}_temp"

elif [ -d $TEMP_DIR ] ; then
    MDYT_g=$( date +"%m-%d-%y---%T" | sed 's/:/_/g' )
    TEMP_DIR="${TEMP_DIR}/${SAMPLE}_temp_${MDYT_g}"

fi

## check filter_seqs
if [ "$FILTER_SEQS" == "True" ] && [ ! -s ${MARKERMAGU_DIR%scripts}filter_seqs/filter_seqs.fna ]; then
    echo "-f True flag requires that this file exists and is not empty: "
    echo "${MARKERMAGU_DIR%scripts}filter_seqs/filter_seqs.fna"
    echo "exiting"
    exit
fi

## check database
if [ "$MM_DB" == "default" ] ; then
    DB_DIR=$( find ${MARKERMAGU_DIR%scripts}DBs/ -type d -name "v*" | tail -n1 )
else
    DB_DIR=${MARKERMAGU_DIR%scripts}DBs/${MM_DB}
fi

if [ ! -d ${DB_DIR} ] ; then
    echo "can't find DB directory. should be a versioned directory here: ${MARKERMAGU_DIR%scripts}DBs/ "
    echo "exiting"
    exit
fi

if [ ! -s ${DB_DIR}/Marker-MAGu_markerDB.fna ] ; then
    echo "can't find marker database for pipeline. should be: ${DB_DIR}/Marker-MAGu_markerDB.fna"
    echo "exiting"
    exit
fi


## check output directories
if [ ! -d ${OUT_DIR} ] ; then
    mkdir ${OUT_DIR}
fi
if [ ! -d ${OUT_DIR}/record/ ] ; then
    mkdir ${OUT_DIR}/record/
fi
if [ ! -d ${TEMP_DIR} ] ; then
    mkdir ${TEMP_DIR}
fi

## arguments file
date > ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Marker-MAGu used arguments" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Sample name:                  $SAMPLE" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Read file(s):                 $READS" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "CPUs:                         $CPUS" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Output directory:             $OUT_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Trim for quality:             $QUAL" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Remove host/spikein seqs:     $FILTER_SEQS" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Temp directory path:          $TEMP_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Keep temp files:              $KEEP" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Marker-MAGu script directory: $MARKERMAGU_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Marker-MAGu database used:    ${DB_DIR}/Marker-MAGu_markerDB.fna" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt


cat ${OUT_DIR}/record/${SAMPLE}.arguments.txt

## Check if reads exist
for READ_FILE in $READS ; do
    if [ -s $READ_FILE ] ; then
        echo $READ_FILE
    else
        echo "$READ_FILE not found."
        echo "exiting"
        exit
    fi
done

## Running quality trimming/host filtering
if [ "$QUAL" == "True" ] && [ "$FILTER_SEQS" == "True" ] ; then
    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: Trimming low quality reads with fastp THEN Aligning reads to filter_seqs.fna to remove host/spike-in @ $MDYT"

    ##trim
    ##filter
    cat ${READS} | \
    fastp --stdin --stdout -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.json | \
    minimap2 -t $CPUS -ax sr ${MARKERMAGU_DIR%scripts}filter_seqs/filter_seqs.fna - | \
    samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.MM_input.fastq

elif [ "$QUAL" == "True" ] ; then
    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: Trimming low quality reads with fastp @ $MDYT"

    ##trim
    cat ${READS} | \
    fastp --stdin -o ${TEMP_DIR}/${SAMPLE}.MM_input.fastq -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.html

elif [ "$FILTER_SEQS" == "True" ] ; then
    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: Aligning reads to filter_seqs.fna to remove host/spike-in @ $MDYT"

    ##filter
    cat ${READS} | 
    minimap2 -t $CPUS -ax sr ${MARKERMAGU_DIR%scripts}filter_seqs/filter_seqs.fna - | samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.MM_input.fastq

else
    ## cat
    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: Concatenating input reads @ $MDYT"

    cat ${READS} > ${TEMP_DIR}/${SAMPLE}.MM_input.fastq
fi


if [ -s ${TEMP_DIR}/${SAMPLE}.MM_input.fastq ] ; then

    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: running seqkit stats on ${SAMPLE} @ $MDYT"

    seqkit stats -T ${TEMP_DIR}/${SAMPLE}.MM_input.fastq > ${OUT_DIR}/${SAMPLE}.MM_input.seq_stats.tsv
    FILTERED_READS=$( tail -n1 ${OUT_DIR}/${SAMPLE}.MM_input.seq_stats.tsv | cut -f4  )

    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: running minimap2 on ${SAMPLE} @ $MDYT"

    minimap2 -t $CPUS -ax sr ${DB_DIR}/Marker-MAGu_markerDB.fna --split-prefix ${TEMP_DIR}/Marker-MAGu_markerDB ${TEMP_DIR}/${SAMPLE}.MM_input.fastq > ${TEMP_DIR}/${SAMPLE}.markermagu.sam

    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: running samtools on ${SAMPLE} @ $MDYT"

    samtools view -@ $CPUS -bSq 1 ${TEMP_DIR}/${SAMPLE}.markermagu.sam | samtools sort -@ $CPUS -o ${TEMP_DIR}/${SAMPLE}.markermagu.sort.bam

    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: running coverm on ${SAMPLE} @ $MDYT"

    coverm contig --bam-files ${TEMP_DIR}/${SAMPLE}.markermagu.sort.bam --min-read-percent-identity 90 --min-read-aligned-percent 50 -m length covered_bases count -o ${TEMP_DIR}/${SAMPLE}.marker-magu.unique_alignment.coverm.tsv -t $CPUS

    MDYT=$( date +"%m-%d-%y---%T" )
    echo "Time Update: running treshold enforcer/abundance calculator Rscript on ${SAMPLE} @ $MDYT"

    Rscript ${MARKERMAGU_DIR}/make_abundance_table.R ${TEMP_DIR}/${SAMPLE}.marker-magu.unique_alignment.coverm.tsv $FILTERED_READS ${OUT_DIR}/${SAMPLE} ${SAMPLE}

else
    echo "${TEMP_DIR}/${SAMPLE}.MM_input.fastq not found"
fi

## delete temp
if [ "$KEEP" == "True" ] ; then
    echo "Keeping temporary files in: ${TEMP_DIR}"
else
    echo "Removing temp files"
    rm ${TEMP_DIR}/*
fi

## finish
echo "Main Output Files Generated in ${OUT_DIR}/"

find ${OUT_DIR}/ -type f -name "${SAMPLE}.detected_species.tsv"
find ${OUT_DIR}/ -type f -name "${SAMPLE}.MM_input.seq_stats.tsv"


MDYT=$( date +"%m-%d-%y---%T" )
echo "Time Update: Finishing Marker-MAGu ${SAMPLE} @ $MDYT"

echo "##################"
echo "##################"
echo "##################"
echo " "
