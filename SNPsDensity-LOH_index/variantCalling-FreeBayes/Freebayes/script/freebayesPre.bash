#!/usr/bin/env bash

###################################################################################################
#                                          freebayesPre                                           #
#                                                                                                 #
#  Description:                                                                                   #
#  This script performs the essential pre-processing steps for the variant calling pipeline       #
#  of Trypanosoma cruzi (CL Brener strain) Nanopore genomic reads. It handles directory           #
#  structuring, parallel creation of input/output folders, and symlinking of raw datasets         #
#  and reference metadata.                                                                        #
#                                                                                                 #
#  Usage:                                                                                         #
#  saveCommand script/freebayesPre.bash                                                           #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.francofundacaobutantan.org.br>                                     #
#    25/04/2024: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir -p  final input job output

# Storing initial processing date and time.
date > log/processingStart.date

# Creating symbolic links to all reads after basecalling
ln -s /project/carol/dnascent/project/tryCru-clb1/00-baseCalling/\
guppy/nonSinc/1/final/*.fastq input/1/
ln -s /project/carol/dnascent/project/tryCru-clb1/00-baseCalling/\
guppy/nonSinc/2/final/*.fastq input/2/
ln -s /project/carol/dnascent/project/tryCru-clb1/00-baseCalling/\
guppy/nonSinc/3/final/*.fastq input/3/
ln -s /project/carol/dnascent/project/tryCru-clb1/00-baseCalling/\
guppy/nonSinc/4/final/*.fastq input/4/
ln -s /project/carol/dnascent/project/tryCru-clb1/00-baseCalling/\
guppy/nonSinc/5/final/*.fastq input/5/
ln -s /project/carol/dnascent/project/tryCru-clb1/00-baseCalling/\
guppy/nonSinc/negControl/final/*.fastq input/negControl/

# Creating symbolic link to reference genome file (fasta)
ln -s /project/carol/dnascent/project/D-NAscent-Origins/\
metaData/tryCru-clb7.fasta   input/genome.fa

# Creating symbolic link to reference genome file (gff)
ln -s /project/carol/dnascent/project/D-NAscent-Origins/\
metaData/tryCru-clb7.gff    input/genome.gff

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/dnascentOriginsProcessing/05-variantCalling/freebayes/* .

# Running freebayes driver script.
saveCommand script/freebayesDriver.bash 2>&1 | tee log/freebayesDriver.out

EOI

exit 0
