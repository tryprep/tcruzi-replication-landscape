###################################################################################################
#                                          freebayesVariables                                     #
#                                                                                                 #
#   This script contains all the variables needed in this computational analysis.                 #
#                                                                                                 #
# Usage: . script/freebayesVariables                                                              #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.francofundacaobutantan.org.br>                                     #
#    25/04/2024: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="100G"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=30

# Reference genome used in this analysis.
export GENOME_FASTA=input/genome.fa

export GENOME_GFF=input/genome.gff

# Directories for the processed files of each process.
export SAMPLES="1 2 3 4 5 negControl"

# Dataset of origin with working.
export TYPE_ORIGIN="nonSynchronized"

# Directories to input and output.
export INPUT_DIR="input"
export OUTPUT_DIR="output"
