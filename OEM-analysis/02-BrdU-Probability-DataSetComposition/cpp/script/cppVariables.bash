###################################################################################################
#                                          cppVariables                                           #
#                                                                                                 #
#   This script contains all the variables needed in this computational analysis.                 #
#                                                                                                 #
# Usage: . script/cppVariables.bash                                                               #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.franco@fundacaobutantan.org.br >                                   #
#                       David da Silva Pires                                                      #
#                      <david.pires@fundacaobutantan.org.br >                                     #
#    23/09/2024: First version.                                                                   #
###################################################################################################

# Number of memory necessary for SLURM
export MEMORY_SIZE="40G"

# Number of CPUs necessary to run job on SLURM
export NUM_OF_CPUS="4"  # Total number of available CPUs
export NUM_OF_CPUS_PARALLEL="10"  # Number of simultaneous tasks
export THREADS_PER_TASK=$((NUM_OF_CPUS / NUM_OF_CPUS_PARALLEL))  # Threads per task


# Genome fasta file
export GENOME_FASTA="input/genome.fa"
export GENOME_GFF="input/genome.gff"

# Subdirectories to samples
export SAMPLES="1 2 3 4 5 6 7 8"
