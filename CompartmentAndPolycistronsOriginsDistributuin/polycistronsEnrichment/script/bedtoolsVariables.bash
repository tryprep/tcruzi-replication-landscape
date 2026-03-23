###################################################################################################
#                                          bedtoolsVariables                                      #
#                                                                                                 #
# This script contains all the variables needed in this computational analysis.                   #
#                                                                                                 #
# Usage: . script/bedtoolsVariables.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #    
#    08/07/2023: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="100M"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=2

# Genome file that will be used in the annotation
export GENOME="input/genome.gff"








