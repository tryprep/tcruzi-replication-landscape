###################################################################################################
#                                          awkVariables                                           #
#                                                                                                 #
#   This script contains all the variables needed in this computational analysis.                 #
#                                                                                                 #
# Usage: . script/awkVariables.bash                                                               #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                     <thiago.franco.esib@esib.butantan.gov.br>                                   #
#    08/03/2023: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="30G"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=30

# bed file from D-NAscet database
export WS=3000
