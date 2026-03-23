#!/usr/bin/env bash

###################################################################################################
#                                            cppDriver                                            #
#                                                                                                 #
#    This script will direct the other pre-processing, processing, and post-processing scripts.   #
# This script In this script, we'll use awk to build databases of DNA replication origins for     #
# synchronized and unsynchronized cells in the S phase of the cycle, as identified by the         #
# "D-Nascent" technique.                                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.franco@fundacaobutantan.org.br >                                   #
#                       David da Silva Pires                                                      #
#                      <david.pires@fundacaobutantan.org.br >                                     #
#    13/09/2024: First version.                                                                   #
###################################################################################################

# Loading environment variables about this experiment.
. script/cppVariables.bash

# Running awk pre-processing, main and post-processing scripts.
saveCommand script/cppPre.bash
saveCommand script/cppMain.bash
