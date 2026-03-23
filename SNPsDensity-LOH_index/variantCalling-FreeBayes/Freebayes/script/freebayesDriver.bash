#!/usr/bin/env bash

###################################################################################################
#                                      freebayesDriver                                            #
#                                                                                                 #
#  Description:                                                                                   #
#  This script acts as the master driver, directing the pre-processing, processing, and           #
#  post-processing steps. It utilizes Freebayes to perform variant calling on Trypanosoma cruzi   #
#  Nanopore reads that do and do not contain DNA replication origins.                             #
#                                                                                                 #
#  Usage:                                                                                         #
#  saveCommand script/freebayesDriver.bash                                                        #
#                                                                                                 #
#  Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                     #
#                       <thiago.franco@fundacaobutantan.org.br>                                   #
#                                                                                                 #
#  25/04/2024: First version.                                                                     #
###################################################################################################
# Loading environment variables about this experiment.
. script/freebayesVariables.bash

# Running awk pre-processing, main and post-processing scripts.
saveCommand script/freebayesPre.bash
saveCommand script/freebayesMain.bash

