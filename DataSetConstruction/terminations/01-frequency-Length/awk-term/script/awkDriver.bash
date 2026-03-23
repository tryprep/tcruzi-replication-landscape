#!/usr/bin/env bash

###################################################################################################
#                                            awkDriver                                            #
#                                                                                                 #
#    This script will direct the other pre-processing, processing, and post-processing scripts.   #
# This script In this script, we'll use awk to build databases of DNA replication origins for     #
# synchronized and unsynchronized cells in the S phase of the cycle, as identified by the         #
# "D-Nascent" technique.                                                                          #
#                                                                                                 #
# Usage: saveCommand script/awkDriver.bash                                                        #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #
#    08/03/2023: First version.                                                                   #
###################################################################################################

# Loading environment variables about this experiment.
. script/awkVariables.bash

# Running awk pre-processing, main and post-processing scripts.
saveCommand script/awkPre.bash
saveCommand script/awkMain.bash
#saveCommand script/awkPost.bash
