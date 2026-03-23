#!/usr/bin/env bash

###################################################################################################
#                               awkPre                                                            #
#                                                                                                 #
# This script is a required pre-processing step to composing the D-NAscent origins  database.     #
# It sets up the directory structure and creates symbolic links for input data.                   #
#                                                                                                 #
# Usage: saveCommand script/awkPre.bash                                                           #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #
#    08/03/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making symbolic link to DNA replication origins database identified by D-NAscent
ln -s /project/carol/dnascent/project/D-NAscent-Origins/00-originsDataComposition/\
awk/final/synchronized.bed  input/
ln -s /project/carol/dnascent/project/D-NAscent-Origins/00-originsDataComposition/\
awk/final/nonSynchronized.bed input/

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/D-NAscent-Origins/01-windowsDefinition/awk/* .

# Running awk driver script.
saveCommand script/awkDriver.bash 2>&1 | tee log/awkDriver.out

# Checking result.
saveCommand script/awkCheck.bash 2>&1 | tee log/awkCheck.out

# Removing intermediate files.
saveCommand script/awkClean.bash 2>&1 | tee log/awkClean.out
EOI
