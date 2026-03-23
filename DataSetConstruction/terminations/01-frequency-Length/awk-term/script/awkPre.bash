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

## First, we will bring the D-NAscent files from the project/tryCru-clb1##
# Making symbolic link to Non-Sincronized files
ln -s /project/carol/dnascent/project/D-NAscent-Origins/00-originsDataComposition/\
awk-term/final/tableOriginsNonSinc.bed input/nonsynchronized.bed

# Making symbolic link to Syncronized Files.
ln -s /project/carol/dnascent/project/D-NAscent-Origins/00-originsDataComposition/\
awk-term/final/tableOriginsSinc.bed input/syncronized.bed

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/D-NAscent-Origins/00-originsDataComposition/awk/* .

# Running awk driver script.
saveCommand script/awkDriver.bash 2>&1 | tee log/awkDriver.out

# Checking result.
saveCommand script/awkCheck.bash 2>&1 | tee log/awkCheck.out

# Removing intermediate files.
saveCommand script/awkClean.bash 2>&1 | tee log/awkClean.out
EOI
