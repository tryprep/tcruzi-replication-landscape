#!/usr/bin/env bash

###################################################################################################
#                                 bedtoolsPre                                                     #
#                                                                                                 #
# This script is a required pre-processing step to annotation the D-NAscent origins database.     #
# It sets up the directory structure and creates symbolic links for input data.                   #                                                                                                 #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsPre.bash                                                      #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                       <thiago.franco.esib@esib.butantan.gov.br>                                 #
#    08/07/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data
# Synchronized and nonSynchronized origins dataset.
ln -s /project/carol/dnascent/project/D-NAscent-Origins/02-windowsDefinition/\
awk/final/nonSynchronized.bed  input/
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/02-windowsDefinition/\
#awk/final/synchronized.bed  input/

# Frequency origins dataset.

## Synchronized.
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03b-frequencyClassification/\
#awk/final/synchronized-lessFrequent.bed input/
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03b-frequencyClassification/\
#awk/final/synchronized-moreFrequent.bed input/

## Nonsynchronized.
ln -s /project/carol/dnascent/project/D-NAscent-Origins/03b-frequencyClassification/\
awk/final/nonSynchronized-lessFrequent.bed input/
ln -s /project/carol/dnascent/project/D-NAscent-Origins/03b-frequencyClassification/\
awk/final/nonSynchronized-moreFrequent.bed input/

# Orc1Cdc6 free origins dataset
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03a-Orc1Cdc6Free/bedtools/\
#final/orc1Cdc6-free-nonSynchronized.bed input/
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03a-Orc1Cdc6Free/bedtools/\
#final/orc1Cdc6-free-synchronized.bed input/

# Nonsynchronized classified late and erly origins.
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03c-earlyLaterClassification/\
#bedtools/final/early.bed input/
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03c-earlyLaterClassification/\
#bedtools/final/late.bed input/

#  Terminations
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/02-windowsDefinition/awk-term/\
#final/nonSynchronized-termination.bed input/
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03b-frequencyClassification/awk-term/\
#final/nonSynchronized-lessFrequent.bed  input/
#ln -s /project/carol/dnascent/project/D-NAscent-Origins/03b-frequencyClassification/awk-term/\
#final/nonSynchronized-moreFrequent.bed input/

# Making symbolic link to reference genome.
ln -s /project/carol/dnascent/project/D-NAscent-Origins/\
metaData/tryCru-clb7.gff input/genome.gff
ln -s /project/carol/dnascent/project/D-NAscent-Origins/\
metaData/tryCru-clb7.fasta input/genome.fasta


# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/D-NAscent-Origins/04-Annotation/bedtools/* .

# Running bedtoolsDriver script.
saveCommand script/bedtoolsDriver.bash\
  2>&1 | tee log/bedtoolsDriver.out

# Checking result.
saveCommand script/bedtoolsCheck.bash\
  2>&1 | tee log/bedtoolsCheck.out

# Removing intermediate files.
saveCommand script/bedtoolsClean.bash\
  2>&1 | tee log/bedtoolsClean.out
EOI

exit 0
