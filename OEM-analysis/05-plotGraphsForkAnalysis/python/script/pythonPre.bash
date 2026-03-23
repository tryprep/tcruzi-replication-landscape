#!/usr/bin/env bash

###################################################################################################
#                               pythonPre                                                         #
#                                                                                                 #
# This script is a required pre-processing step to composing the D-NAscent origins  database.     #
# It sets up the directory structure and creates symbolic links for input data.                   #
#                                                                                                 #
# Usage: saveCommand script/pythonPre.bash                                                        #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.franco@fundacaobutantan.org.br >                                   #
#                       David da Silva Pires                                                      #
#                      <david.pires@fundacaobutantan.org.br >                                     #
#    13/09/2024: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output

# Storing initial processing date and time.
date > log/processingStart.date

# Making symbolic link to genome files
ln -s /project/carol/dnascent/project/D-NAscent-Origins/metaData/tryCru-clb7.gff input/genome.gff
ln -s /project/carol/dnascent/project/D-NAscent-Origins/metaData/tryCru-clb7.fasta input/genome.fa

# Making symbolic link to genome file to process. 
ln -s /project/carol/dnascent/project/IdentificationOfPotentialMachiberyConflictRegions/\
03-ForkSenseDatasetFiltering/awk/output/probabilityForkSense-KeepValues-left-aggregate.bedgraph \
input/whatson.bedgraph
ln -s /project/carol/dnascent/project/IdentificationOfPotentialMachiberyConflictRegions/\
03-ForkSenseDatasetFiltering/awk/output/probabilityForkSense-KeepValues-Rigth-aggregate.bedgraph \
input/click.bedgraph




# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/D-NAscent-Origins/00-originsDataComposition/python/* .

# Running python driver script.
saveCommand script/pythonDriver.bash 2>&1 | tee log/pythonDriver.out

EOI
