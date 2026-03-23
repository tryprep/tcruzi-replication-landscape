#!/usr/bin/env bash

###################################################################################################
#                               cppPre                                                            #
#                                                                                                 #
# This script is a required pre-processing step to composing the D-NAscent origins  database.     #
# It sets up the directory structure and creates symbolic links for input data.                   #
#                                                                                                 #
# Usage: saveCommand script/cppPre.bash                                                           #
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

# Splitting the exported string into an array
IFS=' ' read -ra SAMPLE_ARRAY <<< "$SAMPLES"

# Creating input directories in parallel
for SAMPLE in "${SAMPLE_ARRAY[@]}"
do
  mkdir -p input/${SAMPLE}
  mkdir -p log/${SAMPLE}
  mkdir -p output/${SAMPLE} &
done
wait  # Aguarda a criação de todos os diretórios antes de prosseguir

# Making symbolic links for input data (forkSense)
for SAMPLE in "${SAMPLE_ARRAY[@]}"
do
  for FILE in /project/carol/dnascent/project/tryCru-clb1/03-brduCallingOrigins/DNAscent/sinc/${SAMPLE}/output/*.brduDetect.forkSense
  do
    ln -s ${FILE} input/${SAMPLE}/
  done
done

# Making symbolic links for input data (brduDetect)
for SAMPLE in "${SAMPLE_ARRAY[@]}"
do
  for FILE in /project/carol/dnascent/project/tryCru-clb1/03-brduCallingOrigins/DNAscent/nonSinc/${SAMPLE}/output/*.brduDetect
  do
    ln -s ${FILE} input/${SAMPLE}/
  done
done

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/D-NAscent-Origins/00-originsDataComposition/cpp/* .

# Running cpp driver script.
saveCommand script/cppDriver.bash 2>&1 | tee log/cppDriver.out

EOI
