#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #
#                                                                                                 #
#                                                                                                 #
#    This script allows you to record the origins of DNA replication in the genome. We generated  #
# and analyzed the sections of the genome where the origins of DNA replication are located using  #
# the "intersect" tool of the "bedtools" program. In summary, bedtools will compare the genomic   #
# coordinates of the DNA replication origins to the genomic coordinates of the reference genome   #
# and output a file containing the intersection regions. As a result, we will examine each DNA    #
# replication origins file to determine where genomic areas intersected. It writes the Slurm job  #
# script and submits all of them.                                                                 #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023  Thiago Andrade Franco                                                        #
#                    thiago.franco.esib@esib.butantan.gov.br                                      #
#    08/07/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/bedtoolsMain2.slurm <<EOI
#!/usr/bin/env bash

# Sorting bed files by chromosomes and by starting position.
for i in \`ls output/*.tsv | xargs -i basename {}\`
do
 sort -k6,6  output/\${i} \\
> output/\${i%.tsv}-sorted.tsv
done


# Selecting only lines with "protein coding gene" at the third column. In this step, we will create
# gene File to analysis.
for FILE in \`ls output/*-sorted.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t.*'\\
  output/\${FILE} \\
  > output/\${FILE%-sorted.tsv}-pseudogene.tsv
done

# Select only lines with DGF-1 
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t*DGF-1*'\\
  output/\${FILE} \\
  > output/\${FILE%Overlappinggenome-pseudogene.tsv}-DGF1.tsv
done

# Select only lines with RHS 
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t*RHS*'\\
  output/\${FILE} \\
  > output/\${FILE%Overlappinggenome-pseudogene.tsv}-RHS.tsv
done

# Select only lines with MASP 
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t*MASP*'\\
  output/\${FILE} \\
  > output/\${FILE%Overlappinggenome-pseudogene.tsv}-MASP.tsv
done

# Select only lines with Mucin 
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t*mucin$*'\\
  output/\${FILE} |
  grep -v mucin-associated \\
  > output/\${FILE%Overlappinggenome-pseudogene.tsv}-MUC.tsv
done

# Select only lines with trans-sialidase 
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t*trans-sialidase*'\\
  output/\${FILE} \\
  > output/\${FILE%Overlappinggenome-pseudogene.tsv}-TS.tsv
done

# Select only lines with GP63
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t*GP63*'\\
  output/\${FILE} \\
  > output/\${FILE%Overlappinggenome-pseudogene.tsv}-GP63.tsv
done

# Select only lines with 
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -P '.+\t.*\t.*\t.*\t.*\tpseudogene\t.*\t.*\t.*\t.*\t.*\t*GP63*'\\
  output/\${FILE} \\
  > output/\${FILE%Overlappinggenome-pseudogene}-GP63.tsv
done

# Making others genes
for FILE in \`ls output/*-pseudogene.tsv | xargs -i basename {}\`
do
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'DGF-1$'.*' \\
  output/\${FILE} |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'GP63$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'MASP$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'RHS$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'mucin$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'trans-sialidase$'.*' \\
  > output/\${FILE%Overlappinggenome-pseudogene.tsv}-OTHERS.tsv
done

FEATURES=("DGF1" "GP63" "MASP" "MUC" "TS" "RHS" "OTHERS")

# Creating the file with features enriched for less frequent DNA replication termination
for FEATURE in \${FEATURES[@]}
do
    for NONSYNC in output/nonSynchronized-lessFrequent-\${FEATURE}.tsv
    do
        if [[ -f \$NONSYNC ]]; then
            COUNT=\$(wc -l < "\$NONSYNC")
            echo -e "\${FEATURE}\t\${COUNT}\tLess frequent origins" \
            >> "output/NonSyncPseudogene-lessFrequent.csv"
        fi
    done
done

# Creating the file with features enriched for more frequent DNA replication termination
for FEATURE in \${FEATURES[@]}
do
    for NONSYNC in output/nonSynchronized-moreFrequent-\${FEATURE}.tsv
    do
        if [[ -f \$NONSYNC ]]; then
            COUNT=\$(wc -l < "\$NONSYNC")
            echo -e "\${FEATURE}\t\${COUNT}\tMore frequent origins" \
            >> "output/NonSyncPseudogene-MoreFrequent.csv"
        fi
    done
done

# Creating the file with features enriched for general DNA replication termination (D-NAscent)
for FEATURE in \${FEATURES[@]}
do
    for NONSYNC in output/nonSynchronized-\${FEATURE}.tsv
    do
        if [[ -f \$NONSYNC ]]; then
            COUNT=\$(wc -l < "\$NONSYNC")
            echo -e "\${FEATURE}\t\${COUNT}\tD-NAscent origins" \
            >> "output/NonSyncPseudogene.csv"
        fi
    done
done

exit 0
EOI

MAIN_JOB_ID=`cat log/sbatch-bedtoolsMain.out`

# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-bedtoolsMain2.time memusg \
   --output-to-file log/sbatch-bedtoolsMain2.memusg \
   --time --shell "saveCommand sbatch \
   --dependency=afterany:${MAIN_JOB_ID} \
   --nodes 1 \
   --ntasks ${NUM_OF_CPUS} --mem ${MEMORY_SIZE} \
   -o log/slurm-%A.out \
   -J bedtoolsMain2 job/bedtoolsMain2.slurm \
   2> log/sbatch-bedtoolsMain2.err |
     tee log/sbatch-bedtoolsMain2.out"

# Check how are going all the running awk processes.
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0
