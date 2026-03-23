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
for i in \`ls input/*.bed | xargs -i basename {}\`
do
  sort -k1,1 -k2,2n -k3,3n input/\${i} \\
  > output/\${i%.bed}-sorted.bed
done

# Sorting gff files by chromosomes and by starting position.
for i in \`ls output/*Compartment.gff | xargs -i basename {}\`
do
  sort -k1,1  -k4,4n -k5,5n output/\${i} \\
  > output/\${i%.gff}-sorted.gff
done

# Intersect the bed files from DNA replication Origins with core compartment
for FILE in \`ls output/*.bed | xargs -i basename {}\`
do
  bedtools intersect -wo -a output/\${FILE} \\
    -b output/core-Compartment-sorted.gff \\
    > output/\${FILE%-sorted.bed}OverlappingCore.tsv \\
    2> log/\${FILE%-sorted.bed}OverlappingCore.err \\

done

# Intersect the bed files from DNA replication Origins with disruptive compartment
for FILE in \`ls output/*.bed | xargs -i basename {}\`
do
  bedtools intersect -wo -a output/\${FILE} \\
    -b output/disruptive-Compartment-sorted.gff \\
    > output/\${FILE%-sorted.bed}OverlappingDisruptive.tsv \\
    2> log/\${FILE%-sorted.bed}OverlappingDisruptive.err \\

done

# Intersect the bed files from DNA replication Origins with both compartment
for FILE in \`ls output/*.bed | xargs -i basename {}\`
do
  bedtools intersect -wo -a output/\${FILE} \\
    -b output/both-Compartment-sorted.gff \\
    > output/\${FILE%-sorted.bed}OverlappingBoth.tsv \\
    2> log/\${FILE%-sorted.bed}OverlappingBoth.err \\

done

#######################################################################################################
FEATURES=("Both" "Core" "Disruptive")

# Creating the file with features enriched for less frequent DNA replication termination
for FEATURE in \${FEATURES[@]}
do
    for NONSYNC in output/nonSynchronized-lessFrequentOverlapping\${FEATURE}.tsv
    do
        if [[ -f \$NONSYNC ]]; then
            COUNT=\$(wc -l < "\$NONSYNC")
            echo -e "\${FEATURE}\t\${COUNT}\tLess frequent origins" \
            >> "output/NonSyncGenomeCompartment-lessFrequent.csv"
        fi
    done
done

# Creating the file with features enriched for more frequent DNA replication termination
for FEATURE in \${FEATURES[@]}
do
    for NONSYNC in output/nonSynchronized-moreFrequentOverlapping\${FEATURE}.tsv
    do
        if [[ -f \$NONSYNC ]]; then
            COUNT=\$(wc -l < "\$NONSYNC")
            echo -e "\${FEATURE}\t\${COUNT}\tMore frequent origins" \
            >> "output/NonSyncGenomeCompartment-MoreFrequent.csv"
        fi
    done
done

# Creating the file with features enriched for general DNA replication termination (D-NAscent)
for FEATURE in \${FEATURES[@]}
do
    for NONSYNC in output/nonSynchronizedOverlapping\${FEATURE}.tsv
    do
        if [[ -f \$NONSYNC ]]; then
            COUNT=\$(wc -l < "\$NONSYNC")
            echo -e "\${FEATURE}\t\${COUNT}\tD-NAscent origins" \
            >> "output/NonSyncGenomeCompartment.csv"
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
