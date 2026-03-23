#!/usr/bin/env bash

###################################################################################################
#                                            awkMain                                              #
#                                                                                                 #
#    This script In this script, we'll use awk to build databases of DNA replication origins for  #
# synchronized and unsynchronized cells in the S phase of the cycle, as identified by the         #
# "D-Nascent" technique.                                                                          #
#                                                                                                 #
# Usage: saveCommand script/awkMain.bash                                                        #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #
#    08/03/2023: First version.                                                                   #
###################################################################################################

# Writing the job script.
cat > job/awkMain.slurm << EOI
#!/usr/bin/env bash

# Concatenate the files of the same treatment in a single file
for FILE in \`ls input/tableOrigins-NonSinc*  | xargs -i basename {}\`
do
   cat input/\${FILE} >> output/tableOriginsNonSinc.csv
done

# Concatenate the files of the same treatment in a single file
for FILE in \`ls input/tableOrigins-Sinc* | xargs -i basename {}\`
do
   cat input/\${FILE} >> output/tableOriginsSinc.csv
done

# selecting only DNA replication origns from construction table files
for FILE in \`ls output/*.csv | xargs -i basename {}\`
do
   grep -P '.+\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\tORI\t.*\t.*\t.*' output/\${FILE} \\
   > output/\${FILE%.csv}.tsv
done

# Concatenate the files of the same treatment in a single file
for FILE in \`ls output/*.tsv | xargs -i basename {}\`
do
  awk 'OFS = "\t" {print \$2, \$12, \$13, \$1, 1, \$6 }' output/\${FILE} \\
    > output/\${FILE%.tsv}.bed3+3
  done

# Composing bed3 file to D-NAscent origins
for FILE in \`ls output/*bed3+3 | xargs -i basename {}\`
do
  awk 'OFS = "\t" {print \$1, \$2, \$3 }' output/\${FILE} \\
> output/\${FILE%.bed3+3}.bed
done

# Sorting the bed files to next process. 
for BED_FILES in \`ls output/*.bed | xargs -i basename {}\`
do
  sort -k1,1 -k2,2n -k3,3n output/\${BED_FILES} \\
  > output/\${BED_FILES%.bed}-sorted.bed
done

# Creating genome index.
BASENAME_GENOME=`basename \${GENOME%.fa}`
samtools faidx \${GENOME} \
  -o output/\${BASENAME_GENOME}.fa.fai \
  2> log/index.err > log/index.out

# Creating chrom.sizes.
# https://genomewiki.ucsc.edu/index.php/GBiB:_From_download_to_BLAT_at_assemply_hubs
faToTwoBit \${GENOME} output/\${BASENAME_GENOME}.2bit
twoBitInfo output/\${BASENAME_GENOME}.2bit stdout | sort -k2nr \
  > output/\${BASENAME_GENOME}-chromSizes.txt

# Checking whether the genome coordinates are within the expected size of the reference genome.
for BED_FILE in output/*-sorted.bed; do
    BASENAME=\$(basename "\$BED_FILE")

    # Read the BED file line by line
    while IFS=\$'\t' read -r CHROM START END _; do
        CHROM_SIZE=\$(awk -v chrom="\$CHROM" '\$1 == chrom {print \$2}' "output/\${BASENAME_GENOME}-chromSizes.txt")

        # Check if the start coordinates are greater than or equal to zero.
        if [ "\$START" -lt 0 ]; then
            START=0
        fi

        # Check if the end coordinates are less than or equal to chromosome size.
        if [ "\$END" -gt "\$CHROM_SIZE" ]; then
            END="\$CHROM_SIZE"
        fi

        # Output the corrected coordinates to a new file
        echo -e "\${CHROM}\t\${START}\t\${END}"
    done < "\$BED_FILE" >> "output/\${BASENAME%-sorted.bed}-checked.bed"
done

exit 0
EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/sbatch-awkMain.time memusg \
   --output-to-file log/sbatch-awkMain.memusg --time --shell "saveCommand sbatch --nodes 1 \
   --ntasks ${NUM_OF_CPUS} --mem ${MEMORY_SIZE} -o log/slurm-%A.out -J awkMain job/awkMain.slurm 2> \
   log/sbatch-awkMain.err | tee log/sbatch-awkMain.out"

# Check how are going all the running awk processes.
echo "To check the running processes, execute one of the following commands:"
echo "   $> watch -n 1 squeue"

exit 0
