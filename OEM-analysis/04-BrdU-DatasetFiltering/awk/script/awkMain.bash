#!/usr/bin/env bash

###################################################################################################
#                                            awkMain                                              #
#                                                                                                 #
# In this script we'll use awk to build databases of DNA replication origins for synchronized and #
# unsynchronized cells in the S phase of the cycle, as identified by the "D-NAscent" technique.   #
#                                                                                                 #
# Usage: saveCommand script/awkMain.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.franco@fundacaobutantan.org.br>                                    #
#                      David da Silva Pires                                                       #
#                      <david.pires@fundacaobutantan.org.br>                                      #
#    23/09/2024: First version.                                                                   #
###################################################################################################

# Writing the job script.
cat > job/awkMain.slurm << EOI
#!/usr/bin/env bash

#!/bin/bash

# Create genome index
BASENAME_GENOME=\$(basename "\${GENOME_FASTA%.fa}")
samtools faidx "\${GENOME_FASTA}" \
  2> log/index.err \
  > log/index.out

# Create chrom.sizes file
faToTwoBit "\${GENOME_FASTA}" output/"\${BASENAME_GENOME}.2bit"
twoBitInfo output/"\${BASENAME_GENOME}.2bit" stdout | sort -k2nr \
  > output/"\${BASENAME_GENOME}-chromSizes.txt"

# Generate 20bp windows from the regions in the input BED file
bedtools makewindows -g output/"\${BASENAME_GENOME}-chromSizes.txt" -w "\${WINDOW_SIZE}" \
  > output/genomeSequencesWith-"\${WINDOW_SIZE}"bp.bed

# Keep values and change values less than 0.5 to 0
awk -v prob="\${PROBABILITY}" '{
  if (\$5 < prob) \$5 = 0;
  print
}' input/probability-BrdU.csv > output/probability-BrdU-KeepValues.csv

# Change values >= 0.5 to 1, < 0.5 to 0
awk -v prob="\${PROBABILITY}" '{
  if (\$5 >= prob) \$5 = 1;
  else \$5 = 0;
  print
}' input/probability-BrdU.csv > output/probability-BrdU-ChangesValues.csv

# Create bedgraph files for bedtools mapBed processes
for FILE in output/*.csv; do
  awk 'OFS="\t" {print \$2, \$3, \$4, \$5}' "\${FILE}" |\
  sort -k1,1 -k2,2n -k3,3n -k4,4n \
  > output/"\$(basename "\${FILE%.csv}").bedgraph"
done

# Sort bedgraph files based on the chromosome order defined in genome-chromSizes.txt
for FILES in output/*.bedgraph; do
  bedtools sort -i "\${FILES}" -g output/"\${BASENAME_GENOME}-chromSizes.txt" > output/"\$(basename "\${FILES%.bedgraph}")_sorted.bedgraph"
done

# Create BEDMAP file using sorted files
for FILES in output/*_sorted.bedgraph; do
  bedtools map -a output/genomeSequencesWith-"\${WINDOW_SIZE}"bp.bed -b "\${FILES}" -c 4 -o mean \
  -g output/"\${BASENAME_GENOME}-chromSizes.txt" > output/"\$(basename "\${FILES%_sorted.bedgraph}")-aggregate.bedgraph"
done

EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose \
--output=log/sbatch-awkMain.time memusg \
   --output-to-file log/sbatch-awkMain.memusg \
   --time --shell "saveCommand sbatch --nodes 1 \
   --ntasks \${NUM_OF_CPUS} \
   --mem \${MEMORY_SIZE} \
   -o log/slurm-%A.out \
   -J awkMain job/awkMain.slurm 2> \
   log/sbatch-awkMain.err | \
   tee log/sbatch-awkMain.out"

# Check how are going all the running awk processes.
echo "To check the running processes, execute one of the following commands:"
echo "   \$> watch -n 1 squeue"

exit 0
