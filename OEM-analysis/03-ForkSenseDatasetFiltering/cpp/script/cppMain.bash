#!/usr/bin/env bash

###################################################################################################
#                                            cppMain                                              #
#                                                                                                 #
# In this script we'll use cpp to build databases of DNA replication origins for synchronized and #
# unsynchronized cells in the S phase of the cycle, as identified by the "D-NAscent" technique.   #
#                                                                                                 #
# Usage: saveCommand script/cppMain.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.franco@fundacaobutantan.org.br >                                   #
#                       David da Silva Pires                                                      #
#                      <david.pires@fundacaobutantan.org.br >                                     #
#    13/09/2024: First version.                                                                   #
###################################################################################################

# Writing the job script.
cat > job/cppMain.slurm << EOI
#!/usr/bin/env bash

# Creating chrom.sizes
faToTwoBit \${GENOME_FASTA} output/genome.2bit
twoBitInfo  output/genome.2bit stdout | sort -k2nr \
  > output/genome-chromSizes.txt

# Sort the genomeCromoSize. 
sort -k1,1 -k2,2n output/genome-chromSizes.txt > output/genomeChromoSizes.txt
## remove intermediated file. 
rm output/genome-chromSizes.txt
rm  output/genome.2bit

# Generate 20bp windows from the regions in the input BED file
bedtools makewindows -g output/genomeChromoSizes.txt -w \${WINDOW_SIZE} \
  > output/genomeSequencesWith-\${WINDOW_SIZE}bp.bed 

# Sort the Window genome  file.
sort -k1,1 -k2,2n -k3,3n -T output output/genomeSequencesWith-\${WINDOW_SIZE}bp.bed \
> output/genome-\${WINDOW_SIZE}bp.bed
rm output/genomeSequencesWith-20bp.bed

# Run program in C++
nice -n 19 /usr/bin/time --verbose --output=log/probability_filter.time \
  memusg --output-to-file log/probability_filter.memusg --time --shell "saveCommand script/probability_filter" \
  2> log/probability_filter.err | tee log/probability_filter.out

# Wait for the process to finish
wait

# Sort bedgraph files based on the chromosome order defined in chromSizes.txt
for FILE in \$(ls output/*.bedgraph | xargs -I {} basename {})
do
  sort -k1,1 -k2,2n -k3,3n -T output output/\${FILE} \
> output/\${FILE%.bedgraph}-sorted.bedgraph
done

# Create BEDMAP file using sorted files
for BEDGRAPH in \$( ls output/*-sorted.bedgraph | xargs -I {} basename {})
do
  bedtools map -a output/genome-\${WINDOW_SIZE}bp.bed -b output/\${BEDGRAPH} -c 4 -o mean \
  -g output/genomeChromoSizes.txt \
  > output/\${BEDGRAPH%.bedgraph}-aggregate.bedgraph
done 

# Process aggregate.bedgraph files
for AGGREGATE in \$( ls output/*-aggregate.bedgraph | xargs -I {} basename {})
do
  awk '{ if (\$4 == ".") \$4 = 0; print }' output/\${AGGREGATE} \
  > output/\$(AGGREGATE%-aggregate.bedgraph)-processed.bedgraph
done

exit 0
EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose \
--output=log/sbatch-cppMain.time memusg \
   --output-to-file log/sbatch-cppMain.memusg \
   --time --shell "saveCommand sbatch --nodes 1 \
   --ntasks \${NUM_OF_CPUS} \
   --mem \${MEMORY_SIZE} \
   -o log/slurm-%A.out \
   -J cppMain job/cppMain.slurm 2> \
   log/sbatch-cppMain.err | \
   tee log/sbatch-cppMain.out"

# Check how are going all the running awk processes.
echo "To check the running processes, execute one of the following commands:"
echo "   \$> watch -n 1 squeue"

exit 0
