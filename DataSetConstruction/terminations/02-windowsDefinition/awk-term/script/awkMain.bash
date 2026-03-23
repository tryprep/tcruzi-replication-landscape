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

# defining the DNA replication origins with 1kb of window
for FILE in \`ls input/*.bed | xargs -i basename {}\`
do
  awk 'OFS = "\t" {print \$1, int(((\$3 - \$2) / 2 ) + \$2) - 1000, \\
  int(((\$3 - \$2) / 2 ) +\$2) + 1000 }' input/\${FILE} \\
  > output/\${FILE%.bed}-2kb.bed3
done

# Sorting bed files by chromosomes and by starting position.
for i in \`ls output/*-2kb.bed3 | xargs -i basename {}\`
do
   sort -k1,1 -k2,2n -k3,3n output/\${i} \\
   > output/\${i%-2kb.bed3}-sorted.bed
done

# We cluster all coordinates that have at least 3 kB of intersection.
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
  bedtools merge -i output/\${FILE} -d 1000 -c 1 -o count > \\
    output/\${FILE%-sorted.bed}-merged1kb.csv \\
    2> log/\${FILE%-sorted.bed}-merged1kb.err
done

# Sorting merged1kb.bed files by chromosomes and by starting position.
for i in \`ls output/*.csv | xargs -i basename {}\`
do
   sort -k1,1 -k2,2n -k3,3n output/\${i} \\
   > output/\${i%-merged1kb.csv}-sorted.tsv
done

# defining the DNA replication origins with 1kb of window
for FILE in \`ls output/*-sorted.tsv | xargs -i basename {}\`
do
  awk 'OFS = "\t" {print \$1, int(((\$3 - \$2) / 2 ) + \$2) - 1000, \\
  int(((\$3 - \$2) / 2 ) +\$2) + 1000 }' output/\${FILE} \\
  > output/\${FILE%-sorted.tsv}-1kb.bed
done

# renaming the origins files
for FILE in \`ls output/*-1kb.bed | xargs -i basename {}\`
do
  echo output/\${FILE} | sed -s 's/Amplicon/Origins/g'
done


# defining the center point of genomic coordinates of the DNA replication origins
for FILE in \`ls output/*-1kb.bed | xargs -i basename {}\`
do
   awk 'function abs(v) { return v < 0 ? -v : v } function min(a, b) { return a > b ? b : a } \
     OFS = "\t" {print \$1, int(abs(\$3 - \$2) / 2) + min(\$3, \$2), int(abs(\$3 - \$2) / 2) + \
     min(\$3, \$2) }' output/\${FILE} \\
        > output/\${FILE%-1kb.bed}-CenterPoint.bed
done

# Prepare input data to calculate the frequency of D-nascent DNA replication origins.
for FILE in \`ls output/*-merged1kb.csv | xargs -i basename {}\`
do 
  cut -f4 output/\${FILE} > output/\${FILE%-merged1kb.csv}.txt
done

# prepare input data to plot length graph to clusterized amplicon
for FILE in \`ls output/*-merged1kb.csv | xargs -i basename {}\`
do 
  awk 'OFS="\t" {print \$3 -\$2}' output/\${FILE} \\
  > output/\${FILE%.csv}-Length.txt
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
