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

# Calculating the value of amplicons
for FILE in \`ls input/*.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" {print \$1, \$2, \$3 , (\$3 -\$2)}' \\
  input/\${FILE} \\
  > output/\${FILE}
done

# Sorting bed files by chromosomes and by starting position.
for FILE in \`ls output/*.bed | xargs -i basename {}\`
do
  sort -k1,1 -k2,2n -k3,3n output/\${FILE} \\
  > output/\${FILE%.bed}-sorted.bed
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 0  && \$4 <= 500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-1to500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 501  && \$4 <= 1000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-501to1000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 1001  && \$4 <= 1500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-1001to1500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 1501  && \$4 <= 2000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-1501to2000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 2001  && \$4 <= 2500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-2001to2500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 2501  && \$4 <= 3000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-2501to3000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 3001  && \$4 <= 3500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-3001to3500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 3501  && \$4 <= 4000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-3501to4000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 4001  && \$4 <= 4500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-4001to4500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 4501  && \$4 <= 5000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-4501to5000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 5001  && \$4 <= 5500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-5001to5500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 5501  && \$4 <= 6000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-5501to6000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 6001  && \$4 <= 6500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-6001to6500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 6501  && \$4 <= 7000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-6501to7000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 7001  && \$4 <= 7500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-7001to7500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 7501  && \$4 <= 8000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-7501to8000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 8001  && \$4 <= 8500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-8001to8500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 8501  && \$4 <= 9000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-8501to9000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 9001  && \$4 <= 9500 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-9001to9500.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 9501  && \$4 <= 1000 { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-9501to1000.txt
done

# select only the replicons with the different length
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do 
  awk 'OFS ="\t" && \$4 >= 10001  { print \$0 }'\\
  output/\${FILE} | wc -l \\
  > output/\${FILE%-sorted.bed}-from10001.txt
done

# concatenate all files containing amplicon length information to syncronized cells

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
