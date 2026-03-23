#!/usr/bin/env bash

###################################################################################################
#                                         awkPost                                               #
#                                                                                                 #
#   This script post-process data produced by awk.It writes a job script that making symbolic   #
# links for the final result of the positive fork > 0.7 BrdU.                                      #
#                                                                                                 #
# Usage: saveCommand script/awkPost.bash                                                        #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Marcela de Oliveira Vitarelli and Thiago Andrade Franco                    #
#                    <vitarelli.marcela@gmail.com> and thiago.franco.esib@esib.butantan.gov.br    #
#   28/07/2023: First version.                                                                   #
###################################################################################################


# Creating a Bash script to post-process Scipio result.
cat > job/awkPost.slurm << EOI
#!/usr/bin/env bash

# making symbolic link to final of the DNA replication origins database 
for i in \`ls output/*.bed | xargs -i basename {}\`
do
      ln -s ../output/\${i} final/\${i%.bed}-3kb.bed 
done


EOI

exit 0
# Submitting the job to Slurm.
      /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/sbatch-awkPost.time memusg \
        --output-to-file log/sbatch-awkPost.memusg --time --shell "saveCommand sbatch --nodes 1 \
        --ntasks ${NUM_OF_CPUS} --mem ${MEMORY_SIZE} -o log/slurm-%A.out -J awkPost job/awkPost.slurm \
        2> log/sbatch-awkPost.err | tee log/sbatch-awkPost.out"

# Check how are going all the running awk processes.
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0
