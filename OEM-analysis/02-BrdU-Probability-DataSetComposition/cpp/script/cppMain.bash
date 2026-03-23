#!/usr/bin/env bash

###################################################################################################
#                                            cppMain                                              #
#                                                                                                 #
#    This script In this script, we'll use awk to build databases of DNA replication origins for  #
# synchronized and unsynchronized cells in the S phase of the cycle, as identified by the         #
# "D-Nascent" technique.                                                                          #
#                                                                                                 #
# Usage: saveCommand script/awkMain.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.franco@fundacaobutantan.org.br >                                   #
#                       David da Silva Pires                                                      #
#                      <david.pires@fundacaobutantan.org.br >                                     #
#    23/09/2024: First version.                                                                   #
###################################################################################################

# Writing the job script.
cat > job/cppMain.slurm << EOI
#!/usr/bin/env bash

# Run program in C++.
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/forkSenseProcessor.time \
  memusg --output-to-file log/forkSenseProcessor.memusg --time --shell "script/forkSenseProcessor" \
  2> log/forkSenseProcessor.err | tee log/forkSenseProcessor.out

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
echo " $> watch -n 1 squeue"

exit 0
