#!/usr/bin/env bash

###################################################################################################
#                                            pythonMain                                           #
#                                                                                                 #
# In this script we'll use awk to build databases of DNA replication origins for synchronized and #
# unsynchronized cells in the S phase of the cycle, as identified by the "D-NAscent" technique.   #
#                                                                                                 #
# Usage: saveCommand script/pythonMain.bash                                                       #
#                                                                                                 #
# Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                      #
#                      <thiago.franco@fundacaobutantan.org.br >                                   #
#                       David da Silva Pires                                                      #
#                      <david.pires@fundacaobutantan.org.br >                                     #
#    13/09/2024: First version.                                                                   #
###################################################################################################

# Writing the job script.
cat > job/pythonMain.slurm << EOI
#!/usr/bin/env bash

# Change . by 0 in BEDGRAPH files.

for FILE in \$( ls input/*.bedgraph | xargs -I {} basename {})
do
  awk '{ if (\$4 == ".") \$4 = 0; print }' input/\${FILE} \
  > output/\${FILE}
done

# Run program in Python to OEM
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/OEM-ParseWindowAndStepSizes-NormWholeGen.time \
  memusg --output-to-file log/OEM-ParseWindowAndStepSizes-NormWholeGen.memusg --time \
  python script/OEM-ParseWindowAndStepSizes-NormWholeGen.py \
  output/whatson.bedgraph output/click.bedgraph output/oem.bedgraph -w \${WINDOW_SIZE} -s \${STEP_SIZE} \
  2> log/OEM-ParseWindowAndStepSizes-NormWholeGen.err | tee log/OEM-ParseWindowAndStepSizes-NormWholeGen.out

# Run program in Python to OEM
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/RFD_Calculation.time \
  memusg --output-to-file log/RFD_Calculation.memusg --time \
  python script/RFD_Calculation-A.py \
  output/whatson.bedgraph output/click.bedgraph output/rfd.bedgraph -w \${WINDOW_SIZE} -s \${STEP_SIZE} \
  2> log/RFD_Calculation.err | tee log/RFD_Calculation.out

exit 0

EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose \
--output=log/sbatch-pythonMain.time memusg \
   --output-to-file log/sbatch-pythonMain.memusg \
   --time --shell "saveCommand sbatch --nodes 1 \
   --ntasks \${NUM_OF_CPUS} \
   --mem \${MEMORY_SIZE} \
   -o log/slurm-%A.out \
   -J pythonMain job/pythonMain.slurm 2> \
   log/sbatch-pythonMain.err | \
   tee log/sbatch-pythonMain.out"

# Check how are going all the running awk processes.
echo "To check the running processes, execute one of the following commands:"
echo "   \$> watch -n 1 squeue"

exit 0
