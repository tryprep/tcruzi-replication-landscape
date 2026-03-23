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
cat> job/bedtoolsMain.slurm <<EOI
#!/usr/bin/env bash

# Initially we will create the database for the annotation of the compartments
# Selecting only lines with "protein coding gene" at third column.
# In this step, we will create gene.gff to all analysis.
  grep -P '.+\t.*\tprotein_coding_gene\t.*\t.*\t.*\t.*\t.*\t.*' \\
  input/genome.gff \\
   > output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-protein_coding_gene.gff

# In this step, we will create the files of each protein of the multigene family
for FILE in 'DGF-1' 'GP63' 'RHS' 'trans-sialidase' 'MASP'
do
  /usr/bin/nice -n 19 grep -P '.+\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'\${FILE}'.*' \\
  output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-protein_coding_gene.gff \\
  > output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-\${FILE}.gff &
done

wait

# And then,  we will select the mucin genes
  grep mucin  output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-protein_coding_gene.gff|
  grep -v mucin-associated \\
  > output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-mucin.gff

# Making Core compartment
    grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'DGF-1$'.*' \\
    output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-protein_coding_gene.gff|
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'GP63$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'MASP$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'RHS$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'mucin$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'trans-sialidase$'.*' \\
  > output/core-Compartment.gff

# Making Disruptive compartment
cat output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-mucin.gff \\
  output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-MASP.gff \\
  output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-trans-sialidase.gff \\
  > output/disruptive-Compartment.gff

# Making Both compartment
cat  output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-DGF-1.gff \\
  output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-RHS.gff \\
  output/TriTrypDB-60_TcruziCLBrenerEsmeraldo-like-gene-GP63.gff \\
  > output/both-Compartment.gff

exit 0

EOI

# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-bedtoolsMain.time memusg \
  --output-to-file log/sbatch-bedtoolsMain.memusg \
    --time --shell "saveCommand sbatch --parsable \
    --nodes 1 \
  --ntasks ${NUM_OF_CPUS} \
  --mem ${MEMORY_SIZE} \
  -o log/slurm-%A.out \
  -J bedtoolsMain job/bedtoolsMain.slurm \
  2> log/sbatch-bedtoolsMain.err |
    tee log/sbatch-bedtoolsMain.out"

# Check how are going all the running awk processes.
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0
