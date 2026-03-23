#!/usr/bin/env bash

###################################################################################################
#                                            awkMain                                              #
#                                                                                                 #
#    This script In this script, we'll use awk to build databases of DNA replication origins for  #
# synchronized and unsynchronized cells in the S phase of the cycle, as identified by the         #
# "D-Nascent" technique.                                                                          #
#                                                                                                 #
# Usage: saveCommand script/awkMain.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #
#    08/03/2023: First version.                                                                   #
###################################################################################################

# Writing the job script.
cat > job/awkMain.slurm << EOI
#!/usr/bin/env bash

# Selecionando as origens de replicação menos frequentes de células sincronizadas e não sincronizadas.
for FILE in \$(ls input/*.csv | xargs -i basename {})
do
  awk -F "\t" '\$4 <= 10 {print \$0}' \
  input/\${FILE} | sort -k4,4n \
  > output/\${FILE%.csv}-lessFrequent.tsv
done

# Selecionando as origens de replicação menos frequentes de células sincronizadas e não sincronizadas.
for FILE in \$(ls input/*.csv | xargs -i basename {})
do
  awk -F "\t" '\$4 > 10 {print \$0}' \
  input/\${FILE} | sort -k4,4n \
  > output/\${FILE%.csv}-moreFrequent.tsv
done

# Definindo as origens de replicação do DNA com uma janela de 1kb.
for FILE in \$(ls output/*.tsv | xargs -i basename {})
do
  awk 'BEGIN {OFS = "\t"} {print \$1, int(((\$3 - \$2) / 2 ) + \$2) - 1000, \
  int(((\$3 - \$2) / 2 ) +\$2) + 1000 }' output/\${FILE} \
  > output/\${FILE%.tsv}.bed
done

# Ordenando todos os arquivos bed por cromossomo, início e fim.
for BED_FILE in \$(ls output/*.bed | xargs -i basename {})
do
  sort -k1,1 -k2,2n -k3,3n output/\${BED_FILE} \
  > output/\${BED_FILE%.bed}-sorted.bed
done

# Criando índice do genoma.
BASENAME_GENOME=\$(basename \${GENOME%.fa})
samtools faidx \${GENOME} \
  -o output/\${BASENAME_GENOME}.fa.fai \
  2> log/index.err > log/index.out

# Criando chrom.sizes.
faToTwoBit \${GENOME} output/\${BASENAME_GENOME}.2bit
twoBitInfo output/$\{BASENAME_GENOME}.2bit stdout | sort -k2nr \
  > output/\${BASENAME_GENOME}-chromSizes.txt

# Verificando se as coordenadas do genoma estão dentro do tamanho esperado do genoma de referência.
for BED_FILE in output/*-sorted.bed; do
    BASENAME=\$(basename "\$BED_FILE")

    # Lendo o arquivo BED linha por linha
    while IFS=$'\t' read -r CHROM START END _; do
        CHROM_SIZE=\$(awk -v chrom="\$CHROM" '\$1 == chrom {print \$2}' "output/\${BASENAME_GENOME}-chromSizes.txt")

        # Verificar se as coordenadas de início são maiores ou iguais a zero.
        if [ "\$START" -lt 0 ]; then
            START=0
        fi

        # Verificar se as coordenadas de fim são menores ou iguais ao tamanho do cromossomo.
        if [ "\$END" -gt "\$CHROM_SIZE" ]; then
            END="\$CHROM_SIZE"
        fi

        # Saída das coordenadas corrigidas para um novo arquivo
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
