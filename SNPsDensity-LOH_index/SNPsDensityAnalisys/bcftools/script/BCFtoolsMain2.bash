#!/usr/bin/env bash

###################################################################################################
#                                         BCFtoolsMain2                                            #
#                                                                                                 #
#  Description:                                                                                   #
#  This script executes the primary processing step of the SNP analysis pipeline. It utilizes     #
#  bedtools (intersect and subtract), awk, grep, and sort to process and generate BED files       #
#  defining specific genomic regions of interest.                                                 #
#                                                                                                 #
#  These targeted coordinate regions include DNA replication origins, replication termini,        #
#  genomic compartments, polycistronic transcription units (polycistrons), and areas that         #
#  intersect (or do not intersect) with multigene families.                                       #
#                                                                                                 #
#  The resulting BED files are essential for delimiting strict genomic boundaries for all         #
#  downstream pipeline analyses, enabling accurate assessment of SNP density, distribution,       #
#  and evolutionary metrics across these distinct genomic features.                               #
#                                                                                                 #
#  Usage:                                                                                         #
#  saveCommand script/BCFtoolsMain2.bash                                                           #
#                                                                                                 #
#                                                                                                 #
#  Copyleft (ɔ) 2025 by Thiago Andrade Franco                                                     #
#                       <thiago.franco@fundacaobutantan.org.br>                                   #
#                                                                                                 #
#  23/04/2025: First version.                                                                     #
###################################################################################################
# Define multigenic families
GENES="DGF-1 GP63 RHS MUC TS MASP HIST"

# Create SLURM job script
cat > job/BCFtoolsMain2.slurm <<'EOI'
#!/usr/bin/env bash

# === Module 1: Filter high-quality variants ===
for FILE in input/*.vcf
do
  FILE_BASENAME=$(basename "$FILE")
  bcftools filter -i 'QUAL>20 && INFO/DP>20 && INFO/MQM>50' "$FILE" \
  > "output/${FILE_BASENAME}"
done

# === Module 2: Retain only SNPs ===
for VCF in output/*.vcf
do
  BASENAME=$(basename "$VCF")
  bcftools view -v snps "$VCF" |\
  grep -v -E "TYPE=(del|ins|complex|mnp)" \
  > "output/SNPs-${BASENAME}"
done

# === Module 3: Variant annotation using SnpEff ===
for SNPS in output/SNPs-*.vcf
do
  BASENAME=$(basename "$SNPS")
  snpeff -v trypanosoma_cruzi_clbrener "$SNPS" \
  > "output/${BASENAME%.vcf}-annotated.vcf"
done

# === Module 4: Run C++ filter ===
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/filter_vcf.time \
  memusg --output-to-file=log/filter_vcf.memusg \
  script/filter_vcf \
  > log/filter_vcf.out \
  2> log/filter_vcf.err

# === Module 5: Index filtered SNP VCFs ===
for FILE in output/SNPs-sample*-filtered.vcf
do
  bgzip -c "$FILE" > "${FILE}.gz"
  bcftools index "${FILE}.gz"
done

# === Module 6: General SNP density table ===
for GROUP in output/*/
do
  GROUP_NAME=$(basename "$GROUP")
  OUTPUT_FILE="${GROUP}/SNPs-${GROUP_NAME}-Summary.csv"
  echo -e "Sample\tSNP_Count\tLength\tGroup\tSNPs density" \
  > "$OUTPUT_FILE"

  for BED in "${GROUP}"/*.bed
  do
    [ -s "$BED" ] || continue
    REGION=$(basename "$BED" .bed)
    LENGTH=$(awk '{sum += $3 - $2} END {print sum}' "$BED")

    for VCF in output/SNPs-sample*-filtered.vcf.gz
    do
      SAMPLE=$(basename "$VCF" .vcf.gz | sed 's/SNPs-//' | sed 's/-filtered//')
      SNP_COUNT=$(bcftools view -R "$BED" "$VCF" | grep -vc '^#')
      DENSITY=$(awk -v n="$SNP_COUNT" -v l="$LENGTH" 'BEGIN{printf "%.10f", n/l}')
      echo -e "$SAMPLE\t$SNP_COUNT\t$LENGTH\t$REGION\t$DENSITY" >> "$OUTPUT_FILE"
    done
  done
done

# === Module 7: Summary for multigenic families ===
for FAMILY in ${GENES}
do
  SUMMARY_FILE="output/${FAMILY}/SNPs-${FAMILY}-Summary.csv"
  [ -s "$SUMMARY_FILE" ] && grep -P "\t${FAMILY}\t" "$SUMMARY_FILE" \
  >> "output/SNPsMultigenicSummary.csv"
done

# === Module 8: SNP density per gene ===
for FAMILY in ${GENES}
do
  FAMILY_BED="output/${FAMILY}/${FAMILY}.bed"
  [ -s "$FAMILY_BED" ] || continue

  OUTPUT_FILE="output/${FAMILY}/SNPs-${FAMILY}-genebygene-summary.csv"
  echo -e "Sample\tSNP_Count\tLength\tGene ID\tStrand\tSNPs density" \
  > "$OUTPUT_FILE"

  while read -r CHR START END ID _ STRAND
  do
    BED_TEMP="output/${FAMILY}/${ID}.bed"
    echo -e "$CHR\t$START\t$END\t$ID\t$STRAND" > "$BED_TEMP"
    LENGTH=$((END - START))

    for VCF in output/SNPs-sample*-filtered.vcf.gz
    do
      SAMPLE=$(basename "$VCF" .vcf.gz | sed 's/SNPs-//' | sed 's/-filtered//')
      SNP_COUNT=$(bcftools view -R "$BED_TEMP" "$VCF" | grep -vc '^#')
      DENSITY=$(awk -v n="$SNP_COUNT" -v l="$LENGTH" 'BEGIN{printf "%.10f", n/l}')
      echo -e "$SAMPLE\t$SNP_COUNT\t$LENGTH\t$ID\t$STRAND\t$DENSITY" >> "$OUTPUT_FILE"
    done
  done < "$FAMILY_BED"
done

# === Module 9: Origin summary per gene ===
for FAMILY in ${GENES}
do
  FAMILY_BED="output/${FAMILY}/${FAMILY}.bed"
  [ -s "$FAMILY_BED" ] || continue

  OUTPUT_FILE="output/${FAMILY}/${FAMILY}_origins_Summary.csv"

  ORIGINS=(
    "output/origins/origins.bed"
    "output/origins/miniInitiationZones.bed"
    "output/origins/terminations.bed"
    "output/origins/miniTerminationZones.bed"
    "output/origins/Orc1B.bed"
    "output/origins/Orc1Cdc6.bed"
    "output/origins/genome-without-smAtlas.bed"
    "output/origins/genome-withoutORCs.bed"
  )

  ORIGIN_LABELS=(
    "origins" "miniInitiationZones"
    "terminations" "miniTerminationZones"
    "Orc1B" "Orc1Cdc6" "genome-without-smAtlas.bed"
    "genome-withoutORCs.bed"
  )

  FAMILY_IDS=$(awk '{print $4}' "$FAMILY_BED" | sort -u)
  echo -e "ID\t${ORIGIN_LABELS[*]}" | tr ' ' '\t' > "$OUTPUT_FILE"

  for ID in $FAMILY_IDS
  do
    LINE="$ID"
    for BED in "${ORIGINS[@]}"
    do
      if [[ -s "$BED" ]]; then
        MATCH=$(bedtools intersect -u -a <(awk -v id="$ID" '$4==id' "$FAMILY_BED") -b "$BED")
        [[ -n "$MATCH" ]] && LINE+="\t1" || LINE+="\t0"
      else
        LINE+="\t0"
      fi
    done
    echo -e "$LINE" >> "$OUTPUT_FILE"
  done
done

# === Module 10: SNP density per chromosome ===
OUTPUT_FILE="output/SNPsPerChromosome.csv"
echo -e "Sample\tSNP_Count\tLength\tChromosome\tSNPs_density\tOrigins_density\tOrigin_type" > "$OUTPUT_FILE"

for VCF in output/SNPs-sample*-filtered.vcf.gz
do
  SAMPLE=$(basename "$VCF" | sed 's/SNPs-\(sample[0-9]*\)-filtered.vcf.gz/\1/')

  for CHR in TcChr{1..41}-S
  do
    LENGTH=$(grep -w "$CHR" input/genome-chromoSizes.txt | cut -f2)
    SNP_COUNT=$(bcftools view -r "$CHR" "$VCF" -H | wc -l)
    DENSITY=$(awk -v n="$SNP_COUNT" -v l="$LENGTH" 'BEGIN{printf "%.10f", n/l}')

    for ATLAS in output/origins/*.bed
    do
      TYPE_ORIGINS=$(basename "$ATLAS" .bed)
      ORIGINS_COUNT=$(grep -w "$CHR" "$ATLAS" | wc -l)
      ORIGINS_DENSITY=$(awk -v n="$ORIGINS_COUNT" -v l="$LENGTH" 'BEGIN{printf "%.10f", n/l}')
      echo -e "$SAMPLE\t$SNP_COUNT\t$LENGTH\t$CHR\t$DENSITY\t$ORIGINS_DENSITY\t$TYPE_ORIGINS" >> "$OUTPUT_FILE"
    done
  done
done

# === Module 11: SNPs in origins/terminations summary ===
OUTPUT_FILE="output/SNPs-InOriginsAndTerminationsAndOrcs-Summary.csv"
echo -e "Sample\tSNP_Count\tLength\tGroup\tSNPs_density" > "$OUTPUT_FILE"

for BED in output/origins/*.bed
do
  [ -s "$BED" ] || continue
  REGION=$(basename "$BED" .bed)
  LENGTH=$(awk '{sum += $3 - $2} END {print sum}' "$BED")

  for VCF in output/SNPs-sample*-filtered.vcf.gz
  do
    SAMPLE=$(basename "$VCF" .vcf.gz | sed 's/SNPs-//' | sed 's/-filtered//')
    SNP_COUNT=$(bcftools view -R "$BED" "$VCF" | grep -vc '^#')
    DENSITY=$(awk -v n="$SNP_COUNT" -v l="$LENGTH" 'BEGIN{printf "%.10f", n/l}')
    echo -e "$SAMPLE\t$SNP_COUNT\t$LENGTH\t$REGION\t$DENSITY" >> "$OUTPUT_FILE"
  done
done

# === Module 12: SNPs in multigenic families with Orc1Cdc6 (find + while) ===
# Define output file
OUTPUT="output/SNPsInMultigenicWithOrc1Cdc6.csv"
> "$OUTPUT"  # Clear the output file before appending

# Find all SNP summary CSV files within output subdirectories
find output/*/ -type f -name 'SNPs-*-Summary.csv' | while read -r FILE; do
  echo "Processing $FILE..." >&2  # Print progress to stderr

  # Extract lines that contain both 'with-Orc1Cdc6' and 'without-Orcs',
  # but not 'with-Orc1Cdc6-free' (case-insensitive)
  awk 'BEGIN { IGNORECASE = 1 } (/with-Orc1Cdc6/ || /without-Orcs/) && !/with-Orc1Cdc6-free/ { print }' "$FILE" >> "$OUTPUT"
done

# === Module 13: SNPs in Orc1Cdc6 sites
for GROUP in output/orcs/
do
  GROUP_NAME=$(basename "$GROUP")
  OUTPUT_FILE="${GROUP}/SNPsInOrcsAtMultigenic-Summary.csv"
  echo -e "Sample\tSNP_Count\tLength\tGroup\tSNPs density" \
  > "$OUTPUT_FILE"

  for BED in "${GROUP}"/*.bed
  do
    [ -s "$BED" ] || continue
    REGION=$(basename "$BED" .bed)
    LENGTH=$(awk '{sum += $3 - $2} END {print sum}' "$BED")

    for VCF in output/SNPs-sample*-filtered.vcf.gz
    do
      SAMPLE=$(basename "$VCF" .vcf.gz | sed 's/SNPs-//' | sed 's/-filtered//')
      SNP_COUNT=$(bcftools view -R "$BED" "$VCF" | grep -vc '^#')
      DENSITY=$(awk -v n="$SNP_COUNT" -v l="$LENGTH" 'BEGIN{printf "%.10f", n/l}')
      echo -e "$SAMPLE\t$SNP_COUNT\t$LENGTH\t$REGION\t$DENSITY" >> "$OUTPUT_FILE"
    done
  done
done

# === Module 14: SNP Counting in genome compartment.
# Runs last to include all generated BED files (originals + subtracted versions)

# Create CSV file header
echo -e "SAMPLE\tCOMPARTMENT\tCOUNT\tLENGTH" > output/compartment/compartmentSNPSCount.csv

for COMPARTMENT_BED in output/compartment/*.bed
do
    # Optimization: Calculate compartment size only once (outer loop)
    # FIX: Added missing closing brace '}' in ${COMPARTMENT_BED}
    LENGTH=$(awk -F'\t' '{sum += $3 - $2} END {print sum}' ${COMPARTMENT_BED})
    COMPARTMENT=$(basename ${COMPARTMENT_BED} .bed)

    for SNPS in output/SNPs-sample*-filtered.vcf.gz
    do
        # FIX: Removed space after 'COUNT='.
        # FIX: Corrected 'bcftool' to 'bcftools'.
        COUNT=$(bcftools view -R ${COMPARTMENT_BED} ${SNPS} | grep -v "^#" | wc -l)

        SAMPLE=$(echo ${SNPS} | cut -d '-' -f2)

        # FIX: 'echo' command was missing. Now printing data to file.
        echo -e "${SAMPLE}\t${COMPARTMENT}\t${COUNT}\t${LENGTH}"

    done
done >> output/compartment/compartmentSNPSCount.csv
# FIX: Output redirection moved to the end of the outer loop with '>>' (append) to avoid overwriting.

EOI

# Submit job to SLURM
MAIN_JOB_ID=$(cat log/sbatch-BCFtoolsMain.out)
nice -n 19 /usr/bin/time --verbose --output=log/sbatch-BCFtoolsMain2.time memusg \
--output-to-file log/sbatch-BCFtoolsMain2.memusg \
--time --shell "saveCommand sbatch --parsable \
--dependency=afterany:${MAIN_JOB_ID} \
--nodes 1 \
--ntasks ${NUM_OF_CPUS} \
--mem ${MEMORY_SIZE} \
-o log/slurm-%A.out \
-J BCFtoolsMain2 job/BCFtoolsMain2.slurm \
2> log/sbatch-BCFtoolsMain2.err" | \
tee log/sbatch-BCFtoolsMain2.out

# Info
echo "To monitor jobs: watch -n 1 squeue"

exit 0
