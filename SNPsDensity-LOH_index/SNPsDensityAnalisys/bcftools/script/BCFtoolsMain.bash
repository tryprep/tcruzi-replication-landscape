#!/usr/bin/env bash

###################################################################################################
#                                         BCFtoolsMain                                            #
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
#  saveCommand script/BCFtoolsMain.bash                                                           #
#                                                                                                 #
#                                                                                                 #
#  Copyleft (ɔ) 2025 by Thiago Andrade Franco                                                     #
#                       <thiago.franco@fundacaobutantan.org.br>                                   #
#                                                                                                 #
#  23/04/2025: First version.                                                                     #
###################################################################################################

# Writing the job script to submit to Slurm
cat > job/BCFtoolsMain.slurm << 'EOI'
#!/usr/bin/env bash

###################################################################################################
## STEP 1: Generate BED Files for Multigenic Families                                            ##
###################################################################################################
# Path to the chromosome sizes file (required for sorting and boundary checks)
GENOME_SIZES="input/genome-chromoSizes.txt"

# Process common families
for MULTIGENIC in GP63 RHS DGF-1 MASP
do
    mkdir -p "output/${MULTIGENIC}"
    grep "${MULTIGENIC}" input/genome.gff | \
    awk 'BEGIN {OFS="\t"} $3 == "protein_coding_gene" {match($9, /ID=([^;]*)/, arr); print $1, $4, $5, arr[1], $6, $7}' \
    > "output/${MULTIGENIC}/${MULTIGENIC}.bed"
done

# Process Histone
mkdir -p output/HIST
grep "histone" input/genome.gff | \
awk 'BEGIN {OFS="\t"} $3 == "protein_coding_gene" {match($9, /ID=([^;]*)/, arr); print $1, $4, $5, arr[1], $6, $7}' \
> output/HIST/HIST.bed

# Process Trans-Sialidase
mkdir -p output/TS
grep "trans-sialidase" input/genome.gff | \
awk 'BEGIN {OFS="\t"} $3 == "protein_coding_gene" {match($9, /ID=([^;]*)/, arr); print $1, $4, $5, arr[1], $6, $7}' \
> output/TS/TS.bed

# Process Mucin (excluding mucin_associated)
mkdir -p output/MUC
grep -i "mucin" input/genome.gff | grep -vi "mucin_associated" | \
awk 'BEGIN {OFS="\t"} $3 == "protein_coding_gene" {match($9, /ID=([^;]*)/, arr); print $1, $4, $5, arr[1], $6, $7}' \
> output/MUC/MUC.bed

###################################################################################################
## STEP 2: Preprocess Input BEDs for Origins and Terminations                                    ##
###################################################################################################

# Run C++ program for BED processing
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/processBed.cpp.time \
memusg --output-to-file=log/processBed.cpp.memusg \
script/processBed \
    1> log/processBed.out \
    2> log/processBed.err

###################################################################################################
## STEP 3: Utility Functions for Intersect and Subtract Operations                               ##
###################################################################################################
# Intersect BED with input and process
intersect_and_process() {
    local gene="$1"
    local input_bed="$2"
    local output_suffix="$3"

    bedtools intersect -wo -a "output/${gene}/${gene}.bed" -b "${input_bed}" | \
    awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6}' | \
    sort -u -k1,1 -k2,2n -k3,3n \
    > "output/${gene}/${gene}-${output_suffix}.bed"
}

# Intersect with multiple BEDs
intersect_multiple_and_process() {
    local gene="$1"
    local output_suffix="$2"
    shift 2
    local b_files=("$@")

    bedtools intersect -wo -a "output/${gene}/${gene}.bed" $(printf -- '-b %s ' "${b_files[@]}") | \
    awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6}' | \
    sort -u -k1,1 -k2,2n -k3,3n \
    > "output/${gene}/${gene}-${output_suffix}.bed"
}

# Subtract BED regions
subtract_and_process() {
    local gene="$1"
    local output_suffix="$2"
    shift 2
    local b_files=("$@")

    bedtools subtract -A -a "output/${gene}/${gene}.bed" $(printf -- '-b %s ' "${b_files[@]}") | \
    awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6}' | \
    sort -u -k1,1 -k2,2n -k3,3n \
    > "output/${gene}/${gene}-${output_suffix}.bed"
}

###################################################################################################
## STEP 4: Process Multigenic Families for Replication Features                                  ##
###################################################################################################
# Note: The ${GENES} variable must be previously exported in the environment/driver script.
for MULTIGENIC in ${GENES}; do

    echo -e "\n=== Processing gene family: ${MULTIGENIC} ==="

    ## smAtlas
    intersect_and_process "${MULTIGENIC}" output/origins/origins.bed with-origins
    intersect_and_process "${MULTIGENIC}" output/origins/miniInitiationZones.bed with-miniInitiationZones
    subtract_and_process "${MULTIGENIC}" without-smAtlas \
        "output/${MULTIGENIC}/${MULTIGENIC}-with-origins.bed" \
        "output/${MULTIGENIC}/${MULTIGENIC}-with-miniInitiationZones.bed"

    ## Terminations
    intersect_and_process "${MULTIGENIC}" output/origins/terminations.bed with-terminations
    intersect_and_process "${MULTIGENIC}" output/origins/miniTerminationZones.bed with-miniTerminationZones
    subtract_and_process "${MULTIGENIC}" without-terminations \
        "output/${MULTIGENIC}/${MULTIGENIC}-with-terminations.bed" \
        "output/${MULTIGENIC}/${MULTIGENIC}-with-miniTerminationZones.bed"

    ## ORCs
    intersect_and_process "${MULTIGENIC}" output/origins/Orc1Cdc6.bed with-Orc1Cdc6
    intersect_and_process "${MULTIGENIC}" output/origins/Orc1B.bed with-Orc1B
    intersect_multiple_and_process "${MULTIGENIC}" with-Orcs \
        output/origins/Orc1Cdc6.bed \
        output/origins/Orc1B.bed
    subtract_and_process "${MULTIGENIC}" without-Orcs \
        output/origins/Orc1Cdc6.bed \
        output/origins/Orc1B.bed
done

###################################################################################################
## STEP 5: Generate Random Control Sequences                                                     ##
###################################################################################################
# Generate Base Genome Windows (<10% N)
if [ -z "$GENOME_FASTA" ]; then
    echo "Error: GENOME_FASTA variable is not set."
    exit 1
fi

# Create raw tiling windows (2kb and 6kb)
for SIZE in 2000 6000; do
    echo "Generating ${SIZE}bp windows..."
    bedtools makewindows -g ${GENOME_SIZES} -w ${SIZE} \
    > output/genome-${SIZE}.bed

    # Extract FASTA sequences for N-counting
    bedtools getfasta -fi ${GENOME_FASTA} -bed output/genome-${SIZE}.bed \
    > output/genome-${SIZE}.fasta
done

# Filter FASTA for sequences with <= 10% 'N'
for RANDOM_FASTA in $(ls output/genome-*.fasta | xargs -I {} basename {}); do
    echo "Filtering ${RANDOM_FASTA} for <10% 'N'..."
    awk -v max_N=0.1 '
    BEGIN {RS=">"; FS="\n"}
    NR > 1 {
        header = $1
        seq = substr($0, index($0, "\n") + 1)
        gsub("\n", "", seq)
        N_count = gsub("N", "N", seq) + gsub("n", "n", seq) # Count both N and n
        seq_length = length(seq)

        # Avoid division by zero for empty sequences
        if (seq_length > 0 && (N_count / seq_length) <= max_N) {
            # Reformat sequence for standard FASTA output (60 chars per line)
            print ">" header
            for (i = 1; i <= seq_length; i += 60) {
                print substr(seq, i, 60)
            }
        }
    }' output/${RANDOM_FASTA} \
    > output/${RANDOM_FASTA%.fasta}-filtered.fasta
done

# Convert filtered FASTA headers to BED format
for RANDOM_FASTA_FILTERED in $(ls output/genome*-filtered.fasta | xargs -I {} basename {}); do
    echo "Converting ${RANDOM_FASTA_FILTERED} headers to BED..."

    awk '
    BEGIN { OFS = "\t" }
    /^>/ {
        gsub("\r", "");
        if (match($0, />([^:]+):([0-9]+)-([0-9]+)/, arr)) {
            print arr[1], arr[2], arr[3]
        }
    }' output/${RANDOM_FASTA_FILTERED} \
    > output/${RANDOM_FASTA_FILTERED%.fasta}.bed

done

# Sort filtered BED files
for RANDOM_FILTERED_BED in $(ls output/*-filtered.bed | xargs -I {} basename {}); do
    echo "Sorting ${RANDOM_FILTERED_BED}..."
    bedtools sort -g "${GENOME_SIZES}" \
        -i output/${RANDOM_FILTERED_BED} \
        > output/${RANDOM_FILTERED_BED%-filtered.bed}-sorted.bed
done

## Generate Master Origin-Free Control Sets
echo "Creating Master Origin exclusion list..."
cat output/origins/Orc1B.bed \
    output/origins/Orc1Cdc6.bed \
    output/origins/origins.bed \
    output/origins/miniInitiationZones.bed \
    output/origins/terminations.bed \
    output/origins/miniTerminationZones.bed | \
awk 'OFS="\t" {print $1, $2, $3}' | \
bedtools sort -i - -g "${GENOME_SIZES}" \
    > output/origins/ALL_ORIGINS_MASTER.bed

# Extend master exclusion list by 20kb
bedtools flank -i output/origins/ALL_ORIGINS_MASTER.bed -g "${GENOME_SIZES}" -b 20000 > output/origins/ALL_ORIGINS_MASTER-extend20kb.bed

# Generate final Origin-Free windows
echo "Generating final 2000bp Origin-Free set..."
bedtools subtract -A \
    -a output/genome-2000-sorted.bed \
    -b output/origins/ALL_ORIGINS_MASTER-extend20kb.bed \
    > output/origins/genome-2000-ALL-ORIGIN-FREE.bed

echo "Generating final 6000bp Origin-Free set..."
bedtools subtract -A \
    -a output/genome-6000-sorted.bed \
    -b output/origins/ALL_ORIGINS_MASTER-extend20kb.bed \
    > output/origins/genome-6000-ALL-ORIGIN-FREE.bed

# Cleanup
echo "Cleaning up intermediate files..."

rm output/genome-2000.bed
rm output/genome-6000.bed
rm output/genome-2000.fasta
rm output/genome-6000.fasta
rm output/genome-2000-filtered.fasta
rm output/genome-6000-filtered.fasta
rm output/genome-2000-filtered.bed
rm output/genome-6000-filtered.bed
rm output/origins/ALL_ORIGINS_MASTER.bed
rm output/origins/ALL_ORIGINS_MASTER-extend20kb.bed

echo "Done. Final files:"
echo "output/origins/genome-2000-ALL-ORIGIN-FREE.bed"
echo "output/origins/genome-6000-ALL-ORIGIN-FREE.bed"

###################################################################################################
## STEP 6: Generate Orc1Cdc6 BED files Control Datasets                                          ##
###################################################################################################
# Create output directory if it doesn't exist
OUTPUT_DIR="output/orcs"
mkdir -p "${OUTPUT_DIR}"

# Define Orc1Cdc6 BED file
ORC_BED="output/origins/Orc1Cdc6.bed"

# Temporary file to accumulate all intersections across families
ALL_HITS_TMP="${OUTPUT_DIR}/orc1cdc6_ALL_HITS.tmp"
> "${ALL_HITS_TMP}"  # Empty the file at the beginning

# === Loop over each multigenic family ===
for MULTIGENIC in ${GENES}
do
    # Define input and output paths
    GENE_BED="output/${MULTIGENIC}/${MULTIGENIC}.bed"
    OUTPUT_INTERSECT="${OUTPUT_DIR}/orc1cdc6In${MULTIGENIC}.bed"

    # Check if both required files exist
    if [[ -f "${ORC_BED}" && -f "${GENE_BED}" ]]; then
        echo "[INFO] Processing ${MULTIGENIC}..."

        # Intersect Orc1Cdc6 sites with current multigenic family
        bedtools intersect -wo -a "${ORC_BED}" -b "${GENE_BED}" | \
            awk 'OFS="\t" {print $1, $2, $3}' | \
            sort -k1,1 -k2,2n -k3,3n | uniq > "${OUTPUT_INTERSECT}"

        # Accumulate all intersecting regions
        cat "${OUTPUT_INTERSECT}" >> "${ALL_HITS_TMP}"

    else
        echo "[ERROR] Missing input file(s) for ${MULTIGENIC}. Skipping..."
        [[ ! -f "${ORC_BED}" ]] && echo "[ERROR] Missing file: ${ORC_BED}"
        [[ ! -f "${GENE_BED}" ]] && echo "[ERROR] Missing file: ${GENE_BED}"
    fi
done

# === Identify Orc1Cdc6 regions that do NOT intersect any multigenic family ===

# Sort and remove duplicates from accumulated hits
ALL_HITS_SORTED="${OUTPUT_DIR}/orc1cdc6_ALL_HITS.sorted.bed"
sort -k1,1 -k2,2n -k3,3n "${ALL_HITS_TMP}" | uniq > "${ALL_HITS_SORTED}"

# Subtract all intersecting regions from the full Orc1Cdc6 set
OUTPUT_NOHIT="${OUTPUT_DIR}/orc1cdc6NotInAnyMultigenic.bed"
bedtools intersect -v -a "${ORC_BED}" -b "${ALL_HITS_SORTED}" > "${OUTPUT_NOHIT}"

echo "[INFO] Non-overlapping Orc1Cdc6 regions saved to: ${OUTPUT_NOHIT}"

# Cleanup temporary files
rm -f "${ALL_HITS_TMP}" "${ALL_HITS_SORTED}"

# Orc1Cdc6 in active and inactive origins
for prefix in output/origins; do
    # Active origins: Orc1Cdc6 minus Dormant
    if [[ -f "${prefix}/Orc1Cdc6.bed" && -f "${prefix}/dormant.bed" ]]; then
        bedtools subtract -A \
            -a ${prefix}/Orc1Cdc6.bed \
            -b ${prefix}/dormant.bed \
        | sort -k1,1 -k2,2n -k3,3n \
        > ${prefix}/Orc1Cdc6-onActiveOrigins.bed

        # Inactive origins: just copy Dormant
        cp ${prefix}/dormant.bed ${prefix}/Orc1Cdc6-onInactiveOrigins.bed
    fi
done

###################################################################################################
## STEP 7: Classification of Origins Based on Orc1Cdc6 Overlap                                   ##
###################################################################################################
# 1. Concatenate origins from atlas 2.
cat output/origins/miniInitiationZones.bed output/origins/origins.bed |
awk 'OFS="\t" {print $1, $2, $3}' |\
bedtools sort -i - -g input/genome-chromoSizes.txt \
> output/origins/smAtlas.bed

# 2. Identify Mini Initiation Zones overlapping with Orc1Cdc6 sites.
bedtools intersect -wo \
    -a <(awk 'OFS="\t" {print $1, $2, $3}' output/origins/smAtlas.bed) \
    -b <(awk 'OFS="\t" {print $1, $2, $3}' output/origins/Orc1Cdc6.bed) | \
    sort -u -k1,1 -k2,2n -k3,3n | \
    awk 'OFS="\t" {print $1, $2, $3}' \
    > output/origins/smAtlas-withOrc1cdc6.bed

# 3. Identify Mini Initiation Zones NOT overlapping with Orc1Cdc6 sites.
bedtools subtract -A \
    -a <(awk 'OFS="\t" {print $1, $2, $3}' output/origins/smAtlas.bed) \
    -b output/origins/smAtlas-withOrc1cdc6.bed | \
    sort -u -k1,1 -k2,2n -k3,3n | \
    awk 'OFS="\t" {print $1, $2, $3}' \
    > output/origins/smAtlas-withoutOrc1cdc6.bed

###################################################################################################
## STEP 8: Creating BED files of origins overlapping with genome compartments                    ##
###################################################################################################
# Create output directories if they don't exist
mkdir -p output/compartment
mkdir -p output/origins

# 1. Define the GFF files (Compartments) explicitly
GFF_FILES=("input/GpDR.gff" "input/core.gff" "input/disruptive.gff")

# 2. Capture a static list of the ORIGINAL .bed files currently present.
# This prevents the loop from iterating over the new files created during execution.
ORIGIN_FILES=(output/origins/*.bed)

# Outer Loop: Iterate through the captured list of BED files
for ORIGINS in "${ORIGIN_FILES[@]}"; do

    # Safety check: ensure the file actually exists
    [ -e "$ORIGINS" ] || continue

    # Define basename for the origins file
    BASENAME_ORIGINS=$(basename "${ORIGINS}" .bed)

    # Inner Loop: Iterate through each Compartment GFF for the current BED file
    for COMPARTMENT in "${GFF_FILES[@]}"; do

        # Define basename for the compartment file
        BASENAME_GFF=$(basename "${COMPARTMENT}" .gff)

        echo "Processing: ${BASENAME_ORIGINS} vs ${BASENAME_GFF}"

        # Intersect the current origins file with the current compartment GFF
        bedtools intersect -wo \
            -a "${ORIGINS}" \
            -b <( awk 'OFS="\t" {print $1, $4, $5}' "${COMPARTMENT}" ) | \
            awk 'OFS="\t" {print $1, $2, $3}' | \
            sort -u -k1,1 -k2,2n -k3,3n \
            > "output/compartment/${BASENAME_ORIGINS}-${BASENAME_GFF}.bed"

    done
done

###################################################################################################
## STEP 9: Convert Compartments to BED and Subtract smAtlas Origins                              ##
###################################################################################################

# --- PART 1: Convert GFF to BED ---
for COMPARTMENT in input/core.gff input/disruptive.gff input/GpDR.gff
do
    # Check if input file exists
    if [ ! -f "$COMPARTMENT" ]; then
        echo "Warning: File $COMPARTMENT not found."
        continue
    fi

    COMPARTMENT_BASENAME=$(basename ${COMPARTMENT} .gff)
    COMPARTMENT_DIR="output/compartment"

    awk 'OFS="\t" {print $1, $4, $5}' ${COMPARTMENT} | sort -u -k1,1 -k2,2n -k3,3n \
    > ${COMPARTMENT_DIR}/${COMPARTMENT_BASENAME}.bed
done

# --- PART 2: Bedtools Subtractions (smAtlas) ---

# Check if smAtlas file exists to avoid errors
if [ -f "output/origins/smAtlas.bed" ]; then
    bedtools subtract -A -a output/compartment/core.bed -b output/origins/smAtlas.bed > output/compartment/core-Without-smAtlas.bed
    bedtools subtract -A -a output/compartment/disruptive.bed -b output/origins/smAtlas.bed > output/compartment/disruptive-Without-smAtlas.bed
    bedtools subtract -A -a output/compartment/GpDR.bed -b output/origins/smAtlas.bed > output/compartment/GpDR-Without-smAtlas.bed
fi

echo "All BED file processing completed successfully."

EOI

# Submitting the job to Slurm
/usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-BCFtoolsMain.time memusg \
  --output-to-file log/sbatch-BCFtoolsMain.memusg \
  --time --shell "saveCommand sbatch --parsable \
  --nodes 1 \
  --ntasks ${NUM_OF_CPUS} \
  --mem ${MEMORY_SIZE} \
  -o log/slurm-%A.out \
  -J BCFtoolsMain job/BCFtoolsMain.slurm \
  2> log/sbatch-BCFtoolsMain.err |
  tee log/sbatch-BCFtoolsMain.out"

# Helpful tip to monitor running processes
echo "To check the running processes, execute one of the following commands:"
echo "    $> watch -n 1 sequeue"

exit 0
