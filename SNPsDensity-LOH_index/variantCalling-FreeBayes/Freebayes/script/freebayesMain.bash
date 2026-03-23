#!/usr/bin/env bash

###################################################################################################
#                                          freebayesMain                                          #
#                                                                                                 #
#  Description:                                                                                   #
#  This script generates and submits a Slurm job that executes the core variant calling pipeline. #
#  It maps Nanopore reads to the reference genome using Winnowmap, processes the resulting        #
#  BAM files (sorting, deduplication with Picard, and filtering), and runs Freebayes in           #
#  parallel across Trypanosoma cruzi chromosomes.                                                 #
#                                                                                                 #
#  Usage:                                                                                         #
#  saveCommand script/freebayesMain.bash                                                          #
#                                                                                                 #
#  Copyleft (ɔ) 2024 by Thiago Andrade Franco                                                     #
#                       <thiago.franco@fundacaobutantan.org.br>                                   #
#                                                                                                 #
#  25/04/2024: First version.                                                                     #
###################################################################################################

# Writing the job script
cat > job/freebayesMain.slurm << EOI
#!/usr/bin/env bash

###################################################################################################
##  Module 1: Read Mapping and BAM Processing                                                    ##
##  Prepares the mapping files for the next steps in the variant calling pipeline                ##
###################################################################################################
# Define variables
GENOME="\${GENOME_FASTA}"
PICARD_JAR=\$(find / -name "picard.jar" 2>/dev/null | head -n 1) # Find the path to picard.jar
K_MER=15

# Check if Picard was found
if [ -z "\${PICARD_JAR}" ]; then
    echo "Error: Picard jar not found. Please check the Picard installation."
    exit 1
fi

# Count k-mers and generate repetitive k-mers files
echo "Counting k-mers in the reference genome..."
meryl count k=\${K_MER} output \${OUTPUT_DIR}/merylDB \${GENOME}
meryl print greater-than distinct=0.9998 \${OUTPUT_DIR}/merylDB > \${OUTPUT_DIR}/repetitive_k\${K_MER}.txt

# Process each sample
for SAMPLE_DIR in \${SAMPLES}; do
    SAMPLE_NAME=\$(basename \${SAMPLE_DIR})
    CONCATENATED_FASTQ="\${OUTPUT_DIR}/\${SAMPLE_NAME}.fastq"

    # Concatenate all FASTQ files for the sample
    echo "Concatenating FASTQ files for sample: \${SAMPLE_NAME}..."
    cat \${INPUT_DIR}/\${SAMPLE_DIR}/*.fastq > \${CONCATENATED_FASTQ}

    # Map the concatenated FASTQ file to the reference genome using Winnowmap
    echo "Mapping \${CONCATENATED_FASTQ} to the reference genome..."
    SAM_OUTPUT="\${OUTPUT_DIR}/\${SAMPLE_NAME}_output.sam"
    BAM_OUTPUT="\${OUTPUT_DIR}/\${SAMPLE_NAME}_output.bam"
    SORTED_BAM_OUTPUT="\${OUTPUT_DIR}/\${SAMPLE_NAME}_sorted.bam"
    DEDUP_BAM_OUTPUT="\${OUTPUT_DIR}/\${SAMPLE_NAME}_dedup.bam"
    METRICS_FILE="\${OUTPUT_DIR}/\${SAMPLE_NAME}_dedup.metrics.txt"
    FILTERED_BAM_OUTPUT="\${OUTPUT_DIR}/\${SAMPLE_NAME}_filtered.bam"

    winnowmap -W \${OUTPUT_DIR}/repetitive_k\${K_MER}.txt -ax map-ont \${GENOME} \${CONCATENATED_FASTQ} > \${SAM_OUTPUT}

    # Convert SAM to BAM
    echo "Converting SAM to BAM for sample: \${SAMPLE_NAME}..."
    samtools view -bS \${SAM_OUTPUT} -o \${BAM_OUTPUT}

    # Sort BAM file by coordinates
    echo "Sorting BAM file by coordinates for sample: \${SAMPLE_NAME}..."
    samtools sort \${BAM_OUTPUT} -o \${SORTED_BAM_OUTPUT}

    # Mark duplicates using Picard
    echo "Marking duplicates for sample: \${SAMPLE_NAME}..."
    java -jar \${PICARD_JAR} MarkDuplicates \\
        -I \${SORTED_BAM_OUTPUT} \\
        -O \${DEDUP_BAM_OUTPUT} \\
        -M \${METRICS_FILE} \\
        -REMOVE_DUPLICATES true \\
        -CREATE_INDEX true

    # Filter BAM by quality (MAPQ >= 30) and exclude secondary/unmapped reads (Flag 0x100)
    echo "Filtering BAM by quality and flags for sample: \${SAMPLE_NAME}..."
    samtools view -b -q 30 -F 0x100 \${DEDUP_BAM_OUTPUT} -o \${FILTERED_BAM_OUTPUT}

    # Index the filtered BAM file
    echo "Indexing the filtered BAM file for sample: \${SAMPLE_NAME}..."
    samtools index \${FILTERED_BAM_OUTPUT}

done

echo "Module 1 processing complete."

###################################################################################################
##  Module 2: Variant Calling with Freebayes                                                     ##
##  Runs Freebayes separated by chromosome for the whole genome of T. cruzi                      ##
###################################################################################################
# Genome reference file
GENOME_FASTA="\${GENOME_FASTA}"

# Create genome index
echo "Creating genome index..."
samtools faidx \${GENOME_FASTA} -o output/genome.fa.fai \\
    2> log/index.err > log/index.out

cp \${GENOME_FASTA} output/genome.fa

# Iterate over each filtered BAM file
echo "Starting Freebayes variant calling per chromosome..."
for BAM_FILE in \$(ls output/*filtered.bam | xargs -I {} basename {}); do
    for i in \$(seq 1 41); do  # Iterate from chromosome 1 to 41
        # Run Freebayes for each chromosome in the background
        freebayes -f output/genome.fa -r TcChr\${i} \\
            output/\${BAM_FILE} > output/\${BAM_FILE%filtered.bam}-In-TcChr\${i}.vcf \\
            2> log/\${BAM_FILE%filtered.bam}-In-TcChr\${i}.err &
    done
    wait  # Wait for all background jobs for the current BAM file to finish before proceeding
done

wait
echo "Module 2 processing complete."

EOI

# Submitting the job to Slurm
/usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-freebayesMain.time memusg \
  --output-to-file log/sbatch-freebayesMain.memusg \
    --time --shell "saveCommand sbatch \
    --nodes 1 \
    --ntasks ${NUM_OF_CPUS} \
    --mem ${MEMORY_SIZE} \
    -o log/slurm-%A.out \
    -J freebayesMain job/freebayesMain.slurm \
    2> log/sbatch-freebayesMain.err | \
    tee log/sbatch-freebayesMain.out"

# Check the status of the running processes
echo "To check the running processes, execute the following command:"
echo "    $> watch -n 1 squeue"
