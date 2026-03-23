#!/usr/bin/env bash
###################################################################################################
#                                                                                                 #
#                                        BCFtoolsMain4                                            #
#                                                                                                 #
#   This script automates the analysis of Loss of Heterozygosity (LOH) across defined genomic     #
#   regions. It is designed to run on a high-performance computing (HPC) cluster using SLURM.     #
#                                                                                                 #
#   Pipeline Steps:                                                                               #
#   1. Consensus Map Creation: Merges VCF files from multiple samples and identifies high-        #
#      confidence heterozygous sites that are shared across all samples. This map serves as the   #
#      basis for LOH analysis.                                                                    #
#   2. Per-Sample LOH Calculation: For each sample, the script calculates allele depths at the    #
#      consensus sites and computes a "LOH metric" for each site, which quantifies the deviation  #
#      from the expected 1:1 allelic ratio.                                                       #
#   3. Regional LOH Scoring: The per-site LOH metrics are then mapped onto user-provided          #
#      genomic regions (BED files), and an average LOH score is calculated for each region.       #
#                                                                                                 #
#  Usage:                                                                                         #
#  saveCommand script/BCFtoolsMain3.bash                                                          #
#                                                                                                 #
#                                                                                                 #
#  Copyleft (ɔ) 2025 by Thiago Andrade Franco                                                     #
#                       <thiago.franco@fundacaobutantan.org.br>                                   #
#                                                                                                 #
#  23/04/2025: First version.                                                                     #
###################################################################################################
## -----------------------------------------------------------------------------------------------
## CONFIGURATION
## -----------------------------------------------------------------------------------------------
# Define the samples to be processed
SAMPLES_TO_PROCESS="1 2 3 4"

# Log file from the previous job (BCFtoolsMain3), containing its job ID.
# NOTE: This method is fragile. A better approach is to capture the job ID
# directly upon submission of the previous job.
PREVIOUS_JOB_LOG="log/sbatch-BCFtoolsMain3.out"

## -----------------------------------------------------------------------------------------------
## SLURM SCRIPT GENERATION
## -----------------------------------------------------------------------------------------------
# This section writes the main analysis logic into a SLURM script file.
# Using a 'here document' (<<'EOI') allows for a clean, multi-line script definition.

echo "Generating SLURM job script: job/BCFtoolsMain4.slurm"

cat > job/BCFtoolsMain4.slurm <<'EOI'
#!/usr/bin/env bash

# This command ensures that the script will exit immediately if any command fails.
set -e

echo "===== STEP 1: Create the consensus heterozygous map ====="
# Create a list of VCF files for the specified samples
ls output/SNPs-sample{1..4}-filtered.vcf.gz > output/vcf_list.txt

# Merge the VCF files
# The --force-samples flag is important if samples have different sets of calls
bcftools merge -l output/vcf_list.txt -O z -o output/merged_SNPs.vcf.gz --force-samples

# Index the merged file
bcftools index output/merged_SNPs.vcf.gz

# Filter to find sites that are heterozygous in all 4 samples
bcftools view -i 'COUNT(GT="het")>=4' output/merged_SNPs.vcf.gz -O z -o output/consensus_heterozygous_map.vcf.gz

# Index the final map
bcftools index output/consensus_heterozygous_map.vcf.gz
echo "Consensus map created successfully!"
echo ""

echo "===== STEP 1.5: Index BAM files (Critical Prerequisite) ====="
# This loop ensures that every BAM file has an index (.bai).
# This is required for tools like pysam to perform fast lookups of genomic regions.
for i in $(seq 1 4); do
    echo "--- Indexing BAM file for Sample ${i} ---"
    samtools index input/${i}_filtered.bam
done
echo "All BAM files are indexed."
echo ""

echo "===== STEP 2: Pre-calculate LOH metrics for each sample (Performance Optimization) ====="
# This loop calculates the per-site LOH metrics for each sample ONCE.
# These intermediate files will be reused for each BED file analysis.
for i in $(seq 1 4); do
    echo "--- Pre-calculating metrics for Sample ${i} ---"

    # --- Sub-step A: Count allele depths using the Python script (memory efficient) ---
    python script/count_alleles.py \
        output/consensus_heterozygous_map.vcf.gz \
        input/${i}_filtered.bam \
        output/sample_${i}.allele_depth.tsv \
        input/genome.fasta

    # --- Sub-step B: Calculate the LOH metric per site ---
    awk 'BEGIN{OFS="\t"} {if ($5 != ".") {split($5, ad, ","); ref_d=ad[1]; alt_d=ad[2]; total_d=ref_d+alt_d; if (total_d > 0) {baf=alt_d/total_d; loh_m=sqrt((baf - 0.5)^2); print $1, $2, $2+1, loh_m;}}}' \
        output/sample_${i}.allele_depth.tsv > output/sample_${i}-loh_metric_per_site.bed
done

# Creating the BEDGRAPH file
for BED in output/sample*-loh_metric_per_site.bed
do
    sort -k1,1 -k2,2n -k3,3n ${BED} > ${BED%.bed}.bedgraph
done

echo "All per-sample metrics have been pre-calculated."
echo ""

# Calculate Average LOH (Smart Mean Profile)
# Aligning samples and calculating the mean (ignoring missing data)..."
bedtools unionbedg \
    -i output/sample_1-loh_metric_per_site.bedgraph \
       output/sample_2-loh_metric_per_site.bedgraph \
       output/sample_3-loh_metric_per_site.bedgraph \
       output/sample_4-loh_metric_per_site.bedgraph \
    -filler "NA" | \
    awk 'BEGIN{OFS="\t"} {
        sum=0; count=0;
        # Iterate through columns 4, 5, 6, and 7 (representing the 4 samples)
        for(i=4; i<=7; i++) {
            # Only sum values that are actual data (not "NA")
            if($i != "NA") {
                sum += $i;
                count++;
            }
        }
        # Calculate the mean only if at least one sample has data for this site
        # This creates an adaptive mean: sum / valid_sample_count
        if(count > 0) {
            print $1, $2, $3, sum/count
        }
    }' \
    > output/average_LOH_consensus.bedgraph

echo "Average LOH consensus file created: output/average_LOH_consensus.bedgraph"
echo ""

echo "===== STEP 3: Cleanup ====="
echo "Removing intermediate files..."
rm -f output/sample_*.allele_depth.tsv
rm -f output/sample_*.loh_metric_per_site.bed
echo "Cleanup complete."
echo ""

echo "===== Pipeline finished successfully! ====="

EOI

## -----------------------------------------------------------------------------------------------
## JOB SUBMISSION
## -----------------------------------------------------------------------------------------------
# This section submits the generated script to the SLURM scheduler.

# Retrieve the job ID from the previous step's log file.
MAIN_JOB_ID=$(cat ${PREVIOUS_JOB_LOG})

echo "Submitting BCFtoolsMain4 job to SLURM, dependent on job ${MAIN_JOB_ID}..."

# Submit the job using a combination of system tools for monitoring.
# 'saveCommand' is assumed to be a local wrapper script/alias.
/usr/bin/nice -n 19 /usr/bin/time \
    --verbose \
    --output=log/sbatch-BCFtoolsMain4.time memusg \
    --output-to-file log/sbatch-BCFtoolsMain4.memusg \
    --time --shell "saveCommand sbatch --parsable \
    --dependency=afterany:${MAIN_JOB_ID} \
    --nodes=1 \
    --ntasks=${NUM_OF_CPUS} \
    --mem=${MEMORY_SIZE} \
    -o log/slurm-%A.out \
    -J BCFtoolsMain4 job/BCFtoolsMain4.slurm \
    2> log/sbatch-BCFtoolsMain4.err |
    tee log/sbatch-BCFtoolsMain4.out"

# Provide user feedback on how to monitor the submitted job.
echo ""
echo "Job submitted. To check the running processes, execute one of the following commands:"
echo "    $> watch -n 1 squeue -u \$USER"
echo "    $> squeue -j \$(cat log/sbatch-BCFtoolsMain4.out)"

exit 0
