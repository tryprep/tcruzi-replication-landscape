import pysam
import sys
import gzip

# --- Input and Output files (passed as command-line arguments) ---
vcf_map_file = sys.argv[1]
bam_file = sys.argv[2]
output_file = sys.argv[3]
ref_genome_file = sys.argv[4] # We also need the reference genome

print(f"Reading SNP map from: {vcf_map_file}")
print(f"Reading alignments from: {bam_file}")
print(f"Using reference genome: {ref_genome_file}")
print(f"Saving results to: {output_file}")

# Open the BAM file for reading
samfile = pysam.AlignmentFile(bam_file, "rb")

# Open the output file for writing
with open(output_file, 'w') as out_f:
    # Open the input VCF (handles both .vcf and .vcf.gz)
    opener = gzip.open if vcf_map_file.endswith('.gz') else open
    with opener(vcf_map_file, 'rt') as vcf_f:
        for line in vcf_f:
            # Skip the VCF header lines
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref_allele = fields[3]
            alt_allele = fields[4]

            # Dictionary to count the alleles
            allele_counts = {ref_allele: 0, alt_allele: 0}

            # pysam uses 0-based coordinates, while VCF is 1-based
            # Therefore, the query position must be pos - 1
            query_pos = pos - 1

            # Iterate over the reads covering the position of interest
            # truncate=True ensures that only the base at the exact position is considered
            for pileupcolumn in samfile.pileup(chrom, query_pos, query_pos + 1, truncate=True):
                if pileupcolumn.pos == query_pos:
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            if base == ref_allele:
                                allele_counts[ref_allele] += 1
                            elif base == alt_allele:
                                allele_counts[alt_allele] += 1

            # Write the counts to the output file in the expected format
            ref_count = allele_counts[ref_allele]
            alt_count = allele_counts[alt_allele]

            # Only write to the file if there is coverage at the site
            if ref_count + alt_count > 0:
                out_f.write(f"{chrom}\t{pos}\t{ref_allele}\t{alt_allele}\t{ref_count},{alt_count}\n")

samfile.close()
print("Process completed successfully!")
