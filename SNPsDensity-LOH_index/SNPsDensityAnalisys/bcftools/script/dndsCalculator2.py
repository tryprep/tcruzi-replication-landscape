#!/usr/bin/env python3

# ==========================
# dN/dS Calculator (v2 - múltiplos arquivos FASTA via nargs)
# ==========================

import os
import glob
import csv
import re
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

# ==========================
# Argument Parser
# ==========================

parser = argparse.ArgumentParser(description="Calculate dN/dS ratio for annotated SNP TSV files.")
parser.add_argument("--fasta_glob", nargs='+', required=True, help="List of patterns or multifasta files (use shell expansion)")
parser.add_argument("--tsv_glob", required=True, help="Glob pattern for TSV annotation files (e.g., 'output/TSV/*.tsv')")
parser.add_argument("--output", required=True, help="Output CSV file path")

args = parser.parse_args()

# ==========================
# Function: Parse multifasta
# ==========================

def parse_multifastas(patterns):
    """
    Expand multiple glob patterns and parse multifasta files.
    Only sequences containing 'reference' in the header are used.
    """
    ref_sequences = {}
    files = []
    for pattern in patterns:
        matched = glob.glob(pattern)
        files.extend(matched)

    print(f"[INFO] Found {len(files)} multifasta files.")

    for file in files:
        for record in SeqIO.parse(file, "fasta"):
            if "reference" in record.id.lower():
                region_name = os.path.basename(file).replace('-multifasta.fasta', '')
                ref_sequences[region_name] = str(record.seq).upper()
    return ref_sequences

# ==========================
# Function: Parse TSV files
# ==========================

def parse_tsvs(tsv_files):
    """
    Parse TSV annotation files and count synonymous and nonsynonymous SNPs per region and sample.
    Region is defined as the part before '-SNPs-'.
    Sample is captured as 'sampleN' (N=1..5) after '-SNPs-'.
    """
    snp_counts = defaultdict(lambda: defaultdict(lambda: {
        'synonymous_variant': 0,
        'nonsynonymous_variant': 0
    }))

    for tsv in tsv_files:
        basename = os.path.basename(tsv).replace('.tsv', '')

        # Extract Region: everything before '-SNPs-'
        region_match = re.match(r"^(.*)-SNPs-", basename, re.IGNORECASE)
        region = region_match.group(1) if region_match else 'UnknownRegion'

        # Extract Sample: pattern 'sampleN' after '-SNPs-'
        sample_match = re.search(r"-SNPs-(sample\d+)", basename, re.IGNORECASE)
        sample = sample_match.group(1) if sample_match else 'UnknownSample'

        line_count = 0
        snp_line_count = 0

        with open(tsv) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line_count += 1
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue

                ann = fields[4]
                ann_parts = ann.split('|')
                if len(ann_parts) < 2:
                    continue

                effect = ann_parts[1]
                if 'synonymous_variant' in effect:
                    snp_counts[region][sample]['synonymous_variant'] += 1
                    snp_line_count += 1
                elif any(term in effect for term in ['missense_variant', 'start_lost', 'stop_gained', 'stop_lost']):
                    snp_counts[region][sample]['nonsynonymous_variant'] += 1
                    snp_line_count += 1

        print(f"[DEBUG] {os.path.basename(tsv)}: {line_count} lines, {snp_line_count} SNPs counted.")

    return snp_counts

# ==========================
# Function: Compute sites
# ==========================

def compute_sites(seq):
    """
    Compute the number of synonymous and nonsynonymous sites in a coding DNA sequence.
    """
    table = CodonTable.unambiguous_dna_by_name["Standard"]

    seq = seq.upper().replace('T', 'U')
    syn_sites = 0.0
    nonsyn_sites = 0.0

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) != 3 or 'N' in codon:
            continue

        try:
            ref_aa = Seq(codon.replace('U', 'T')).translate(table=table)
        except:
            continue

        syn_count = 0
        nonsyn_count = 0

        for pos in range(3):
            for nt in ['A', 'U', 'C', 'G']:
                if nt == codon[pos]:
                    continue

                mutated_codon = list(codon)
                mutated_codon[pos] = nt
                mutated_codon = ''.join(mutated_codon)

                try:
                    mut_aa = Seq(mutated_codon.replace('U', 'T')).translate(table=table)
                except:
                    continue

                if ref_aa == mut_aa:
                    syn_count += 1
                else:
                    nonsyn_count += 1

        syn_sites += syn_count / 9.0
        nonsyn_sites += nonsyn_count / 9.0

    return syn_sites, nonsyn_sites

# ==========================
# Main Execution
# ==========================

# Load input files
tsv_files = glob.glob(args.tsv_glob)
print(f"[INFO] Found {len(tsv_files)} TSV files.")

ref_sequences = parse_multifastas(args.fasta_glob)
snp_counts = parse_tsvs(tsv_files)

print("[INFO] Calculating dN/dS ratios...\n")

output_rows = []

for region, seq in ref_sequences.items():
    S_sites, N_sites = compute_sites(seq)

    if region not in snp_counts:
        print(f"[WARNING] No SNP data found for region: {region}")
        continue

    for sample, counts in snp_counts[region].items():
        syn = counts['synonymous_variant']
        nonsyn = counts['nonsynonymous_variant']

        dN = nonsyn / N_sites if N_sites > 0 else 'NA'
        dS = syn / S_sites if S_sites > 0 else 'NA'

        if isinstance(dN, str) or isinstance(dS, str) or dS == 0:
            dnds = 'NA' if dS == 0 else 'Inf'
        else:
            dnds = round(dN / dS, 5)

        output_rows.append({
            'Region': region,
            'Sample': sample,
            'Syn_SNPs': syn,
            'NonSyn_SNPs': nonsyn,
            'S_sites': round(S_sites, 3),
            'N_sites': round(N_sites, 3),
            'dN': round(dN, 5) if isinstance(dN, float) else dN,
            'dS': round(dS, 5) if isinstance(dS, float) else dS,
            'dN/dS': dnds
        })

# ==========================
# Save Output CSV
# ==========================

print(f"[INFO] Writing output to {args.output}...")

with open(args.output, 'w', newline='') as csvfile:
    fieldnames = ['Region', 'Sample', 'Syn_SNPs', 'NonSyn_SNPs', 'S_sites', 'N_sites', 'dN', 'dS', 'dN/dS']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for row in output_rows:
        writer.writerow(row)

print("[INFO] dN/dS calculation complete.\n")
