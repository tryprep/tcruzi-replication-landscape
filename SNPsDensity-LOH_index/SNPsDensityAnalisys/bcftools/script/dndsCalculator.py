#!/usr/bin/env python3

# =====================
# dN/dS Calculator (Region-specific and Sample-specific)
# =====================

import os
import glob
import csv
import re
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

# =====================
# ARGUMENT PARSER
# =====================

parser = argparse.ArgumentParser(description="Calculate dN/dS ratio for origin/termination regions.")
parser.add_argument("--fasta_glob", required=True, help="Glob pattern for multifasta files (e.g., 'output/origins/*multifasta.fasta')")
parser.add_argument("--tsv_glob", required=True, help="Glob pattern for TSV annotation files (e.g., 'output/TSV/*.tsv')")
parser.add_argument("--output", required=True, help="Output CSV file to store results")

args = parser.parse_args()

# =====================
# FUNCTIONS
# =====================

def parse_multifastas(multifasta_files):
    """
    Parse multifasta files to extract reference sequences for each region.
    """
    ref_sequences = {}
    for file in multifasta_files:
        for record in SeqIO.parse(file, "fasta"):
            if "reference" in record.id.lower():
                region_name = os.path.basename(file).replace('-multifasta.fasta', '')
                ref_sequences[region_name] = str(record.seq).upper()
    return ref_sequences


def parse_tsvs(tsv_files):
    """
    Parse TSV annotation files and count synonymous and nonsynonymous SNPs per region and sample.
    Region is defined as everything before '-SNPs', and sample is captured from 'sampleN' after '-SNPs'.
    """
    snp_counts = defaultdict(lambda: defaultdict(lambda: {
        'synonymous_variant': 0,
        'nonsynonymous_variant': 0
    }))

    for tsv in tsv_files:
        basename = os.path.basename(tsv).replace('.tsv', '')

        # Extract Region: everything before '-SNPs'
        region_match = re.match(r"^(.*)-SNPs-", basename, re.IGNORECASE)
        region = region_match.group(1) if region_match else 'UnknownRegion'

        # Extract Sample: pattern 'sampleN' (N=1..5) after '-SNPs-'
        sample_match = re.search(r"-SNPs-(sample\d+)", basename, re.IGNORECASE)
        sample = sample_match.group(1) if sample_match else 'UnknownSample'

        with open(tsv) as f:
            for line in f:
                if line.startswith('#'):
                    continue
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
                elif any(term in effect for term in ['missense_variant', 'start_lost', 'stop_gained', 'stop_lost']):
                    snp_counts[region][sample]['nonsynonymous_variant'] += 1

    return snp_counts


def compute_sites(seq):
    """
    Compute the number of synonymous and nonsynonymous sites in a reference sequence.
    This is done by simulating all possible point mutations in each codon.
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

# =====================
# MAIN EXECUTION
# =====================

# Load input files
multifasta_files = glob.glob(args.fasta_glob)
tsv_files = glob.glob(args.tsv_glob)
output_file = args.output

print("\n[INFO] Parsing reference multifasta sequences...")
ref_sequences = parse_multifastas(multifasta_files)

print("[INFO] Parsing TSV files for SNP counts...")
snp_counts = parse_tsvs(tsv_files)

print("[INFO] Calculating dN/dS ratios...\n")

output_rows = []

# Compute dN/dS for each region and sample
for region, seq in ref_sequences.items():
    S_sites, N_sites = compute_sites(seq)

    if region not in snp_counts:
        print(f"[WARNING] No SNP data found for region: {region}")
        continue

    for sample, counts in snp_counts[region].items():
        syn = counts['synonymous_variant']
        nonsyn = counts['nonsynonymous_variant']

        # Compute dN and dS
        dN = nonsyn / N_sites if N_sites > 0 else 'NA'
        dS = syn / S_sites if S_sites > 0 else 'NA'

        # Compute dN/dS ratio
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

# =====================
# SAVE OUTPUT CSV
# =====================

print(f"[INFO] Writing output to {output_file}...")

with open(output_file, 'w', newline='') as csvfile:
    fieldnames = ['Region', 'Sample', 'Syn_SNPs', 'NonSyn_SNPs', 'S_sites', 'N_sites', 'dN', 'dS', 'dN/dS']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for row in output_rows:
        writer.writerow(row)

print("[INFO] dN/dS calculation complete.\n")
