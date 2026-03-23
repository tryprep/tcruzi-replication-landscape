#!/usr/bin/env python3
import sys
import bisect
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List

# --- Data Structures ---

@dataclass
class Region:
    chrom: str
    start: int
    end: int
    strand: str

# --- File Reading Functions ---

def read_chrom_sizes(path: str) -> Dict[str, int]:
    """Reads chromosome sizes file (format: chrom \t size)."""
    chrom_sizes = {}
    try:
        with open(path, 'r') as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
    except Exception as e:
        print(f"Error reading chromosome sizes: {e}")
        sys.exit(1)
    return chrom_sizes

def read_bed_to_dict(path: str) -> Dict[str, List[Region]]:
    """Reads a BED file and organizes regions by chromosome for fast lookup."""
    data = {}
    try:
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 6:
                    chrom = parts[0]
                    if chrom not in data:
                        data[chrom] = []
                    data[chrom].append(Region(
                        chrom=chrom,
                        start=int(parts[1]),
                        end=int(parts[2]),
                        strand=parts[5]
                    ))
    except Exception as e:
        print(f"Warning reading {path}: {e}")
    return data

def read_snps_by_sample(snp_files: List[str]) -> Dict[str, Dict[str, List[int]]]:
    """
    Reads multiple SNP files and stores them per sample and per chromosome.
    Returns: { 'sample1': { 'chr1': [pos1, pos2], ... }, 'sample2': ... }
    """
    all_samples_snps = {}
    for i, path in enumerate(snp_files):
        sample_name = f"sample{i+1}"
        all_samples_snps[sample_name] = {}
        try:
            with open(path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        chrom, pos = parts[0], int(parts[1])
                        if chrom not in all_samples_snps[sample_name]:
                            all_samples_snps[sample_name][chrom] = []
                        all_samples_snps[sample_name][chrom].append(pos)

            # Sorting is essential for binary search (O(log N))
            for chrom in all_samples_snps[sample_name]:
                all_samples_snps[sample_name][chrom].sort()
            print(f"[*] Loaded SNPs for {sample_name} from {path}")
        except Exception as e:
            print(f"Error loading SNP file {path}: {e}")
    return all_samples_snps

# --- Processing Logic ---

def get_strand_from_poly(chrom: str, start: int, end: int, poly_ref: Dict[str, List[Region]]) -> str:
    """
    Assigns strand if the window is fully contained within a polycistron.
    Otherwise returns 'GAP'.
    """
    if chrom in poly_ref:
        for poly in poly_ref[chrom]:
            # Rule: Window must be strictly inside the polycistron limits
            if start >= poly.start and end <= poly.end:
                return poly.strand
    return "GAP"

def count_snps(sample_snps: Dict[str, List[int]], chrom: str, start: int, end: int) -> int:
    """Counts SNPs in range [start, end) using efficient binary search."""
    if chrom not in sample_snps:
        return 0
    positions = sample_snps[chrom]
    idx_start = bisect.bisect_left(positions, start)
    idx_end = bisect.bisect_left(positions, end)
    return idx_end - idx_start

def write_line(out_f, region_id, sample, r_type, chrom, start, end, group, strand, snp_count):
    """Calculates density and writes a tab-separated line to the output file."""
    length = end - start
    density = (snp_count / length) if length > 0 else 0.0
    line = (f"{region_id}\t{sample}\t{r_type}\t{chrom}\t"
            f"{start}\t{end}\t{length}\t{snp_count}\t"
            f"{density:.6f}\t{group}\t{strand}\n")
    out_f.write(line)

def process_windows(out_f, sample_name, sample_snps, site, poly_ref, chrom_sizes, group):
    """Generates the central site and 30kb flanking windows (upstream/downstream)."""
    start, end = site.start, site.end
    chrom = site.chrom
    chr_limit = chrom_sizes.get(chrom, 0)
    region_id = f"ORC_{chrom}_{start}_{end}_{group}"

    WINDOW_SIZE = 2000
    MAX_DISTANCE = 30000

    # 1. Central Site
    s_strand = get_strand_from_poly(chrom, start, end, poly_ref)
    s_count = count_snps(sample_snps, chrom, start, end)
    write_line(out_f, region_id, sample_name, "site", chrom, start, end, group, s_strand, s_count)

    # 2. Upstream Flanks (15 windows of 2kb)
    for offset in range(WINDOW_SIZE, MAX_DISTANCE + 1, WINDOW_SIZE):
        w_start = max(0, start - offset)
        w_end = start - (offset - WINDOW_SIZE)
        if w_end <= 0: break # Out of chromosome bounds

        u_count = count_snps(sample_snps, chrom, w_start, w_end)
        u_strand = get_strand_from_poly(chrom, w_start, w_end, poly_ref)
        write_line(out_f, region_id, sample_name, f"upstream_{offset}bp", chrom, w_start, w_end, group, u_strand, u_count)

    # 3. Downstream Flanks (15 windows of 2kb)
    for offset in range(WINDOW_SIZE, MAX_DISTANCE + 1, WINDOW_SIZE):
        w_start = end + (offset - WINDOW_SIZE)
        w_end = min(chr_limit, end + offset)
        if w_start >= chr_limit: break # Out of chromosome bounds

        d_count = count_snps(sample_snps, chrom, w_start, w_end)
        d_strand = get_strand_from_poly(chrom, w_start, w_end, poly_ref)
        write_line(out_f, region_id, sample_name, f"downstream_{offset}bp", chrom, w_start, w_end, group, d_strand, d_count)

# --- Main Execution ---

def main():
    if len(sys.argv) < 6:
        print("Usage: python3 orc_analyzer.py <chromSizes> <Poly_BED> <Input_BED_Dir_or_File> <Output_TSV> <SNP_File1> <SNP_File2> ...")
        sys.exit(1)

    # Input arguments
    chrom_sizes_path = sys.argv[1]
    poly_bed_path = sys.argv[2]
    input_path = Path(sys.argv[3])
    output_path = sys.argv[4]
    snp_files = sys.argv[5:] # Files for sample1, sample2, etc.

    print("[*] Initializing analysis...")

    # Load reference data
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    poly_ref = read_bed_to_dict(poly_bed_path)
    all_snps_data = read_snps_by_sample(snp_files)

    # Determine which BED files to process (Origins, Initiation Zones, etc.)
    if input_path.is_dir():
        bed_files = sorted(list(input_path.glob("*.bed")))
    else:
        bed_files = [input_path]

    # Process and write results
    with open(output_path, 'w') as out_file:
        # TSV Header
        out_file.write("ID\tsample\tregion\tchrom\tstart\tend\tlength\tSNPs_count\tSNPs_density\tGroup\tStrand\n")

        for b_file in bed_files:
            group_name = b_file.stem
            print(f"[*] Processing group: {group_name}")

            # Load ORC sites for the current group
            sites_to_analyze = read_bed_to_dict(str(b_file))

            for chrom, sites in sites_to_analyze.items():
                for site in sites:
                    # Iterate through each sample to ensure independent SNP counting
                    for sample_name, sample_snps in all_snps_data.items():
                        process_windows(out_file, sample_name, sample_snps, site, poly_ref, chrom_sizes, group_name)

    print(f"[*] Done! Results saved to: {output_path}")

if __name__ == "__main__":
    main()
