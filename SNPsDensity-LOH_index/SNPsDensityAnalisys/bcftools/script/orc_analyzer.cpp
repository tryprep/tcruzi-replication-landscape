#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <climits>

using SNPMap = std::unordered_map<std::string, std::vector<int>>;

// Structure representing a genomic region (e.g., from a BED file)
struct Region {
    std::string chrom;   // Chromosome name
    int start;           // Start coordinate (0-based)
    int end;             // End coordinate
    std::string id;      // Unique identifier
    std::string strand;  // Strand information (+, -, or .)
};

// Reads chromosome names and their total lengths from a file
std::unordered_map<std::string, int> read_chrom_sizes(const std::string& path) {
    std::unordered_map<std::string, int> chrom_sizes;
    std::ifstream infile(path);
    std::string chrom;
    int size;
    while (infile >> chrom >> size) {
        chrom_sizes[chrom] = size;
    }
    return chrom_sizes;
}

// Reads SNP positions from multiple files and sorts them for efficient searching
SNPMap read_snps(const std::vector<std::string>& snpFiles) {
    SNPMap snps;
    for (const auto& path : snpFiles) {
        std::ifstream infile(path);
        std::string chrom;
        int pos;
        while (infile >> chrom >> pos) {
            snps[chrom].push_back(pos);
        }
    }
    // Sorting allows for binary search (O(log n)) instead of linear search
    for (auto& [chr, positions] : snps) {
        std::sort(positions.begin(), positions.end());
    }
    return snps;
}

// Reads genomic regions from a BED-formatted file
std::vector<Region> read_bed(const std::string& path) {
    std::vector<Region> regions;
    std::ifstream infile(path);
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string chrom, strand, score;
        int start, end;
        // Expecting standard BED columns
        if (!(iss >> chrom >> start >> end >> strand >> score >> strand)) continue;
        regions.push_back({chrom, start, end, "", strand});
    }
    return regions;
}

// Checks which strand a specific window belongs to based on a reference GFF/BED
std::string get_strand(const std::vector<Region>& genome_bed, const std::string& chrom, int start, int end) {
    for (const auto& region : genome_bed) {
        if (region.chrom == chrom && start >= region.start && end <= region.end) {
            return region.strand;
        }
    }
    return "GAP"; // Returns GAP if the region is not covered in the reference
}

// Counts SNPs within a range [start, end) using binary search for performance
int count_snps(const SNPMap& snps, const std::string& chrom, int start, int end) {
    if (snps.find(chrom) == snps.end()) return 0;
    const auto& positions = snps.at(chrom);
    auto it_start = std::lower_bound(positions.begin(), positions.end(), start);
    auto it_end = std::lower_bound(positions.begin(), positions.end(), end);
    return std::distance(it_start, it_end);
}

// Formats and writes the calculated data for a single region/window to the output file
void write_region(std::ofstream& out, const std::string& id, const std::string& sample, const std::string& region_type,
                  const std::string& chrom, int start, int end, const std::string& group,
                  const std::string& strand, int snp_count) {
    int length = end - start;
    // SNP density is calculated as $Density = \frac{count}{length}$
    double snp_density = (length > 0) ? static_cast<double>(snp_count) / length : 0.0;

    out << id << "\t" << sample << "\t" << region_type << "\t" << chrom << "\t"
    << start << "\t" << end << "\t" << length << "\t" << snp_count << "\t"
    << snp_density << "\t" << group << "\t" << strand << "\n";
                  }

                  // Generates 2000bp windows up to a 30kb limit around a target site
                  void create_windows(std::ofstream& out, const SNPMap& snps, const std::vector<Region>& genome_bed,
                                      const std::unordered_map<std::string, int>& chrom_sizes,
                                      const Region& site, const std::string& sample_name, const std::string& group) {

                      int start = site.start;
                      int end = site.end;
                      int chr_limit = chrom_sizes.at(site.chrom);
                      std::string id = "ORC_" + site.chrom + "_" + std::to_string(start) + "_" + std::to_string(end) + "_" + group;

                      const int WINDOW_SIZE = 2000;
                      const int MAX_DISTANCE = 30000; // Limit counting to 30kb from the site

                      // 1. Process the central site itself
                      int site_snps = count_snps(snps, site.chrom, start, end);
                      std::string site_strand = get_strand(genome_bed, site.chrom, start, end);
                      write_region(out, id, sample_name, "site", site.chrom, start, end, group, site_strand, site_snps);

                      // 2. Process Upstream regions (moving backwards from 'start')
                      for (int offset = WINDOW_SIZE; offset <= MAX_DISTANCE; offset += WINDOW_SIZE) {
                          int w_start = start - offset;
                          int w_end = start - (offset - WINDOW_SIZE);

                          if (w_start < 0) {
                              // Handle cases near the beginning of the chromosome
                              if (w_end > 0) {
                                  int count = count_snps(snps, site.chrom, 0, w_end);
                                  std::string strand = get_strand(genome_bed, site.chrom, 0, w_end);
                                  write_region(out, id, sample_name, "upstream_" + std::to_string(offset) + "bp", site.chrom, 0, w_end, group, strand, count);
                              }
                              break; // Stop if we reach the start of the chromosome
                          } else {
                              int count = count_snps(snps, site.chrom, w_start, w_end);
                              std::string strand = get_strand(genome_bed, site.chrom, w_start, w_end);
                              write_region(out, id, sample_name, "upstream_" + std::to_string(offset) + "bp", site.chrom, w_start, w_end, group, strand, count);
                          }
                      }

                      // 3. Process Downstream regions (moving forwards from 'end')
                      for (int offset = WINDOW_SIZE; offset <= MAX_DISTANCE; offset += WINDOW_SIZE) {
                          int w_start = end + (offset - WINDOW_SIZE);
                          int w_end = end + offset;

                          if (w_end > chr_limit) {
                              // Handle cases near the end of the chromosome
                              if (w_start < chr_limit) {
                                  int count = count_snps(snps, site.chrom, w_start, chr_limit);
                                  std::string strand = get_strand(genome_bed, site.chrom, w_start, chr_limit);
                                  write_region(out, id, sample_name, "downstream_" + std::to_string(offset) + "bp", site.chrom, w_start, chr_limit, group, strand, count);
                              }
                              break; // Stop if we reach the end of the chromosome
                          } else {
                              int count = count_snps(snps, site.chrom, w_start, w_end);
                              std::string strand = get_strand(genome_bed, site.chrom, w_start, w_end);
                              write_region(out, id, sample_name, "downstream_" + std::to_string(offset) + "bp", site.chrom, w_start, w_end, group, strand, count);
                          }
                      }
                                      }

                                      int main(int argc, char* argv[]) {
                                          if (argc < 6) {
                                              std::cerr << "Usage: " << argv[0] << " <chromSizes> <GFF_BED> <BED_DIR> <OUTPUT_FILE> <SNP_TSV_FILES...>" << std::endl;
                                              return 1;
                                          }

                                          // Initialize file paths from command line arguments
                                          std::string chrom_sizes_file = argv[1];
                                          std::string gff_file = argv[2];
                                          std::string bed_dir = argv[3];
                                          std::string output_file = argv[4];
                                          std::vector<std::string> snp_files(argv + 5, argv + argc);

                                          // Load data into memory
                                          auto chrom_sizes = read_chrom_sizes(chrom_sizes_file);
                                          auto snps = read_snps(snp_files);
                                          auto genome_bed = read_bed(gff_file);

                                          // Prepare output file and write header
                                          std::ofstream outfile(output_file);
                                          outfile << "ID\tsample\tregion\tchrom\tstart\tend\tlength\tSNPs_count\tSNPs_density\tGroup\tStrand\n";

                                          // Iterate through all BED files in the specified directory
                                          for (const auto& bed_file : std::filesystem::directory_iterator(bed_dir)) {
                                              auto regions = read_bed(bed_file.path().string());
                                              std::string group = bed_file.path().stem().string();

                                              for (const auto& site : regions) {
                                                  // Processing for 4 hypothetical samples as per original logic
                                                  for (int i = 1; i <= 4; ++i) {
                                                      create_windows(outfile, snps, genome_bed, chrom_sizes, site, "sample" + std::to_string(i), group);
                                                  }
                                              }
                                          }

                                          return 0;
                                      }
