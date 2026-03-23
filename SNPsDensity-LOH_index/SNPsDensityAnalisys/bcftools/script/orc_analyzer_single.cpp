#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <filesystem>

/**
 * Structure to store SNP positions per chromosome.
 * Using a vector of ints for sorted positions to allow binary search.
 */
using SNPMap = std::unordered_map<std::string, std::vector<int>>;

/**
 * Structure for genomic regions (BED format style).
 */
struct Region {
    std::string chrom;
    int start;
    int end;
    std::string strand;
};

/**
 * Reads chromosome sizes from a tab-delimited file (chrom\tsize).
 * Returns a map for O(1) average lookup.
 */
std::unordered_map<std::string, int> read_chrom_sizes(const std::string& path) {
    std::unordered_map<std::string, int> chrom_sizes;
    std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open chromosome sizes file: " << path << std::endl;
        return chrom_sizes;
    }

    std::string chrom;
    int size;
    while (infile >> chrom >> size) {
        chrom_sizes[chrom] = size;
    }
    return chrom_sizes;
}

/**
 * Reads SNPs from multiple files and sorts positions for efficient binary search.
 */
SNPMap read_snps(const std::vector<std::string>& snpFiles) {
    SNPMap snps;
    for (const auto& path : snpFiles) {
        std::ifstream infile(path);
        if (!infile.is_open()) {
            std::cerr << "Warning: Could not open SNP file: " << path << std::endl;
            continue;
        }

        std::string chrom;
        int pos;
        while (infile >> chrom >> pos) {
            snps[chrom].push_back(pos);
        }
    }

    // Sorting is required for std::lower_bound (binary search) efficiency
    for (auto& [chr, positions] : snps) {
        std::sort(positions.begin(), positions.end());
    }
    return snps;
}

/**
 * Reads BED file. Handles cases where strand or other columns might be missing.
 * Default strand is set to "." if not provided.
 */
std::vector<Region> read_bed(const std::string& path) {
    std::vector<Region> regions;
    std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open BED file: " << path << std::endl;
        return regions;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip empty lines and headers

        std::istringstream iss(line);
        std::string chrom, name, score, strand_val = ".";
        int start, end;

        // BED format requires at least chrom, start, end
        if (!(iss >> chrom >> start >> end)) continue;

        // Try to read optional columns; if they don't exist, strand stays as "."
        iss >> name >> score >> strand_val;

        regions.push_back({chrom, start, end, strand_val});
    }
    return regions;
}

/**
 * Checks for intersection (overlap) with annotations to determine the strand.
 */
std::string get_strand(const std::vector<Region>& genome_bed, const std::string& chrom, int start, int end) {
    for (const auto& region : genome_bed) {
        if (region.chrom == chrom) {
            // Intersection logic: (StartA < EndB) and (EndA > StartB)
            if (start < region.end && end > region.start) {
                return region.strand;
            }
        }
    }
    return "GAP";
}

/**
 * Counts SNPs within a range [start, end) using binary search (O(log N)).
 */
int count_snps(const SNPMap& snps, const std::string& chrom, int start, int end) {
    auto it = snps.find(chrom);
    if (it == snps.end()) return 0;

    const auto& positions = it->second;
    auto it_start = std::lower_bound(positions.begin(), positions.end(), start);
    auto it_end = std::lower_bound(positions.begin(), positions.end(), end);

    return static_cast<int>(std::distance(it_start, it_end));
}

/**
 * Writes a single row to the output TSV file.
 */
void write_region(std::ofstream& out, const std::string& id, const std::string& sample, const std::string& region_type,
                  const std::string& chrom, int start, int end, const std::string& strand, int snp_count) {
    int length = end - start;
    double snp_density = (length > 0) ? static_cast<double>(snp_count) / length : 0.0;

    out << id << "\t" << sample << "\t" << region_type << "\t" << chrom << "\t"
    << start << "\t" << end << "\t" << length << "\t" << snp_count << "\t"
    << snp_density << "\t" << strand << "\n";
                  }

                  /**
                   * Generates sliding windows around a central site (e.g., origin).
                   * Handles chromosome boundaries safely to avoid out_of_range errors.
                   */
                  void create_windows(std::ofstream& out, const SNPMap& snps, const std::vector<Region>& genome_bed,
                                      const std::unordered_map<std::string, int>& chrom_sizes,
                                      const Region& site, const std::string& sample_name) {

                      // BUG FIX: Use .find() instead of .at() to prevent crashes if chrom is missing
                      auto it_chrom = chrom_sizes.find(site.chrom);
                      if (it_chrom == chrom_sizes.end()) {
                          static bool warned = false;
                          if (!warned) {
                              std::cerr << "Warning: Chromosome " << site.chrom << " not found in sizes file. Skipping site." << std::endl;
                              warned = true; // Avoid flooding the console
                          }
                          return;
                      }

                      int start = site.start;
                      int end = site.end;
                      int chr_max = it_chrom->second;
                      std::string id = "ORC_" + site.chrom + "_" + std::to_string(start) + "_" + std::to_string(end);

                      // List to store windows: {start_position, label}
                      std::vector<std::pair<int, std::string>> window_list = {{start, "site"}};

                      const int MAX_DIST = 30000;   // 30kb limit
                      const int WINDOW_SIZE = 2000; // 2kb window size

                      // --- UPSTREAM WINDOWS ---
                      for (int offset = WINDOW_SIZE; offset <= MAX_DIST; offset += WINDOW_SIZE) {
                          int w_start = start - offset;
                          if (w_start >= 0) {
                              window_list.push_back({w_start, "upstream_" + std::to_string(offset) + "bp"});
                          } else {
                              // Handle remaining space if start of chromosome is reached
                              int rem_up = start % WINDOW_SIZE;
                              if (rem_up > 0) {
                                  window_list.push_back({0, "upstream_rem_" + std::to_string(rem_up) + "bp"});
                              }
                              break;
                          }
                      }

                      // --- DOWNSTREAM WINDOWS ---
                      for (int offset = WINDOW_SIZE; offset <= MAX_DIST; offset += WINDOW_SIZE) {
                          int w_start = end + (offset - WINDOW_SIZE);
                          int w_end = w_start + WINDOW_SIZE;
                          if (w_end <= chr_max) {
                              window_list.push_back({w_start, "downstream_" + std::to_string(offset) + "bp"});
                          } else {
                              // Handle remaining space if end of chromosome is reached
                              int dist_to_end = chr_max - end;
                              int rem_down = dist_to_end % WINDOW_SIZE;
                              if (rem_down > 0) {
                                  window_list.push_back({chr_max - rem_down, "downstream_rem_" + std::to_string(rem_down) + "bp"});
                              }
                              break;
                          }
                      }

                      // --- PROCESS AND OUTPUT ---
                      for (const auto& [w_start, w_type] : window_list) {
                          int w_end;
                          if (w_type == "site") {
                              w_end = end;
                          } else if (w_type.find("rem") != std::string::npos) {
                              if (w_type.find("upstream") != std::string::npos)
                                  w_end = w_start + (start % WINDOW_SIZE);
                              else
                                  w_end = w_start + ((chr_max - end) % WINDOW_SIZE);
                          } else {
                              w_end = w_start + WINDOW_SIZE;
                          }

                          // Final boundary safety check
                          if (w_end > chr_max) w_end = chr_max;
                          if (w_start < 0) continue;

                          int snp_count = count_snps(snps, site.chrom, w_start, w_end);
                          std::string strand = get_strand(genome_bed, site.chrom, w_start, w_end);
                          write_region(out, id, sample_name, w_type, site.chrom, w_start, w_end, strand, snp_count);
                      }
                                      }

                                      int main(int argc, char* argv[]) {
                                          if (argc < 6) {
                                              std::cerr << "Usage: " << argv[0] << " <chromSizes> <Annotations_BED> <Origins_BED> <OutputFile> <SNP_Files...>" << std::endl;
                                              return 1;
                                          }

                                          // Loading reference data
                                          auto chrom_sizes = read_chrom_sizes(argv[1]);
                                          auto genome_bed = read_bed(argv[2]);
                                          auto origins = read_bed(argv[3]);

                                          // Collecting SNP files from arguments
                                          std::vector<std::string> snp_files;
                                          for (int i = 5; i < argc; ++i) snp_files.push_back(argv[i]);
                                          auto snps = read_snps(snp_files);

                                          std::ofstream outfile(argv[4]);
                                          if (!outfile.is_open()) {
                                              std::cerr << "Error: Could not open output file " << argv[4] << std::endl;
                                              return 1;
                                          }

                                          // TSV Header
                                          outfile << "ID\tsample\tregion\tchrom\tstart\tend\tlength\tSNPs_count\tSNPs_density\tStrand\n";

                                          // Process each origin site for 4 samples
                                          for (const auto& site : origins) {
                                              for (int i = 1; i <= 4; ++i) {
                                                  create_windows(outfile, snps, genome_bed, chrom_sizes, site, "sample" + std::to_string(i));
                                              }
                                          }

                                          std::cout << "Processing complete. Results saved to: " << argv[4] << std::endl;
                                          return 0;
                                      }
