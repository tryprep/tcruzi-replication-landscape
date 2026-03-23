#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

struct BEDRecord {
    std::string chrom;
    int start;
    int end;
};

struct GFFRecord {
    std::string chrom;
    int start;
    int end;
    std::string id;
    char strand;
};

// Function to split strings by delimiter
std::vector<std::string> split(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Load GFF as BED
std::vector<GFFRecord> loadGFF(const std::string &filename) {
    std::vector<GFFRecord> records;
    std::ifstream infile(filename);
    std::string line;

    while (std::getline(infile, line)) {
        auto fields = split(line, '\t');
        if (fields.size() >= 6) {
            std::string id = ".";
            size_t pos = fields[8].find("ID=");
            if (pos != std::string::npos) {
                size_t end_pos = fields[8].find(";", pos);
                id = fields[8].substr(pos + 3, end_pos - pos - 3);
            }
            records.push_back({fields[0], std::stoi(fields[3]), std::stoi(fields[4]), id, fields[6][0]});
        }
    }
    infile.close();
    return records;
}

// Find the closest GFF record when no overlap is found
GFFRecord findClosestRecord(const std::string &chrom, int start, int end, const std::vector<GFFRecord> &gffRecords) {
    GFFRecord closest = {".", -1, -1, ".", '.'};
    int minDistance = std::numeric_limits<int>::max();

    for (const auto &record : gffRecords) {
        if (record.chrom == chrom) {
            int distance = std::min(abs(record.start - end), abs(record.end - start));
            if (distance < minDistance) {
                minDistance = distance;
                closest = record;
            }
        }
    }
    return closest;
}

// Process BED files
void processBED(const std::string &bedFile, const std::vector<GFFRecord> &gffRecords, const std::string &outputFile) {
    std::ifstream infile(bedFile);
    std::ofstream outfile(outputFile);
    std::string line;

    while (std::getline(infile, line)) {
        auto fields = split(line, '\t');
        if (fields.size() >= 3) {
            std::string chrom = fields[0];
            int start = std::stoi(fields[1]);
            int end = std::stoi(fields[2]);

            int bestOverlap = 0;
            std::string bestId = ".";
            char bestStrand = '.';

            for (const auto &record : gffRecords) {
                if (record.chrom == chrom && start <= record.end && end >= record.start) {
                    int overlap = std::min(end, record.end) - std::max(start, record.start);
                    if (overlap > bestOverlap) {
                        bestOverlap = overlap;
                        bestId = record.id;
                        bestStrand = record.strand;
                    }
                }
            }

            // If no overlap is found, find the closest record
            if (bestOverlap == 0) {
                GFFRecord closest = findClosestRecord(chrom, start, end, gffRecords);
                bestId = closest.id;
                bestStrand = closest.strand;
            }

            outfile << chrom << '\t' << start << '\t' << end << '\t' << bestId << '\t' << "0" << '\t' << bestStrand << '\n';
        }
    }

    infile.close();
    outfile.close();
}

int main() {
    // Create output directory
    fs::create_directories("output/origins");

    // Load GFF data
    std::vector<GFFRecord> gffRecords = loadGFF("input/genome.gff");

    // Process each BED file
    for (const auto &entry : fs::directory_iterator("input")) {
        if (entry.path().extension() == ".bed") {
            std::string bedFile = entry.path().string();
            std::string outputFile = "output/origins/" + entry.path().filename().string();
            processBED(bedFile, gffRecords, outputFile);
        }
    }

    return 0;
}
