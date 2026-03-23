/*
 * filter_vcf.cpp
 *
 * This program processes annotated VCF files (output from variant annotation tools such as SnpEff)
 * and filters variant lines to retain only those where all annotations (ANN=) share the same effect type.
 * Header lines are preserved. The program reads from multiple VCF files (sample1 to sample4),
 * performs the filtering, and writes the output to new files.
 *
 * This implementation improves performance over the equivalent shell script by using C++.
 * Developed for large-scale VCF processing in bioinformatics pipelines.
 *
 * Author: [Thiago Franco]
 * Date: [15/04/2025]
 */

// filter_vcf.cpp
// This program filters annotated VCF files by retaining only variants with consistent annotation effects.
// It keeps entries with a single effect type, or multiple effects if they are all in a defined whitelist or share the same gene.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

using namespace std;

const set<string> WHITELIST = {
    "intergenic_region",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
    "non_coding_transcript_variant",
    "synonymous_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant"
};

vector<string> split(const string &s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

bool should_keep_variant(const string &ann_values) {
    set<string> effects;
    set<string> gene_ids;

    vector<string> annotations = split(ann_values, ',');
    for (const string &ann : annotations) {
        vector<string> fields = split(ann, '|');
        if (fields.size() < 5) continue;

        string effect = fields[1];
        string gene_id = fields[4];

        effects.insert(effect);
        gene_ids.insert(gene_id);
    }

    if (effects.size() == 1) return true; // Only one effect type

    bool all_in_whitelist = all_of(effects.begin(), effects.end(), [](const string &eff) {
        return WHITELIST.count(eff) > 0;
    });

    if (all_in_whitelist) return true;
    if (gene_ids.size() == 1) return true; // Same gene

    return false;
}

void process_vcf(const string &input_path, const string &output_path) {
    ifstream infile(input_path);
    ofstream outfile(output_path);
    string line;

    while (getline(infile, line)) {
        if (line.rfind("#", 0) == 0) {
            outfile << line << '\n';
            continue;
        }

        vector<string> fields = split(line, '\t');
        if (fields.size() < 8) continue;

        string info = fields[7];
        size_t ann_pos = info.find("ANN=");
        if (ann_pos == string::npos) continue;

        size_t end_pos = info.find(';', ann_pos);
        string ann_values = info.substr(ann_pos + 4, end_pos - ann_pos - 4);

        if (should_keep_variant(ann_values)) {
            outfile << line << '\n';
        }
    }

    infile.close();
    outfile.close();
}

int main() {
    for (int i = 1; i <= 4; ++i) {
        string input_file = "output/SNPs-sample" + to_string(i) + "-annotated.vcf";
        string output_file = "output/SNPs-sample" + to_string(i) + "-filtered.vcf";

        cout << "🔍 Processing " << input_file << "..." << endl;
        process_vcf(input_file, output_file);
        cout << "✅ Filtered VCF written to " << output_file << endl;
    }

    return 0;
}
