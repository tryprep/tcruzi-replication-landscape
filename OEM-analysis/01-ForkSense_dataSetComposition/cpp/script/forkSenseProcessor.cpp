#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <regex>
#include <cstdlib>

// Function to split a string by spaces
std::vector<std::string> split_string(const std::string& str) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;
    while (ss >> token) {
        result.push_back(token);
    }
    return result;
}

// Function to process forkSense files
void process_forkSense_files(const std::string& sub_directory) {
    int count = 1;
    for (const auto& entry : std::filesystem::directory_iterator("input/" + sub_directory)) {
        if (entry.path().extension() == ".forkSense") {
            std::ifstream infile(entry.path());
            std::ofstream outfile;
            std::string line;
            std::string outdir = "output/" + sub_directory;
            std::filesystem::create_directories(outdir);

            while (std::getline(infile, line)) {
                if (line[0] == '>') {
                    if (outfile.is_open()) {
                        outfile.close();
                    }
                    outfile.open(outdir + "/FRAGMENT_" + std::to_string(count) + ".forkSense");
                    count++;
                }
                if (outfile.is_open()) {
                    outfile << line << "\n";
                }
            }
            if (outfile.is_open()) {
                outfile.close();
            }
        }
    }
}

// Function to create probability fork dataset
void create_fork_probability_dataset(const std::string& sub_directory) {
    std::ofstream outFile("output/probabilityForkSense.csv", std::ios_base::app);
    std::regex header_regex(R"(>(\S+)\s(\S+)\s.*\s(\S+)$)");

    for (const auto& entry : std::filesystem::directory_iterator("output/" + sub_directory)) {
        if (entry.path().stem().string().find("FRAGMENT_") != std::string::npos) {
            std::ifstream processed_file(entry.path());
            std::string line, header_line;

            if (std::getline(processed_file, header_line)) {
                std::smatch match;
                if (std::regex_search(header_line, match, header_regex)) {
                    std::string id = match[1].str();
                    std::string chromo = match[2].str();
                    std::string strand = match[3].str();

                    while (std::getline(processed_file, line)) {
                        std::stringstream ss(line);
                        std::string start_str, fork_left_str, fork_right_str;
                        ss >> start_str >> fork_left_str >> fork_right_str;

                        int start = std::stoi(start_str);
                        int end = start + 1;
                        double probability_fork_left = std::stod(fork_left_str);
                        double probability_fork_right = std::stod(fork_right_str);

                        outFile << id << "\t" << chromo << "\t" << start << "\t" << end << "\t"
                                << probability_fork_left << "\t" << probability_fork_right << "\t"
                                << strand << "\n";
                    }
                }
            }
        }
    }
    outFile.close();
}

int main() {
    // List of sample directories
    std::string samples = std::getenv("SAMPLES"); // Assuming SAMPLES is an environment variable
    std::vector<std::string> sample_array = split_string(samples);

    // Process forkSense files for each sub-directory
    for (const auto& sub_directory : sample_array) {
        process_forkSense_files(sub_directory);
    }

    // Create the fork probability dataset
    for (const auto& sub_directory : sample_array) {
        create_fork_probability_dataset(sub_directory);
    }

    return 0;
}
