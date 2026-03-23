#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <regex>

namespace fs = std::filesystem;

void split(const std::string& s, std::vector<std::string>& tokens) {
    std::stringstream ss(s);
    std::string token;
    while (ss >> token) {
        tokens.push_back(token);
    }
}

void processFile(const std::string& subDirectory, const std::string& fileName) {
    std::ifstream inputFile("input/" + subDirectory + "/" + fileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: input/" << subDirectory << "/" << fileName << std::endl;
        return;
    }

    std::string line;
    int count = 1;
    std::ofstream outFile;

    while (std::getline(inputFile, line)) {
        if (line.rfind(">", 0) == 0) { // Line starts with '>'
            if (outFile.is_open()) {
                outFile.close();
            }
            std::string outputFileName = "output/" + subDirectory + "/FRAGMENT_" + std::to_string(count) + ".brduDetect";
            outFile.open(outputFileName);
            count++;
        }
        if (outFile.is_open()) {
            outFile << line << std::endl;
        }
    }

    if (outFile.is_open()) {
        outFile.close();
    }
}

void createForkProbabilityDataset(const std::vector<std::string>& sampleArray) {
    std::string line;
    std::ofstream outputFile("output/probability-BrdU.csv", std::ios::app);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file" << std::endl;
        return;
    }

    for (const auto& outSubDirectory : sampleArray) {
        for (const auto& entry : fs::directory_iterator("output/" + outSubDirectory)) {
            std::string processedFile = entry.path().filename().string();
            std::ifstream inputFile("output/" + outSubDirectory + "/" + processedFile);
            std::string headerLine;

            // Get the header line (first line)
            std::getline(inputFile, headerLine);
            std::regex headerRegex(R"(>([^ ]+) ([^ ]+) .* (rev|fwd))");
            std::smatch match;

            if (std::regex_search(headerLine, match, headerRegex)) {
                std::string id = match[1].str();
                std::string chromo = match[2].str();
                std::string strand = match[3].str();

                while (std::getline(inputFile, line)) {
                    std::istringstream lineStream(line);
                    std::string startStr, probabilityBrdU;
                    lineStream >> startStr >> probabilityBrdU;

                    // Trim spaces
                    startStr.erase(std::remove_if(startStr.begin(), startStr.end(), ::isspace), startStr.end());
                    probabilityBrdU.erase(std::remove_if(probabilityBrdU.begin(), probabilityBrdU.end(), ::isspace), probabilityBrdU.end());

                    // Check if START is a valid number
                    if (std::regex_match(startStr, std::regex(R"(^[0-9]+$)"))) {
                        int start = std::stoi(startStr);
                        int end = start + 1;

                        // Print the result to the output file
                        outputFile << id << "\t" << chromo << "\t" << start << "\t" << end << "\t" << probabilityBrdU << "\t" << strand << std::endl;
                    } else {
                        std::cerr << "Warning: Invalid START value in line: " << line << std::endl;
                    }
                }
            }
        }
    }

    outputFile.close();
}

int main() {
    std::string samples = std::getenv("SAMPLES"); // Assuming SAMPLES is an environment variable
    std::vector<std::string> sampleArray;

    // Splitting the exported string into an array
    split(samples, sampleArray);

    // Process each sample directory
    for (const auto& subDirectory : sampleArray) {
        for (const auto& entry : fs::directory_iterator("input/" + subDirectory)) {
            std::string fileName = entry.path().filename().string();
            processFile(subDirectory, fileName);
        }
    }

    // Create the fork probability dataset
    createForkProbabilityDataset(sampleArray);

    return 0;
}
