#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

// Função para processar cada linha de entrada e aplicar as transformações
void processKeepValues(const std::string &line, std::ofstream &leftBedgraph, std::ofstream &rightBedgraph) {
    std::stringstream ss(line);
    std::string col1, col2, col3, col4, col5, col6;
    ss >> col1 >> col2 >> col3 >> col4 >> col5 >> col6;

    // Converte colunas 5 e 6 para double
    double probLeft = std::stod(col5);
    double probRight = std::stod(col6);

    // Muda valores < 0.5 para 0
    probLeft = (probLeft < 0.5 ? 0 : probLeft);
    probRight = (probRight < 0.5 ? 0 : probRight);

    // Escreve nos arquivos bedgraph
    leftBedgraph << col2 << "\t" << col3 << "\t" << col4 << "\t" << std::fixed << std::setprecision(6) << probLeft << "\n";
    rightBedgraph << col2 << "\t" << col3 << "\t" << col4 << "\t" << std::fixed << std::setprecision(6) << probRight << "\n";
}

// Função para processar valores de mudança
void processChangeValues(const std::string &line, std::ofstream &leftBedgraph, std::ofstream &rightBedgraph) {
    std::stringstream ss(line);
    std::string col1, col2, col3, col4, col5, col6;
    ss >> col1 >> col2 >> col3 >> col4 >> col5 >> col6;

    // Converte colunas 5 e 6 para double
    double probLeft = std::stod(col5);
    double probRight = std::stod(col6);

    // Muda valores < 0.5 para 0 e >= 0.5 para 1
    probLeft = (probLeft >= 0.5 ? 1 : 0);
    probRight = (probRight >= 0.5 ? 1 : 0);

    // Escreve nos arquivos bedgraph
    leftBedgraph << col2 << "\t" << col3 << "\t" << col4 << "\t" << std::fixed << std::setprecision(6) << probLeft << "\n";
    rightBedgraph << col2 << "\t" << col3 << "\t" << col4 << "\t" << std::fixed << std::setprecision(6) << probRight << "\n";
}

// Função para criar os arquivos bedgraph
void createBedgraphFiles(const std::string &inputFile) {
    std::ifstream inFile(inputFile);
    if (!inFile) {
        std::cerr << "Erro ao abrir o arquivo de entrada: " << inputFile << std::endl;
        return;
    }

    // Cria arquivos para KeepValues
    std::ofstream leftKeepBedgraph("output/probabilityForkSense-KeepValues-left.bedgraph");
    std::ofstream rightKeepBedgraph("output/probabilityForkSense-KeepValues-Rigth.bedgraph");

    // Processa e gera arquivos KeepValues
    std::string line;
    while (std::getline(inFile, line)) {
        processKeepValues(line, leftKeepBedgraph, rightKeepBedgraph);
    }

    inFile.clear(); // Limpa o estado do fluxo
    inFile.seekg(0); // Volta ao início do arquivo

    // Cria arquivos para ChangeValues
    std::ofstream leftChangeBedgraph("output/probabilityForkSense-ChangeValues-Left.bedgraph");
    std::ofstream rightChangeBedgraph("output/probabilityForkSense-ChangeValues-Rigth.bedgraph");

    // Processa e gera arquivos ChangeValues
    while (std::getline(inFile, line)) {
        processChangeValues(line, leftChangeBedgraph, rightChangeBedgraph);
    }

    inFile.close();
    leftKeepBedgraph.close();
    rightKeepBedgraph.close();
    leftChangeBedgraph.close();
    rightChangeBedgraph.close();
}

int main() {
    const std::string inputFile = "input/probabilityForkSense.csv";
    createBedgraphFiles(inputFile);
    return 0;
}
