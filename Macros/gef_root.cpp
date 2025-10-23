#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

int main() {
    // Input file
    std::string inputFileName = "/Users/nico/Desktop/Tese/GEF/events_U238_100MeV.lmd";
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return 1;
    }

    // Output ROOT file
    std::string outputFileName = "fission_events_100MeV.root";
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree *tree = new TTree("FissionTree", "Tree storing fission event data");

    // Variables to store data
    Int_t Z1 = 0, Z2 = 0;
    Int_t A1post = 0, A2post = 0;
    Double_t tke = 0;
    Double_t tke1 = 0, tke2 = 0;

    

    // Branches
    tree->Branch("Z1", &Z1, "Z1/I");
    tree->Branch("Z2", &Z2, "Z2/I");
    tree->Branch("A1", &A1post, "A1/I");
    tree->Branch("A2", &A2post, "A2/I");
    tree->Branch("TKE", &tke, "TKE/D");
    tree->Branch("KE1", &tke1, "KE1/D");
    tree->Branch("KE2", &tke2, "KE2/D");

    // Read and parse the input file
    std::string line;
    std::vector<std::string> tokens;
    bool isEvent = false;

    while (std::getline(inputFile, line)) {
        if (line.empty() || line[0] == '*') {
            isEvent = false; // Reset event flag on empty or comment lines
            continue;
        }

        std::istringstream iss(line);
        std::string token;

        // Tokenize the line
        while (iss >> token) {
            tokens.push_back(token);
        }

        // Check if the line starts a new event
        if (!isEvent && tokens.size() > 2 && tokens[1] == "92" && tokens[2] == "238") {
            isEvent = true; // Start of a new event
        }

        // If the line ends but the data is incomplete, continue reading
        if (tokens.size() < 18) continue;

        // Extract data
        try {
            Z1 = std::stoi(tokens[4]);
            Z2 = std::stoi(tokens[5]);
            A1post = std::stoi(tokens[8]);
            A2post = std::stoi(tokens[9]);
            tke = std::stod(tokens[25]);
            tke1 = static_cast<double>(tke) / (1.0 + static_cast<double>(A1post) / static_cast<double>(A2post));
            tke2 = static_cast<double>(tke) / (1.0 + static_cast<double>(A2post) / static_cast<double>(A1post));
           

        } catch (const std::exception &e) {
            std::cerr << "Error parsing line: " << line << "\n" << e.what() << std::endl;
            tokens.clear();
            continue;
        }

        // Fill the tree
        tree->Fill();
        tokens.clear(); // Clear tokens for the next line
    }

    // Write the tree to the ROOT file
    tree->Write();
    outputFile->Close();

    // Close the input file
    inputFile.close();

    std::cout << "Data successfully saved to " << outputFileName << std::endl;
    return 0;
}