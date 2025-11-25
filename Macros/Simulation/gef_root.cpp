#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cctype>
#include "TFile.h"
#include "TTree.h"

int main() {
    // Input file
    std::string inputFileName = "/Users/nico/Desktop/Tese/Macros/Macros/n_tof_cerium/GEF_data/Z92_A238_n_E3MeV.lmd";
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return 1;
    }

    // Output ROOT file
    std::string outputFileName = "/Users/nico/Desktop/Tese/Macros/Macros/n_tof_cerium/ROOT_files/fission_events_uranium_lowenergy.root";
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree *tree = new TTree("FissionTree", "Tree storing fission event data");

    // Variables to store data
    Double_t neutronEnergy = 0;
    Int_t Z1 = 0, Z2 = 0;
    Int_t A1post = 0, A2post = 0;
    Double_t TKEpost = 0;
    Double_t KE_1 = 0, KE_2 = 0;

    // Branches
    tree->Branch("NeutronEnergy", &neutronEnergy, "NeutronEnergy/D");
    tree->Branch("Z1", &Z1, "Z1/I");
    tree->Branch("Z2", &Z2, "Z2/I");
    tree->Branch("A1", &A1post, "A1/I");
    tree->Branch("A2", &A2post, "A2/I");
    tree->Branch("TKEpost", &TKEpost, "TKEpost/D");
    tree->Branch("KE_1", &KE_1, "KE_1/D");
    tree->Branch("KE_2", &KE_2, "KE_2/D");

    // Read and parse the input file
    std::string line;
    int eventCount = 0;
    Double_t currentNeutronEnergy = 0;
    bool foundNeutronEnergy = false;

    while (std::getline(inputFile, line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Check if this is a neutron energy header line (even if it's a comment)
        if (line.find("formed by (n,f) with En =") != std::string::npos) {
            // Extract neutron energy from the line
            size_t enPos = line.find("En =");
            if (enPos != std::string::npos) {
                // Start after "En ="
                size_t startPos = enPos + 4;
                
                // Extract the number - it's after spaces
                std::string energyStr;
                for (size_t i = startPos; i < line.length(); i++) {
                    if (std::isdigit(line[i]) || line[i] == '.') {
                        // Found start of number
                        while (i < line.length() && (std::isdigit(line[i]) || line[i] == '.')) {
                            energyStr += line[i];
                            i++;
                        }
                        break;
                    }
                }
                
                if (!energyStr.empty()) {
                    try {
                        currentNeutronEnergy = std::stod(energyStr);
                        foundNeutronEnergy = true;
                        std::cout << "Found neutron energy: " << currentNeutronEnergy << " MeV" << std::endl;
                    } catch (const std::exception &e) {
                        std::cerr << "Error parsing neutron energy from: '" << energyStr << "'" << std::endl;
                    }
                }
            }
        }
        // Skip comment lines that don't contain neutron energy
        else if (line[0] == '*') {
            continue;
        }
        // Process data lines (lines that contain numbers and are not comments)
        else {
            // Skip lines that are too short
            if (line.length() < 10) continue;
            
            // Check if line contains digits (data lines start with spaces followed by digits)
            bool hasDigit = false;
            for (char c : line) {
                if (std::isdigit(c)) {
                    hasDigit = true;
                    break;
                }
            }
            
            if (hasDigit) {
                std::istringstream iss(line);
                std::vector<std::string> tokens;
                std::string token;
                
                // Tokenize the line
                while (iss >> token) {
                    tokens.push_back(token);
                }

                // Check if we have enough tokens (at least 26 for TKEpost)
                if (tokens.size() >= 26) {
                    try {
                        // Extract the required variables
                        // According to the description:
                        // tokens[4] = Z1, tokens[5] = Z2
                        // tokens[8] = A1post, tokens[9] = A2post
                        // tokens[25] = TKEpost
                        Z1 = std::stoi(tokens[4]);
                        Z2 = std::stoi(tokens[5]);
                        A1post = std::stoi(tokens[8]);
                        A2post = std::stoi(tokens[9]);
                        TKEpost = std::stod(tokens[25]);
                        KE_1 = TKEpost * A2post / (A1post + A2post);
                        KE_2 = TKEpost * A1post / (A1post + A2post);

                        // Use the current neutron energy for this event
                        neutronEnergy = currentNeutronEnergy;

                        // Fill the tree
                        tree->Fill();
                        eventCount++;

                        // Debug: print first few events
                        if (eventCount <= 5) {
                            std::cout << "Event " << eventCount << ": Z1=" << Z1 << ", Z2=" << Z2 
                                      << ", A1=" << A1post << ", A2=" << A2post 
                                      << ", TKEpost=" << TKEpost 
                                      << ", En=" << neutronEnergy << " MeV" << std::endl;
                        }

                    } catch (const std::exception &e) {
                        std::cerr << "Error parsing event data from line: " << line << std::endl;
                        std::cerr << "Exception: " << e.what() << std::endl;
                        continue;
                    }
                } else {
                    // Only show error for lines that look like they should be data
                    if (tokens.size() > 10) {
                        std::cerr << "Insufficient tokens in line: " << line << std::endl;
                        std::cerr << "Expected at least 26 tokens, found " << tokens.size() << std::endl;
                        std::cerr << "First few tokens: ";
                        for (size_t i = 0; i < tokens.size() && i < 5; ++i) {
                            std::cerr << "'" << tokens[i] << "' ";
                        }
                        std::cerr << std::endl;
                    }
                }
            }
        }
    }

    if (!foundNeutronEnergy) {
        std::cerr << "Warning: Could not find neutron energy in file" << std::endl;
    }

    // Write the tree to the ROOT file
    tree->Write();
    outputFile->Close();

    // Close the input file
    inputFile.close();

    std::cout << "Successfully processed " << eventCount << " events" << std::endl;
    std::cout << "Data saved to " << outputFileName << std::endl;

    return 0;
}

  