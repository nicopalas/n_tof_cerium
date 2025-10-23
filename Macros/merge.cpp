#include <TChain.h>
#include <TSystem.h>
#include <iostream>

void merge() {
    TChain chain("coincidences");
    int filesAdded = 0;
    
    for (int run = 118558; run <= 118794; run++) {
        TString fileName = Form("/Users/nico/Desktop/Tese/Results/anodes/coincidences_final/threshold_400/out_run%d.root", run);
        
        if (gSystem->AccessPathName(fileName)) {
            std::cout << "Missing: " << fileName << std::endl;
            continue;
        }
        
        // Try to add file and check if tree exists
        Long64_t entriesBefore = chain.GetEntries();
        chain.Add(fileName);
        Long64_t entriesAfter = chain.GetEntries();
        
        if (entriesAfter > entriesBefore) {
            std::cout << "Added: " << fileName << std::endl;
            filesAdded++;
        } else {
            std::cout << "No coincidences tree in: " << fileName << std::endl;
        }
    }
    
    if (filesAdded > 0) {
        std::cout << "\nMerging " << filesAdded << " files..." << std::endl;
        chain.Merge("/Users/nico/Desktop/Tese/Results/anodes/out_anodes_400.root", "fast");
        std::cout << "Done! Output: out_anodes_400.root" << std::endl;
    } else {
        std::cout << "No files to merge!" << std::endl;
    }
}