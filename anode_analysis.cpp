
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <vector>
#include <tuple>
#include <iostream>
#include <cstdlib>
Int_t mult=1; //mult will be the number of signals in coincidence (each detector may have multiple signals)
Int_t max_mul = 5; //max number of signals in coincidence in a single detector to limit the size
Int_t mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
Double_t tof0[5],tof1[5],tof2[5],tof3[5],tof4[5],tof5[5],tof6[5],tof7[5],tof8[5],tof9[5]; //must be hardcoded if arrays are defined at file scope
Float_t amp0[5],amp1[5],amp2[5],amp3[5],amp4[5],amp5[5],amp6[5],amp7[5],amp8[5],amp9[5];  //must be hardcoded if arrays are defined at file scope
Int_t detn_all[10]; 
 // we have two variables that represent the same, but one is an array (with the signals in coincidence in each detector) 
 //and the other is int but both count the number of signals in each coincidence
Bool_t addCoincidences(Int_t the_detn) {
    if(the_detn<0 || the_detn>9) {        
        std::cout << "ERROR in addTof: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(the_detn==0) { mult0++; detn_all[the_detn]++; }
    if(the_detn==1) { mult1++; detn_all[the_detn]++; }
    if(the_detn==2) { mult2++; detn_all[the_detn]++; }
    if(the_detn==3) { mult3++; detn_all[the_detn]++; }
    if(the_detn==4) { mult4++; detn_all[the_detn]++; }
    if(the_detn==5) { mult5++; detn_all[the_detn]++; }
    if(the_detn==6) { mult6++; detn_all[the_detn]++; }
    if(the_detn==7) { mult7++; detn_all[the_detn]++; }
    if(the_detn==8) { mult8++; detn_all[the_detn]++; }
    if(the_detn==9) { mult9++; detn_all[the_detn]++; }
    return true;
}

Bool_t addTof(Int_t the_detn, Double_t the_tof) {
    if(the_detn<0 || the_detn>9) {        
        std::cout << "ERROR in addTof: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(detn_all[the_detn] < 1 || detn_all[the_detn] > max_mul) {         //check the mutliplicity is valid
        std::cout << "ERROR in addTof: detn_all[" << the_detn <<"] value: " << detn_all[the_detn] << " not valid" << std::endl;
        for (int i = 0; i < detn_all[the_detn]; ++i) {
            if (the_detn == 0) std::cout << "tof0[" << i << "] = " << tof0[i] << std::endl;
            if (the_detn == 1) std::cout << "tof1[" << i << "] = " << tof1[i] << std::endl;
            if (the_detn == 2) std::cout << "tof2[" << i << "] = " << tof2[i] << std::endl;
            if (the_detn == 3) std::cout << "tof3[" << i << "] = " << tof3[i] << std::endl;
            if (the_detn == 4) std::cout << "tof4[" << i << "] = " << tof4[i] << std::endl;
            if (the_detn == 5) std::cout << "tof5[" << i << "] = " << tof5[i] << std::endl;
            if (the_detn == 6) std::cout << "tof6[" << i << "] = " << tof6[i] << std::endl;
            if (the_detn == 7) std::cout << "tof7[" << i << "] = " << tof7[i] << std::endl;
            if (the_detn == 8) std::cout << "tof8[" << i << "] = " << tof8[i] << std::endl;
            if (the_detn == 9) std::cout << "tof9[" << i << "] = " << tof9[i] << std::endl;
        }
        return false;
    }
    if(the_detn==0) tof0[mult0 - 1] = the_tof;
    if(the_detn==1) tof1[mult1 - 1] = the_tof;
    if(the_detn==2) tof2[mult2 - 1] = the_tof;
    if(the_detn==3) tof3[mult3 - 1] = the_tof;
    if(the_detn==4) tof4[mult4 - 1] = the_tof;
    if(the_detn==5) tof5[mult5 - 1] = the_tof;
    if(the_detn==6) tof6[mult6 - 1] = the_tof;
    if(the_detn==7) tof7[mult7 - 1] = the_tof;
    if(the_detn==8) tof8[mult8 - 1] = the_tof;
    if(the_detn==9) tof9[mult9 - 1] = the_tof;
    return true;
}

Bool_t addAmp(Int_t the_detn, Float_t the_amp) {
    if(the_detn<0 || the_detn>9) {        
        std::cout << "ERROR in addAmp: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(detn_all[the_detn] < 1 || detn_all[the_detn] > max_mul) {        
        std::cout << "ERROR in addAmp: detn_all[" << the_detn <<"] value: " << detn_all[the_detn] << " not valid" << std::endl;
        return false;
    }
    if(the_detn==0) amp0[mult0 - 1] = the_amp;
    if(the_detn==1) amp1[mult1 - 1] = the_amp;
    if(the_detn==2) amp2[mult2 - 1] = the_amp;
    if(the_detn==3) amp3[mult3 - 1] = the_amp;
    if(the_detn==4) amp4[mult4 - 1] = the_amp;
    if(the_detn==5) amp5[mult5 - 1] = the_amp;
    if(the_detn==6) amp6[mult6 - 1] = the_amp;
    if(the_detn==7) amp7[mult7 - 1] = the_amp;
    if(the_detn==8) amp8[mult8 - 1] = the_amp;
    if(the_detn==9) amp9[mult9 - 1] = the_amp;
    return true;
}
void anode_analysis(int run_number, float threshold, double time_for_coincidence) {

    std::cout << "INFO: Processing run_number " << run_number << std::endl;
    std::vector<TFile*> files;
    std::vector<TTree*> trees;
    for (int run = run_number; run <= run_number; ++run) { //select runs in case of multiple runs at once, but usually will be just one.
        std::string filename = "/nucl_lustre/n_tof_INTC_P_665/DATA/run" + std::to_string(run) + ".root";
        TFile *file = TFile::Open(filename.c_str());
        if (file && !file->IsZombie()) {
            files.push_back(file);
            TTree *tree;
            file->GetObject("PPAN", tree);
            
            if (tree) {
                trees.push_back(tree);
            }
        }
    }
    Long64_t nentries = 0, nbytes = 0, nb = 0;
    Long64_t jentry = 0, ientry = 0;
    std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>> signals;
    std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>> signals_pickup;
    std::vector<std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>>> configurations;

    Double_t psTime, tof;
    Int_t RunNumber, time, BunchNumber, PSpulse, detn;
    Float_t PulseIntensity, amp;

    TTree *tree = nullptr;
    for (auto tree : trees) {
        tree->SetBranchStatus("*", 0); // turn off all branches
        tree->SetBranchStatus("RunNumber", 1); // Int_t
        tree->SetBranchStatus("time", 1); // Int_t
        tree->SetBranchStatus("psTime", 1); // Double_t
        tree->SetBranchStatus("BunchNumber", 1); // Int_t
        tree->SetBranchStatus("PSpulse", 1); // Int_t
        tree->SetBranchStatus("PulseIntensity", 1); // Float_t
        tree->SetBranchStatus("detn", 1); // Int_t
        tree->SetBranchStatus("tof", 1); // Double_t
        tree->SetBranchStatus("amp", 1); // Float_t
        tree->SetBranchAddress("RunNumber", &RunNumber);
        tree->SetBranchAddress("time", &time);
        tree->SetBranchAddress("psTime", &psTime);
        tree->SetBranchAddress("BunchNumber", &BunchNumber);
        tree->SetBranchAddress("PSpulse", &PSpulse);
        tree->SetBranchAddress("PulseIntensity", &PulseIntensity);
        tree->SetBranchAddress("detn", &detn);
        tree->SetBranchAddress("tof", &tof);
        tree->SetBranchAddress("amp", &amp);
        
        std::cout << "INFO: new tree found and configured" << std::endl;
        nentries = tree->GetEntriesFast();
        std::cout << "INFO: tree entries: " << nentries << std::endl;
        for (jentry = 0; jentry < nentries; jentry++) {
            ientry = tree->LoadTree(jentry);
            if (ientry < 0) break;
            nb = tree->GetEntry(jentry); nbytes += nb;
            if (detn <7){
                tof+=15; 
                // adding offset to the time in the first detectors (taking uranium as reference)
                //the reason for this offset is the different time references between the 2 DAQ systems used.
                // we add this offset so coincidences can be found between all detectors within a time window.
                // because if they were to be found without this offset, the time window would have to be very large
                // and this would lead to a lot of accidental coincidences.
            }
            signals.emplace_back(RunNumber, time, psTime, BunchNumber, PSpulse, PulseIntensity, detn, tof, amp); //insert tuple at the end
        }
    }
        
        std::cout <<"INFO: number of signals: "<< signals.size() << std::endl; // print number of signals
        std::vector<bool> used(signals.size(), false); // track used signals to avoid double counting
    
        for (size_t i = 0; i < signals.size(); ++i) {
            if (used[i] == false && std::get<8>(signals[i]) > threshold) { 
                std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>> coincidences;
                 // start in the first detector and not used

                for (size_t j = i + 1; j < signals.size(); ++j) { // search in the next signals in next detectors
                    if (std::get<8>(signals[j])>threshold  && used[j] == false){
                        if (std::get<0>(signals[j]) == std::get<0>(signals[i]) && std::get<5>(signals[j]) == std::get<5>(signals[i]) && std::get<4>(signals[j]) == std::get<4>(signals[i])) { //check everything coincides
                            if (std::get<3>(signals[j]) == std::get<3>(signals[i]) && std::get<2>(signals[j]) == std::get<2>(signals[i]) && std::get<1>(signals[j]) == std::get<1>(signals[i])) {
                                if (abs(std::get<7>(signals[j]) - std::get<7>(signals[i])) < time_for_coincidence) {
                                    used[j] = true;
                                    used[i] = true;
                                    //mark as used only after the coincidence is found
                                    if (coincidences.size() == 0) {
                                        coincidences.push_back(signals[i]); //add only once the first signal
                                    
                                    }
                                    coincidences.push_back(signals[j]);
                                }
                            }
                        }
                    }
                }
                if (!coincidences.empty()) {
                    configurations.push_back(coincidences); 
                                                            // save coincidences into a vector
                }
                
                
            }
        }
        std::cout << "INFO: number of coincidences: " << configurations.size() << std::endl;
        //now we want to search for those so-called lonely coincidences. we have to loop over the signals again. this time, we will treat differently those who
        //have been used and those who have not. if they were used, we will search over the other coincidence. if they coincide we would merge them.
        //if they were not used, we will search for the closest coincidence and merge them...
        Int_t lonelyCoincidence = 0;
        for (size_t i = 0; i < signals.size(); ++i) {
            if (used[i] == false && std::get<8>(signals[i]) > threshold) {
                for (size_t j=0; j<configurations.size();j++){
                    for (size_t k = 0; k < configurations[j].size(); ++k) {
                        if (std::get<0>(signals[i]) == std::get<0>(configurations[j][k]) && std::get<8>(configurations[j][k]) >threshold && std::get<5>(signals[i]) == std::get<5>(configurations[j][k]) && std::get<4>(signals[i]) == std::get<4>(configurations[j][k])) {
                            if (std::get<1>(signals[i])==std::get<1>(configurations[j][k]) && std::get<3>(signals[i]) == std::get<3>(configurations[j][k]) && std::get<4>(signals[i]) == std::get<4>(configurations[j][k]) && std::get<2>(signals[i]) == std::get<2>(configurations[j][k])) {
                                if (abs(std::get<7>(signals[i]) - std::get<7>(configurations[j][k])) < time_for_coincidence) {
                                    lonelyCoincidence++;
                                    break;
    
                                    configurations[j].push_back(signals[i]);
                                    used[i] = true;
                                }
                            }
                        }
                    }
                }
                }
            if (used[i]==true && std::get<8>(signals[i])>threshold) {
                Int_t coincidence_index=-1;
                for (Int_t j = 0; j < configurations.size(); ++j) {
                    for (Int_t k = 0; k < configurations[j].size(); ++k) {
                        if (std::get<0>(signals[i]) == std::get<0>(configurations[j][k]) && std::get<1>(configurations[j][k]) == std::get<1>(signals[i]) && std::get<2>(configurations[j][k]) == std::get<2>(signals[i])) {
                            if (std::get<3>(signals[i]) == std::get<3>(configurations[j][k]) && std::get<4>(signals[i]) == std::get<4>(configurations[j][k]) && std::get<5>(signals[i]) == std::get<5>(configurations[j][k])) {
                                if (std::get<6>(signals[i]) == std::get<6>(configurations[j][k]) && std::get<7>(signals[i]) == std::get<7>(configurations[j][k]) && std::get<8>(signals[i]) == std::get<8>(configurations[j][k])) {
                                    coincidence_index = j;
                                    break;
                                    }
                                }
                            }
                        }
                        if (coincidence_index != -1) { //if you find the index stop the loop
                            break;
                        }
                    }
                    for (size_t k=0; k<configurations.size();k++){ //loop to find coincidences
                        if (k!= coincidence_index){
                        for (size_t l=0; l<configurations[k].size(); l++){
                            if (std::get<8>(configurations[k][l])>threshold && std::get<0>(signals[i]) == std::get<0>(configurations[k][l]) && std::get<1>(configurations[k][l]) == std::get<1>(signals[i]) && std::get<2>(configurations[k][l]) == std::get<2>(signals[i])) {
                                if (std::get<3>(signals[i]) == std::get<3>(configurations[k][l]) && std::get<4>(signals[i]) == std::get<4>(configurations[k][l]) && std::get<5>(signals[i]) == std::get<5>(configurations[k][l]) && std::get<6>(signals[i])== std::get<6>(configurations[k][l])) {
                                        if (abs(std::get<7>(signals[i])- std::get<7>(configurations[k][l])) < time_for_coincidence) {
                                                for (size_t m=0; m<configurations[k].size(); m++){
                                                    configurations[coincidence_index].push_back(configurations[k][m]); //add the signals from the other coincidence to the first one
                                                }
                                            configurations.erase(configurations.begin()+k);
                                            break; //remove the element merged to avoid double counting
                                        }
                                    }
                                }
                            }
                            }
                        }
                    }
                }
    std::cout << "INFO: number of coincidences (after merging): " << configurations.size() << std::endl;
    std::cout << "INFO: number of lonely coincidences: " << lonelyCoincidence << std::endl; //printing to ensure everything is working
    for(Int_t index=0; index<10; index++) detn_all[index] = 0; //initialize
    // Output file
    std::string output_filename = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/coincidences_raw/Threshold=400/output_run" + std::to_string(run_number) + ".root";
    TFile* outputFile = TFile::Open(output_filename.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "ERROR: Could not create output file " << output_filename << std::endl;
    }
    TTree *outtree = new TTree("nTOF_coincidences","nTOF_coincidences");
    outtree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    outtree->Branch("time",&time,"time/I");
    outtree->Branch("psTime",&psTime,"psTime/D");
    outtree->Branch("BunchNumber", &BunchNumber, "BunchNumber/I");
    outtree->Branch("PSpulse", &PSpulse, "PSpulse/I");
    outtree->Branch("PulseIntensity", &PulseIntensity, "PulseIntensity/F");
    outtree->Branch("detn",detn_all,"detn_all[10]/I");
    outtree->Branch("mult0",&mult0,"mult0/I");
    outtree->Branch("mult1",&mult1,"mult1/I");
    outtree->Branch("mult2",&mult2,"mult2/I");
    outtree->Branch("mult3",&mult3,"mult3/I");
    outtree->Branch("mult4",&mult4,"mult4/I");
    outtree->Branch("mult5",&mult5,"mult5/I");
    outtree->Branch("mult6",&mult6,"mult6/I");
    outtree->Branch("mult7",&mult7,"mult7/I");
    outtree->Branch("mult8",&mult8,"mult8/I");
    outtree->Branch("mult9",&mult9,"mult9/I");
    outtree->Branch("tof0",tof0,"tof0[mult0]/D");
    outtree->Branch("tof1",tof1,"tof1[mult1]/D");
    outtree->Branch("tof2",tof2,"tof2[mult2]/D");
    outtree->Branch("tof3",tof3,"tof3[mult3]/D");
    outtree->Branch("tof4",tof4,"tof4[mult4]/D");
    outtree->Branch("tof5",tof5,"tof5[mult5]/D");
    outtree->Branch("tof6",tof6,"tof6[mult6]/D");
    outtree->Branch("tof7",tof7,"tof7[mult7]/D");
    outtree->Branch("tof8",tof8,"tof8[mult8]/D");
    outtree->Branch("tof9",tof9,"tof9[mult9]/D");
    outtree->Branch("amp0",amp0,"amp0[mult0]/F");
    outtree->Branch("amp1",amp1,"amp1[mult1]/F");
    outtree->Branch("amp2",amp2,"amp2[mult2]/F");
    outtree->Branch("amp3",amp3,"amp3[mult3]/F");
    outtree->Branch("amp4",amp4,"amp4[mult4]/F");
    outtree->Branch("amp5",amp5,"amp5[mult5]/F");
    outtree->Branch("amp6",amp6,"amp6[mult6]/F");
    outtree->Branch("amp7",amp7,"amp7[mult7]/F");
    outtree->Branch("amp8",amp8,"amp8[mult8]/F");
    outtree->Branch("amp9",amp9,"amp9[mult9]/F");
    outtree->Branch("mult",&mult,"mult/I");
    

    for (size_t i = 0; i < configurations.size(); ++i) {
            RunNumber = std::get<0>(configurations[i][0]);
            time = std::get<1>(configurations[i][0]);
            psTime = std::get<2>(configurations[i][0]);
            BunchNumber = std::get<3>(configurations[i][0]);
            PSpulse = std::get<4>(configurations[i][0]);
            PulseIntensity = std::get<5>(configurations[i][0]); //store every coincidence in the tree
        for (size_t j = 0; j < configurations[i].size(); ++j) {//loop over configurations and over each signal in the configurations
            
            detn = std::get<6>(configurations[i][j]);
            tof = std::get<7>(configurations[i][j]);
            amp = std::get<8>(configurations[i][j]);
            addCoincidences(detn);
            addTof(detn, tof);
            addAmp(detn, amp);
        }
        mult = configurations[i].size(); //mult will be the number of signals in coincidence

        outtree->Fill(); //fill the tree
    
        mult=1;
        mult0 = 0; mult1 = 0; mult2 = 0; mult3 = 0; mult4 = 0; mult5 = 0; mult6 = 0; mult7 = 0; mult8 = 0; mult9 = 0;
    
        for(Int_t index=0; index<10; index++) detn_all[index] = 0;
    }
        outtree->Write();
        outputFile->Close();
    }
