#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <TF1.h>
#include <TMath.h>
#include <vector>
// this first part is the same as in anode_analysis.cpp
Int_t mult=1;
Int_t mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
Int_t RunNumber, eventTime, BunchNumber, PSpulse;
Double_t psTime;
Double_t tof;
Float_t PulseIntensity;
Double_t tof0[20],tof1[20],tof2[20],tof3[20],tof4[20],tof5[20],tof6[20],tof7[20],tof8[20],tof9[20]; //must be hardcoded if arrays are defined at file scope
Float_t amp0[20],amp1[20],amp2[20],amp3[20],amp4[20],amp5[20],amp6[20],amp7[20],amp8[20],amp9[20];
Double_t new_tof0[20], new_tof1[20], new_tof2[20], new_tof3[20], new_tof4[20], new_tof5[20], new_tof6[20], new_tof7[20], new_tof8[20], new_tof9[20];
 //must be hardcoded if arrays are defined at file scope
Double_t gamma_flash0[20], gamma_flash1[20], gamma_flash2[20], gamma_flash3[20], gamma_flash4[20], gamma_flash5[20], gamma_flash6[20], gamma_flash7[20], gamma_flash8[20], gamma_flash9[20];
Int_t detn_all[10];
Double_t c=299792458.0; //speed of light in m/s
Double_t neutron_mass= 939.565420;// mev/c^2
Double_t l=185.0; // flight path in meters
// we load the variables that we need to compute the neutron energy

void gamma_flash(int run_number) {

    std::string filename = "/Users/nico/Desktop/Tese/Results/anodes/coincidences_raw/toy/output_run" + std::to_string(run_number) + ".root";
    TFile *file = TFile::Open(filename.c_str()) ;//we need raw data for the pickup
    TTree *intree = (TTree*)file->Get("nTOF_coincidences");


    intree->SetBranchStatus("*", 1); // turn on all branches
    intree->SetBranchAddress("RunNumber", &RunNumber);
    intree->SetBranchAddress("time", &eventTime);
    intree->SetBranchAddress("psTime", &psTime);
    intree->SetBranchAddress("BunchNumber", &BunchNumber);
    intree->SetBranchAddress("PSpulse", &PSpulse);
    intree->SetBranchAddress("PulseIntensity", &PulseIntensity);
    intree->SetBranchAddress("detn", detn_all);
    intree->SetBranchAddress("mult0", &mult0);
    intree->SetBranchAddress("mult1", &mult1);
    intree->SetBranchAddress("mult2", &mult2);
    intree->SetBranchAddress("mult3", &mult3);
    intree->SetBranchAddress("mult4", &mult4);
    intree->SetBranchAddress("mult5", &mult5);
    intree->SetBranchAddress("mult6", &mult6);
    intree->SetBranchAddress("mult7", &mult7);
    intree->SetBranchAddress("mult8", &mult8);
    intree->SetBranchAddress("mult9", &mult9);
    intree->SetBranchAddress("tof0", tof0);
    intree->SetBranchAddress("tof1", tof1);
    intree->SetBranchAddress("tof2", tof2);
    intree->SetBranchAddress("tof3", tof3);
    intree->SetBranchAddress("tof4", tof4);
    intree->SetBranchAddress("tof5", tof5);
    intree->SetBranchAddress("tof6", tof6);
    intree->SetBranchAddress("tof7", tof7);
    intree->SetBranchAddress("tof8", tof8);
    intree->SetBranchAddress("tof9", tof9);
    intree->SetBranchAddress("amp0", amp0);
    intree->SetBranchAddress("amp1", amp1);
    intree->SetBranchAddress("amp2", amp2);
    intree->SetBranchAddress("amp3", amp3);
    intree->SetBranchAddress("amp4", amp4);
    intree->SetBranchAddress("amp5", amp5);
    intree->SetBranchAddress("amp6", amp6);
    intree->SetBranchAddress("amp7", amp7);
    intree->SetBranchAddress("amp8", amp8);
    intree->SetBranchAddress("amp9", amp9);
    intree->SetBranchAddress("mult", &mult);

    Long64_t nentries = intree->GetEntries(); // getting the number of entries to loop over them
    std::string out_file = "/Users/nico/Desktop/Tese/Results/anodes/gamma_flash/gamma_flash_" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(out_file.c_str(), "RECREATE");
    TTree *outtree = new TTree("coincidences", "Coincidences Tree");
    outtree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    outtree->Branch("BunchNumber", &BunchNumber, "BunchNumber/I");
    outtree->Branch("PSpulse", &PSpulse, "PSpulse/I");
    outtree->Branch("time", &eventTime, "time/I");
    outtree->Branch("psTime", &psTime, "psTime/D");
    outtree->Branch("PulseIntensity", &PulseIntensity, "PulseIntensity/F");
    outtree->Branch("mult0", &mult0, "mult0/I");
    outtree->Branch("mult1", &mult1, "mult1/I");
    outtree->Branch("mult2", &mult2, "mult2/I");
    outtree->Branch("mult3", &mult3, "mult3/I");
    outtree->Branch("mult4", &mult4, "mult4/I");
    outtree->Branch("mult5", &mult5, "mult5/I");
    outtree->Branch("mult6", &mult6, "mult6/I");
    outtree->Branch("mult7", &mult7, "mult7/I");
    outtree->Branch("mult8", &mult8, "mult8/I");
    outtree->Branch("mult9", &mult9, "mult9/I");
    outtree->Branch("amp0", amp0, "amp0[mult0]/F");
    outtree->Branch("amp1", amp1, "amp1[mult1]/F");
    outtree->Branch("amp2", amp2, "amp2[mult2]/F");
    outtree->Branch("amp3", amp3, "amp3[mult3]/F");
    outtree->Branch("amp4", amp4, "amp4[mult4]/F");
    outtree->Branch("amp5", amp5, "amp5[mult5]/F");
    outtree->Branch("amp6", amp6, "amp6[mult6]/F");
    outtree->Branch("amp7", amp7, "amp7[mult7]/F");
    outtree->Branch("amp8", amp8, "amp8[mult8]/F");
    outtree->Branch("amp9", amp9, "amp9[mult9]/F");
    outtree->Branch("detn", detn_all, "detn[10]/I");
    outtree->Branch("mult", &mult, "mult/I");
    outtree->Branch("gamma_flash0", gamma_flash0, "gamma_flash0[mult0]/D");
    outtree->Branch("gamma_flash1", gamma_flash1, "gamma_flash1[mult1]/D");
    outtree->Branch("gamma_flash2", gamma_flash2, "gamma_flash2[mult2]/D");
    outtree->Branch("gamma_flash3", gamma_flash3, "gamma_flash3[mult3]/D");
    outtree->Branch("gamma_flash4", gamma_flash4, "gamma_flash4[mult4]/D");
    outtree->Branch("gamma_flash5", gamma_flash5, "gamma_flash5[mult5]/D");
    outtree->Branch("gamma_flash6", gamma_flash6, "gamma_flash6[mult6]/D");
    outtree->Branch("gamma_flash7", gamma_flash7, "gamma_flash7[mult7]/D");
    outtree->Branch("gamma_flash8", gamma_flash8, "gamma_flash8[mult8]/D");
    outtree->Branch("gamma_flash9", gamma_flash9, "gamma_flash9[mult9]/D");
    // now we will loop over the entries and search  the gamma flash, i.e. coincidence in all detectors
    for (Long64_t j = 0; j < nentries; j++) {
        intree->GetEntry(j);
                if (mult > 8 && mult0>0 && mult1>0 && mult2>0 && mult3>0 && mult4>0 && mult5>0 && mult6>0 && mult7>0 && mult8>0 && mult9>0) { //trigger in every detector
                    for(int k = 0; k < mult0; k++){gamma_flash0[k] = tof0[k];} //readjusting the gamma flash
                    for(int k = 0; k < mult1; k++){gamma_flash1[k] = tof1[k];}
                    for(int k = 0; k < mult2; k++){gamma_flash2[k] = tof2[k];}
                    for(int k = 0; k < mult3; k++){gamma_flash3[k] = tof3[k];}
                    for(int k = 0; k < mult4; k++){gamma_flash4[k] = tof4[k] ;}
                    for(int k = 0; k < mult5; k++){gamma_flash5[k] = tof5[k] ;}
                    for(int k = 0; k < mult6; k++){gamma_flash6[k] = tof6[k] ;}
                    for(int k = 0; k < mult7; k++){gamma_flash7[k] = tof7[k] ;}
                    for(int k = 0; k < mult8; k++){gamma_flash8[k] = tof8[k] ;}
                    for(int k = 0; k < mult9; k++){gamma_flash9[k] = tof9[k] ;}
                    outtree->Fill();
                    std::cout << "Filled entry " << outtree->GetEntries() << std::endl;

                }
            }
outtree->Write();
outfile->Close();
file->Close();
file2->Close();
std::cout << "Finished writing ROOT file!" << std::endl;

}
