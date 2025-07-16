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
Int_t mult=1;
Int_t mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
Int_t RunNumber, eventTime, BunchNumber, PSpulse;
Double_t psTime;
Double_t tof;
Float_t PulseIntensity;
Double_t tof0[5],tof1[5],tof2[5],tof3[5],tof4[5],tof5[5],tof6[5],tof7[5],tof8[5],tof9[5]; //must be hardcoded if arrays are defined at file scope
Float_t amp0[5],amp1[5],amp2[5],amp3[5],amp4[5],amp5[5],amp6[5],amp7[5],amp8[5],amp9[5];
Double_t new_tof0[5], new_tof1[5], new_tof2[5], new_tof3[5], new_tof4[5], new_tof5[5], new_tof6[5], new_tof7[5], new_tof8[5], new_tof9[5];
 //must be hardcoded if arrays are defined at file scope
Double_t gamma_flash0[5], gamma_flash1[5], gamma_flash2[5], gamma_flash3[5], gamma_flash4[5], gamma_flash5[5], gamma_flash6[5], gamma_flash7[5], gamma_flash8[5], gamma_flash9[5];
Int_t detn_all[10];
Double_t c=299792458.0; //speed of light in m/s
Double_t neutron_mass= 939.565420;// mev/c^2
Double_t l=185.0; // flight path in meters
// we load the variables that we need to compute the neutron energy

void gamma_flash(int run_number) {

    std::string filename = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/coincidences_raw/Threshold=400/output_run" + std::to_string(run_number) + ".root";
    TFile *file = TFile::Open(filename.c_str()) ;//we need raw data for the pickup
    TTree *intree = (TTree*)file->Get("nTOF_coincidences");

    // Open the raw run118560.root file and load the PKUP tree
    std::string filename2 = "/nucl_lustre/n_tof_INTC_P_665/DATA/run" + std::to_string(run_number) + ".root";
    TFile *file2 = TFile::Open(filename2.c_str());//we need raw data for the pickup
    TTree *pkupTree = (TTree*)file2->Get("PKUP");
    pkupTree->SetBranchStatus("*", 0); // turn off all branches
    pkupTree->SetBranchStatus("RunNumber", 1); // Int_t
    pkupTree->SetBranchStatus("time", 1); // Int_t
    pkupTree->SetBranchStatus("psTime", 1); // Double_t
    pkupTree->SetBranchStatus("BunchNumber", 1); // Int_t
    pkupTree->SetBranchStatus("PSpulse", 1); // Int_t
    pkupTree->SetBranchStatus("PulseIntensity", 1); // Float_t
    pkupTree->SetBranchStatus("detn", 1); // Int_t
    pkupTree->SetBranchStatus("tof", 1); // Double_t
    pkupTree->SetBranchStatus("amp", 1); // Float_t

    Int_t pkup_RunNumber, pkup_time, pkup_BunchNumber, pkup_PSpulse;
    Double_t pkup_psTime;
    Float_t pkup_PulseIntensity;
    Float_t amp;

    pkupTree->SetBranchAddress("RunNumber", &pkup_RunNumber);
    pkupTree->SetBranchAddress("time", &pkup_time);
    pkupTree->SetBranchAddress("psTime", &pkup_psTime);
    pkupTree->SetBranchAddress("BunchNumber", &pkup_BunchNumber);
    pkupTree->SetBranchAddress("PSpulse", &pkup_PSpulse);
    pkupTree->SetBranchAddress("PulseIntensity", &pkup_PulseIntensity);
    pkupTree->SetBranchAddress("tof", &tof);
    pkupTree->SetBranchAddress("amp", &amp);

    intree->SetBranchStatus("*", 1); // turn on all branches
    intree->SetBranchAddress("RunNumber", &RunNumber);
    intree->SetBranchAddress("time", &eventTime);
    intree->SetBranchAddress("psTime", &psTime);
    intree->SetBranchAddress("BunchNumber", &BunchNumber);
    intree->SetBranchAddress("PSpulse", &PSpulse);
    intree->SetBranchAddress("PulseIntensity", &PulseIntensity);
    intree->SetBranchAddress("detn", &detn_all);
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
    intree->SetBranchAddress("tof0", &tof0);
    intree->SetBranchAddress("tof1", &tof1);
    intree->SetBranchAddress("tof2", &tof2);
    intree->SetBranchAddress("tof3", &tof3);
    intree->SetBranchAddress("tof4", &tof4);
    intree->SetBranchAddress("tof5", &tof5);
    intree->SetBranchAddress("tof6", &tof6);
    intree->SetBranchAddress("tof7", &tof7);
    intree->SetBranchAddress("tof8", &tof8);
    intree->SetBranchAddress("tof9", &tof9);
    intree->SetBranchAddress("amp0", &amp0);
    intree->SetBranchAddress("amp1", &amp1);
    intree->SetBranchAddress("amp2", &amp2);
    intree->SetBranchAddress("amp3", &amp3);
    intree->SetBranchAddress("amp4", &amp4);
    intree->SetBranchAddress("amp5", &amp5);
    intree->SetBranchAddress("amp6", &amp6);
    intree->SetBranchAddress("amp7", &amp7);
    intree->SetBranchAddress("amp8", &amp8);
    intree->SetBranchAddress("amp9", &amp9);
    intree->SetBranchAddress("mult", &mult);

    Long64_t nentries = intree->GetEntries(); // getting the number of entries to loop over them
    Long64_t npkupEntries = pkupTree->GetEntries();
    std::string out_file = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/gamma_flash/threshold=400/gamma_flash_" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(out_file.c_str(), "RECREATE");
    TTree *outtree = new TTree("coincidences", "Coincidences Tree");
    outtree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    outtree->Branch("BunchNumber", &BunchNumber, "BunchNumber/I");
    outtree->Branch("PSpulse", &PSpulse, "PSpulse/I");
    outtree->Branch("time", &eventTime, "time/I");
    outtree->Branch("psTime", &psTime, "psTime/D");
    outtree->Branch("PulseIntensity", &PulseIntensity, "PulseIntensity/F");
    outtree->Branch("amp0", &amp0, "amp0/F");
    outtree->Branch("amp1", &amp1, "amp1/F");
    outtree->Branch("amp2", &amp2, "amp2/F");
    outtree->Branch("amp3", &amp3, "amp3/F");
    outtree->Branch("amp4", &amp4, "amp4/F");
    outtree->Branch("amp5", &amp5, "amp5/F");
    outtree->Branch("amp6", &amp6, "amp6/F");
    outtree->Branch("amp7", &amp7, "amp7/F");
    outtree->Branch("amp8", &amp8, "amp8/F");
    outtree->Branch("amp9", &amp9, "amp9/F");
    outtree->Branch("detn", detn_all, "detn[10]/I");
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
    outtree->Branch("mult", &mult, "mult/I");
    outtree->Branch("gamma_flash0", &gamma_flash0, "gamma_flash0[mult0]/D");
    outtree->Branch("gamma_flash1", &gamma_flash1, "gamma_flash1[mult1]/D");
    outtree->Branch("gamma_flash2", &gamma_flash2, "gamma_flash2[mult2]/D");
    outtree->Branch("gamma_flash3", &gamma_flash3, "gamma_flash3[mult3]/D");
    outtree->Branch("gamma_flash4", &gamma_flash4, "gamma_flash4[mult4]/D");
    outtree->Branch("gamma_flash5", &gamma_flash5, "gamma_flash5[mult5]/D");
    outtree->Branch("gamma_flash6", &gamma_flash6, "gamma_flash6[mult6]/D");
    outtree->Branch("gamma_flash7", &gamma_flash7, "gamma_flash7[mult7]/D");
    outtree->Branch("gamma_flash8", &gamma_flash8, "gamma_flash8[mult8]/D");
    outtree->Branch("gamma_flash9", &gamma_flash9, "gamma_flash9[mult9]/D");
    // now we will loop over the entries and search  the gamma flash, i.e. coincidence in all detectors
    //we loop over pickup data and check if parameters coincide with the gamma_flash for a proper adjustment of the times
    //we dont have to worry about tof or amp since there is only one tof and amp per bunchNumber so there will never be ambivalence
    for (Long64_t i = 0; i < npkupEntries; i++) {
        pkupTree->GetEntry(i);
        for (Long64_t j = 0; j < nentries; j++) {
            intree->GetEntry(j);
                if (RunNumber == pkup_RunNumber && PulseIntensity == pkup_PulseIntensity && BunchNumber == pkup_BunchNumber && PSpulse == pkup_PSpulse && eventTime == pkup_time && psTime == pkup_psTime) { // they must coincide in order to compute gamma flash
                    if (mult > 8 && mult0>0 && mult1>0 && mult2>0 && mult3>0 && mult4>0 && mult5>0 && mult6>0 && mult7>0 && mult8>0 && mult9>0) { //trigger in every detector
                        for(int k = 0; k < mult0; k++){gamma_flash0[k] = tof0[k] - tof;} //readjusting the gamma flash
                        for(int k = 0; k < mult1; k++){gamma_flash1[k] = tof1[k] - tof;}
                        for(int k = 0; k < mult2; k++){gamma_flash2[k] = tof2[k] - tof;}
                        for(int k = 0; k < mult3; k++){gamma_flash3[k] = tof3[k] - tof;}
                        for(int k = 0; k < mult4; k++){gamma_flash4[k] = tof4[k] - tof;}
                        for(int k = 0; k < mult5; k++){gamma_flash5[k] = tof5[k] - tof;}
                        for(int k = 0; k < mult6; k++){gamma_flash6[k] = tof6[k] - tof;}
                        for(int k = 0; k < mult7; k++){gamma_flash7[k] = tof7[k] - tof;}
                        for(int k = 0; k < mult8; k++){gamma_flash8[k] = tof8[k] - tof;}
                        for(int k = 0; k < mult9; k++){gamma_flash9[k] = tof9[k] - tof;}
                        outtree->Fill();      // fill only if the conditions are met to avoid unnecessary data
			std::cout << outtree->GetEntries() << std::endl;

                        }
                }
            }
            }
outtree->Write();
outfile->Close();
file->Close();
file2->Close();
std::cout << "Finished writing ROOT file!" << std::endl;

}
