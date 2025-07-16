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
void neutron_energy(int run_number) {
    std::vector<Double_t> gamma_flash_means; // now we will compute the mean of the gamma flash
    std::string infile_name = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/gamma_flash/threshold=400/gamma_flash_" + std::to_string(run_number) + ".root";
    TFile *infile = TFile::Open(infile_name.c_str()); //input files
      // loading data from previous step
    TTree *intree = (TTree*)infile->Get("coincidences");

    intree->SetBranchAddress("gamma_flash0", gamma_flash0);
    intree->SetBranchAddress("gamma_flash1", gamma_flash1);
    intree->SetBranchAddress("gamma_flash2", gamma_flash2);
    intree->SetBranchAddress("gamma_flash3", gamma_flash3);
    intree->SetBranchAddress("gamma_flash4", gamma_flash4);
    intree->SetBranchAddress("gamma_flash5", gamma_flash5);
    intree->SetBranchAddress("gamma_flash6", gamma_flash6);
    intree->SetBranchAddress("gamma_flash7", gamma_flash7);
    intree->SetBranchAddress("gamma_flash8", gamma_flash8);
    intree->SetBranchAddress("gamma_flash9", gamma_flash9);
    intree->SetBranchAddress("detn", detn_all);

    Long64_t nentries = intree->GetEntries();
    if (nentries == 0) {
        std::cout << "gamma_flash empty" << std::endl;
        return;
    }
    std::vector<TH1D*> hists; //vector of hists to see everything works graphically
    for (int i = 0; i < 10; i++) {
        hists.push_back(new TH1D(Form("hist_%d", i), Form("Gamma Flash %d", i), 100, -740, -720)); // restrict to the range of interest of gamma flash
    }

    for (Long64_t j = 0; j < nentries; j++) {
        intree->GetEntry(j);
        // we make the histograms in order to perform a gaussian fit. we will use the mean of the gaussian fit as the gamma flash to avoid the tails.
        // maybe this is not the best approach but it is the one that I have thought of...
        hists[0]->Fill(gamma_flash0[0]);
        hists[1]->Fill(gamma_flash1[0]);
        hists[2]->Fill(gamma_flash2[0]);
        hists[3]->Fill(gamma_flash3[0]);
        hists[4]->Fill(gamma_flash4[0]);
        hists[5]->Fill(gamma_flash5[0]);
        hists[6]->Fill(gamma_flash6[0]);
        hists[7]->Fill(gamma_flash7[0]);
        hists[8]->Fill(gamma_flash8[0]);
        hists[9]->Fill(gamma_flash9[0]);
    }

    for (int i = 0; i < 10; i++) {
        TCanvas *c = new TCanvas(Form("c_%d", i), Form("Gamma Flash Histogram %d", i), 800, 600);
        hists[i]->Draw();
        TF1 *gausFit = new TF1("gausFit", "gaus", -740, -720);
        hists[i]->Fit(gausFit, "Q"); // fit the histogram with a gaussian
        gamma_flash_means.push_back(gausFit->GetParameter(1)); // Get the mean of the Gaussian fit

        TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend->AddEntry(hists[i], Form("Entries: %d", (int)hists[i]->GetEntries()), "l");
        legend->AddEntry(gausFit, Form("Mean: %.2f", gausFit->GetParameter(1)), "l");
        legend->AddEntry(gausFit, Form("Sigma: %.2f", gausFit->GetParameter(2)), "l");
        legend->AddEntry(gausFit, Form("Chi2/NDF: %.2f", gausFit->GetChisquare() / gausFit->GetNDF()), "l");
        legend->Draw();
       
    }

    infile->Close();

    // Printing
    std::cout <<"INFO Gamma Flash Means:" << std::endl;
    for (size_t i = 0; i < gamma_flash_means.size(); i++) {
        std::cout << "Mean of detector " << i << ": " << gamma_flash_means[i] << std::endl;
    }
    std::vector<Double_t> offsets_gamma_flash(10); //for synchronization we will use the last detector as reference (uranium) and we will calibrate the rest of the detectors.
    for (size_t i = 0; i < gamma_flash_means.size(); i++) {
        offsets_gamma_flash[i] = gamma_flash_means[i] - gamma_flash_means[0]; //offset
    }

    std::string outfile2 = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/anodes_final/threshold=400/out_run" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(outfile2.c_str(), "RECREATE"); //output file
    TTree *outtree = new TTree("coincidences", "Corrected Coincidences Tree");

    std::string file_name = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/coincidences_raw/Threshold=400/output_run" + std::to_string(run_number) + ".root";
    TFile *inputFile = TFile::Open(file_name.c_str()); //input files
    TTree *inputTree = (TTree*)inputFile->Get("nTOF_coincidences");  
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
    Double_t tof;

    pkupTree->SetBranchAddress("RunNumber", &pkup_RunNumber);
    pkupTree->SetBranchAddress("time", &pkup_time);
    pkupTree->SetBranchAddress("psTime", &pkup_psTime);
    pkupTree->SetBranchAddress("BunchNumber", &pkup_BunchNumber);
    pkupTree->SetBranchAddress("PSpulse", &pkup_PSpulse);
    pkupTree->SetBranchAddress("PulseIntensity", &pkup_PulseIntensity);
    pkupTree->SetBranchAddress("tof", &tof);
    pkupTree->SetBranchAddress("amp", &amp);
    inputTree->SetBranchAddress("RunNumber", &RunNumber);
    inputTree->SetBranchAddress("BunchNumber", &BunchNumber);
    inputTree->SetBranchAddress("PSpulse", &PSpulse);
    inputTree->SetBranchAddress("time", &eventTime);
    inputTree->SetBranchAddress("psTime", &psTime);
    inputTree->SetBranchAddress("PulseIntensity", &PulseIntensity);
    inputTree->SetBranchAddress("amp0", amp0);
    inputTree->SetBranchAddress("amp1", amp1);
    inputTree->SetBranchAddress("amp2", amp2);
    inputTree->SetBranchAddress("amp3", amp3);
    inputTree->SetBranchAddress("amp4", amp4);
    inputTree->SetBranchAddress("amp5", amp5);
    inputTree->SetBranchAddress("amp6", amp6);
    inputTree->SetBranchAddress("amp7", amp7);
    inputTree->SetBranchAddress("amp8", amp8);
    inputTree->SetBranchAddress("amp9", amp9);
    inputTree->SetBranchAddress("detn", &detn_all);
    inputTree->SetBranchAddress("mult0", &mult0);
    inputTree->SetBranchAddress("mult1", &mult1);
    inputTree->SetBranchAddress("mult2", &mult2);
    inputTree->SetBranchAddress("mult3", &mult3);
    inputTree->SetBranchAddress("mult4", &mult4);
    inputTree->SetBranchAddress("mult5", &mult5);
    inputTree->SetBranchAddress("mult6", &mult6);
    inputTree->SetBranchAddress("mult7", &mult7);
    inputTree->SetBranchAddress("mult8", &mult8);
    inputTree->SetBranchAddress("mult9", &mult9);
    inputTree->SetBranchAddress("mult", &mult);
    inputTree->SetBranchAddress("tof0", tof0);
    inputTree->SetBranchAddress("tof1", tof1);
    inputTree->SetBranchAddress("tof2", tof2);
    inputTree->SetBranchAddress("tof3", tof3);
    inputTree->SetBranchAddress("tof4", tof4);
    inputTree->SetBranchAddress("tof5", tof5);
    inputTree->SetBranchAddress("tof6", tof6);
    inputTree->SetBranchAddress("tof7", tof7);
    inputTree->SetBranchAddress("tof8", tof8);
    inputTree->SetBranchAddress("tof9", tof9);
    Double_t gamma_flash0[5], gamma_flash1[5], gamma_flash2[5], gamma_flash3[5], gamma_flash4[5];
    Double_t gamma_flash5[5], gamma_flash6[5], gamma_flash7[5], gamma_flash8[5], gamma_flash9[5];
    Double_t neutron_energy=0; //initially fixed to 0
    //the new tree will combine gamma_flash (constant and equal to the mean for every entry) with the rest of the data
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
    outtree->Branch("detn", &detn_all, "detn[10]/I");
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
    outtree->Branch("mult", &mult, "mult/I");
    outtree->Branch("tof0", tof0, "tof0[mult0]/D");
    outtree->Branch("tof1", tof1, "tof1[mult1]/D");
    outtree->Branch("tof2", tof2, "tof2[mult2]/D");
    outtree->Branch("tof3", tof3, "tof3[mult3]/D");
    outtree->Branch("tof4", tof4, "tof4[mult4]/D");
    outtree->Branch("tof5", tof5, "tof5[mult5]/D");
    outtree->Branch("tof6", tof6, "tof6[mult6]/D");
    outtree->Branch("tof7", tof7, "tof7[mult7]/D");
    outtree->Branch("tof8", tof8, "tof8[mult8]/D");
    outtree->Branch("tof9", tof9, "tof9[mult9]/D");
    //gamma_flash of the same size of tof in order to easily check the time calibration (just check tof-gamma_flash for mult>8 and see the mean is near 0)
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
    outtree-> Branch("neutron_energy", &neutron_energy, "neutron_energy/D");

    Long64_t kentries = inputTree->GetEntries();
    Long64_t npkupEntries = pkupTree->GetEntries();
    for (Long64_t i = 0; i < kentries; i++) {
        inputTree->GetEntry(i);
        for (Long64_t j = 0; j < npkupEntries; j++) {
            pkupTree->GetEntry(j); //again adjusting with the pickup data
            if (RunNumber == pkup_RunNumber && PulseIntensity == pkup_PulseIntensity && BunchNumber == pkup_BunchNumber && PSpulse == pkup_PSpulse && eventTime == pkup_time && psTime == pkup_psTime) {
                // calibrate the tof with the offset and the pickup times and store them.
               for (int j = 0; j < mult0; j++) { tof0[j] = tof0[j] - tof - offsets_gamma_flash[0]; gamma_flash0[j] = gamma_flash_means[0]; }
                for (int j = 0; j < mult1; j++) { tof1[j] = tof1[j] - tof - offsets_gamma_flash[1]+0.0704/299792458.0*1E9; gamma_flash1[j] = gamma_flash_means[0]+0.03/299792458.0*10E9; } //30 mm distance between PPACS
                for (int j = 0; j < mult2; j++) { tof2[j] = tof2[j] - tof - offsets_gamma_flash[2]+2*0.0704/299792458.0*1E9; gamma_flash2[j] = gamma_flash_means[0]+2*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult3; j++) { tof3[j] = tof3[j] - tof - offsets_gamma_flash[3]+3*0.0704/299792458.0*1E9; gamma_flash3[j] = gamma_flash_means[0]+3*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult4; j++) { tof4[j] = tof4[j] - tof - offsets_gamma_flash[4]+4*0.0704/299792458.0*1E9; gamma_flash4[j] = gamma_flash_means[0]+4*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult5; j++) { tof5[j] = tof5[j] - tof - offsets_gamma_flash[5]+5*0.0704/299792458.0*1E9; gamma_flash5[j] = gamma_flash_means[0]+5*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult6; j++) { tof6[j] = tof6[j] - tof - offsets_gamma_flash[6]+6*0.0704/299792458.0*1E9; gamma_flash6[j] = gamma_flash_means[0]+6*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult7; j++) { tof7[j] = tof7[j] - tof - offsets_gamma_flash[7]+7*0.0704/299792458.0*1E9; gamma_flash7[j] = gamma_flash_means[0]+7*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult8; j++) { tof8[j] = tof8[j] - tof - offsets_gamma_flash[8]+8*0.0704/299792458.0*1E9; gamma_flash8[j] = gamma_flash_means[0]+8*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult9; j++) { tof9[j] = tof9[j] - tof - offsets_gamma_flash[9]+9*0.0704/299792458.0*1E9; gamma_flash9[j] = gamma_flash_means[0]+9*(0.03/299792458.0*1E9); }
 //the key in finding the neutron energy is that it must be computed only if a real coincidence appears.
                // a fission trigger can only involve 2,3 or a maximum of 4 detectors (D.Tarrio thesis)
                //so what we will do is the following:
                //1. search for the first detector that triggers the coincidence (detn_all[j]==1)
                //2. compute the neutron energy in each case via E_n=m_n*c^2*(beta-1). add 185/0.3 to account for the distance travelled by photons
                //3. there may be some coincidences that are not real because they involve non-adjacent detectors.
                // For example, if the coincidence is between detectors 0 and 9, the neutron energy will be NaN because the tof difference is too small (random coincidence).
                // this was tested at first and i decided to set this type of coincidences to 0 energy since they are not relevant (at least their energy isn't)
                for (int j = 0; j < 9; j++) {
 // we only need to go to 9 because if detn_all[9] is the first one to trigger then the coincidence will not involve adjacent detectors, hence it is not real.
                    if (detn_all[j] == 1 && detn_all[j+1]==1) {
                        switch (j) {
                case 0:
                    neutron_energy = neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof0[0] - gamma_flash0[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 1:
                    neutron_energy =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof1[0] - gamma_flash1[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 2:
                    neutron_energy =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof2[0] - gamma_flash2[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 3:
                    neutron_energy =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof3[0] - gamma_flash3[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 4:
                    neutron_energy =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof4[0] - gamma_flash4[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 5:
                    neutron_energy = neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof5[0] - gamma_flash5[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 6:
                    neutron_energy =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof6[0] - gamma_flash6[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 7:
                    neutron_energy =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof7[0] - gamma_flash7[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 8:
                    neutron_energy =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof8[0] - gamma_flash8[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 9:
                    neutron_energy = neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof9[0] - gamma_flash9[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
            }
            break; // once you find the first anode signal break the loop
        }
                if (std::isnan(neutron_energy)) {
                    neutron_energy = 0; 
                //sometimes the energy is NaN because the tof difference is too small. This happens due to random coincidences betweeen non-adjacent detectors
                }
                }
            outtree->Fill(); //we only fill the hist if the if statement(the pickup coincides with the signal).)
            neutron_energy=0; //set the neutron energy again to 0.
                    }
                }
    }
std::cout << "INFO: Writing output file with" << outtree->GetEntries() << "entries" << std::endl;

outfile->Write();
outfile->Close();
file2->Close();
inputFile->Close();
}
