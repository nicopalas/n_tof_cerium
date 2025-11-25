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
#include<TSpectrum.h>
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
Double_t neutron_energy_0, neutron_energy_1, neutron_energy_2, neutron_energy_3, neutron_energy_4, neutron_energy_5;
Double_t neutron_energy_6, neutron_energy_7, neutron_energy_8;
// we load the variables that we need to compute the neutron energy
void neutron_energy(int run_number) {
    std::vector<Double_t> gamma_flash_means; // now we will compute the mean of the gamma flash
    std::string infile_name = "gamma_flash_" + std::to_string(run_number) + ".root";
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
    intree->SetBranchAddress("detn_all", detn_all);
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

    Long64_t nentries = intree->GetEntries();
    if (nentries == 0) {
        std::cout << "[INFO] gamma_flash tree is empty." << std::endl;
        return;
    }

    std::cout << "[INFO] Processing gamma_flash tree with " << nentries << " entries." << std::endl;

    // Create histograms for each gamma_flash channel
    std::vector<TH1D*> hists;
    for (int i = 0; i < 10; i++) {
        hists.push_back(new TH1D(Form("h_flash%d", i), Form("Gamma Flash %d", i), 200, -760, -710));
        hists[i]->GetXaxis()->SetTitle("TOF (ns)");
        hists[i]->GetYaxis()->SetTitle("Counts");
    }

    // Fill histograms with all multiplicities for each channel
    for (Long64_t j = 0; j < nentries; j++) {
        intree->GetEntry(j);

        double *flashes[10] = {gamma_flash0, gamma_flash1, gamma_flash2, gamma_flash3, gamma_flash4,
                               gamma_flash5, gamma_flash6, gamma_flash7, gamma_flash8, gamma_flash9};
        int mults[10] = {mult0, mult1, mult2, mult3, mult4, mult5, mult6, mult7, mult8, mult9};

        for (int i = 0; i < 10; i++) {
            for (int k = 0; k < mults[i]; k++) {
                hists[i]->Fill(flashes[i][k]);
            }
        }
    }

    // Analyze histograms: find main gamma flash peak via TSpectrum
    for (int i = 0; i < 10; i++) {
        TCanvas *c = new TCanvas(Form("c_flash_%d", i), Form("Gamma Flash %d", i), 900, 700);
        hists[i]->Draw();

        TSpectrum spectrum(10);  // up to 10 possible peaks
        Int_t nPeaks = spectrum.Search(hists[i], 1, "nodraw", 0.4); // σ=1 bins smoothing, threshold 40%

        double *peaks = spectrum.GetPositionX();
        double bestPeakX = 0;
        double maxY = 0;

        for (int p = 0; p < nPeaks; ++p) {
            double valY = hists[i]->GetBinContent(hists[i]->FindBin(peaks[p]));
            if (valY > maxY) {
                maxY = valY;
                bestPeakX = peaks[p];
            }
        }

        if (nPeaks > 0) {
            gamma_flash_means.push_back(bestPeakX);
            std::cout << "[INFO] Channel " << i << ": found peak at TOF = " << bestPeakX
                      << " (" << nPeaks << " candidates)\n";

            // Draw marker line
            TLine *line = new TLine(bestPeakX, 0, bestPeakX, hists[i]->GetMaximum());
            line->SetLineColor(kRed);
            line->SetLineWidth(2);
            line->Draw("same");

            TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
            leg->AddEntry(hists[i], Form("Entries: %d", (int)hists[i]->GetEntries()), "l");
            leg->AddEntry((TObject*)0, Form("Peak: %.2f ns", bestPeakX), "");
            leg->Draw();
        } else {
            gamma_flash_means.push_back(0.0);
            std::cout << "[WARN] Channel " << i << ": no peak found.\n";
        }
    }

    std::cout << "\n==== Gamma Flash Peak Summary ====\n";
    for (size_t i = 0; i < gamma_flash_means.size(); ++i) {
        std::cout << "Channel " << i << " → peak = " << gamma_flash_means[i] << " ns\n";
    }



    std::vector<Double_t> offsets_gamma_flash(10); //for synchronization we will use the last detector as reference (uranium) and we will calibrate the rest of the detectors.
    for (size_t i = 0; i < gamma_flash_means.size(); i++) {
        offsets_gamma_flash[i] = gamma_flash_means[i] - gamma_flash_means[0]; //offset
    }

    std::string outfile2 = "out_run" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(outfile2.c_str(), "RECREATE"); //output file
    TTree *outtree = new TTree("coincidences", "Corrected Coincidences Tree");

    std::string file_name = "output_run" + std::to_string(run_number) + ".root";
    TFile *inputFile = TFile::Open(file_name.c_str()); //input files
    TTree *inputTree = (TTree*)inputFile->Get("nTOF_coincidences");  


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
    inputTree->SetBranchAddress("detn_all", detn_all);
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
    Double_t gamma_flash0[20], gamma_flash1[20], gamma_flash2[20], gamma_flash3[20], gamma_flash4[20];
    Double_t gamma_flash5[20], gamma_flash6[20], gamma_flash7[20], gamma_flash8[20], gamma_flash9[20];
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
    outtree->Branch("detn_all", detn_all, "detn[10]/I");
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
    outtree-> Branch("neutron_energy_0", &neutron_energy_0, "neutron_energy_0/D");
    outtree-> Branch("neutron_energy_1", &neutron_energy_1, "neutron_energy_1/D");
    outtree-> Branch("neutron_energy_2", &neutron_energy_2, "neutron_energy_2/D");
    outtree-> Branch("neutron_energy_3", &neutron_energy_3, "neutron_energy_3/D");
    outtree-> Branch("neutron_energy_4", &neutron_energy_4, "neutron_energy_4/D");
    outtree-> Branch("neutron_energy_5", &neutron_energy_5, "neutron_energy_5/D");
    outtree-> Branch("neutron_energy_6", &neutron_energy_6, "neutron_energy_6/D");
    outtree-> Branch("neutron_energy_7", &neutron_energy_7, "neutron_energy_7/D");
    outtree-> Branch("neutron_energy_8", &neutron_energy_8, "neutron_energy_8/D");

    Long64_t kentries = inputTree->GetEntries();
    for (Long64_t i = 0; i < kentries; i++) {
        inputTree->GetEntry(i);
                // 50 mm distance between PPACS (we use 0.0704 because it is the distance between the center of the detectors having into account the inclination))
                for (int j = 0; j < mult0; j++) { tof0[j] = tof0[j] - offsets_gamma_flash[0]; gamma_flash0[j] = gamma_flash_means[0]; }
                for (int j = 0; j < mult1; j++) { tof1[j] = tof1[j] - offsets_gamma_flash[1]+0.0704/299792458.0*1E9; gamma_flash1[j] = gamma_flash_means[0]+0.03/299792458.0*10E9; }
                for (int j = 0; j < mult2; j++) { tof2[j] = tof2[j] - offsets_gamma_flash[2]+2*0.0704/299792458.0*1E9; gamma_flash2[j] = gamma_flash_means[0]+2*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult3; j++) { tof3[j] = tof3[j] - offsets_gamma_flash[3]+3*0.0704/299792458.0*1E9; gamma_flash3[j] = gamma_flash_means[0]+3*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult4; j++) { tof4[j] = tof4[j] - offsets_gamma_flash[4]+4*0.0704/299792458.0*1E9; gamma_flash4[j] = gamma_flash_means[0]+4*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult5; j++) { tof5[j] = tof5[j] - offsets_gamma_flash[5]+5*0.0704/299792458.0*1E9; gamma_flash5[j] = gamma_flash_means[0]+5*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult6; j++) { tof6[j] = tof6[j] - offsets_gamma_flash[6]+6*0.0704/299792458.0*1E9; gamma_flash6[j] = gamma_flash_means[0]+6*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult7; j++) { tof7[j] = tof7[j] - offsets_gamma_flash[7]+7*0.0704/299792458.0*1E9; gamma_flash7[j] = gamma_flash_means[0]+7*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult8; j++) { tof8[j] = tof8[j] - offsets_gamma_flash[8]+8*0.0704/299792458.0*1E9; gamma_flash8[j] = gamma_flash_means[0]+8*(0.03/299792458.0*1E9); }
                for (int j = 0; j < mult9; j++) { tof9[j] = tof9[j] - offsets_gamma_flash[9]+9*0.0704/299792458.0*1E9; gamma_flash9[j] = gamma_flash_means[0]+9*(0.03/299792458.0*1E9); }
 //the key in finding the neutron energy is that it must be computed only if a real coincidence appears.
                // a fission trigger can only involve 2,3 or a maximum of 4 detectors (D.Tarrio thesis)
                //so what we will do is the following:
                //1. search for the first detector that triggers the coincidence (detn_all[j]==1)
                //2. compute the neutron energy in each case via E_n=m_n*c^2*(beta-1). add 185/0.3 to account for the distance travelled by photons
                //3. there may be some coincidences that are not real because they involve non-adjacent detectors.
                // For example, if the coincidence is between detectors 0 and 9, the neutron energy will be NaN because the tof difference is too small (random coincidence).
                // this was tested at first and i decided to set this type of coincidences to 0 energy since they are not relevant (at least their energy isn't)
                for (int j = 0; j <= 8; j++) {
 // we only need to go to 8 because if detn_all[8] is the first one to trigger then the coincidence will not involve adjacent detectors, hence it is not real.
                    if (detn_all[j] >= 1 && detn_all[j+1] >= 1) {
                        switch (j) {
                case 0:
                    neutron_energy_0 = neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof0[0] - gamma_flash0[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 1:
                    neutron_energy_1 =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof1[0] - gamma_flash1[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 2:
                    neutron_energy_2 =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof2[0] - gamma_flash2[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 3:
                    neutron_energy_3 =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof3[0] - gamma_flash3[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 4:
                    neutron_energy_4 =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof4[0] - gamma_flash4[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 5:
                    neutron_energy_5 = neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof5[0] - gamma_flash5[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 6:
                    neutron_energy_6 =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof6[0] - gamma_flash6[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 7:
                    neutron_energy_7 =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof7[0] - gamma_flash7[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
                case 8:
                    neutron_energy_8 =  neutron_mass * (1/ TMath::Power(1-TMath::Power(l / (c*(tof8[0] - gamma_flash8[0]+185/0.3)*1E-9), 2),0.5)-1);
                    break;
            }
            break; // once you find the first anode signal break the loop
        }
                    
            //sometimes the energy is NaN because the tof difference is too small. This happens due to random coincidences betweeen non-adjacent detectors
        if (TMath::IsNaN(neutron_energy_0) || neutron_energy_0>5000.) neutron_energy_0 = -1000.;
        if (TMath::IsNaN(neutron_energy_1) || neutron_energy_1>5000.) neutron_energy_1 = -1000.;
        if (TMath::IsNaN(neutron_energy_2) || neutron_energy_2>5000.) neutron_energy_2 = -1000.;
        if (TMath::IsNaN(neutron_energy_3) || neutron_energy_3>5000.) neutron_energy_3 = -1000.;
        if (TMath::IsNaN(neutron_energy_4) || neutron_energy_4>5000.) neutron_energy_4 = -1000.;
        if (TMath::IsNaN(neutron_energy_5) || neutron_energy_5>5000.) neutron_energy_5 = -1000.;
        if (TMath::IsNaN(neutron_energy_6) || neutron_energy_6>5000.) neutron_energy_6 = -1000.;
        if (TMath::IsNaN(neutron_energy_7) || neutron_energy_7>5000.) neutron_energy_7 = -1000.;
        if (TMath::IsNaN(neutron_energy_8) || neutron_energy_8>5000.) neutron_energy_8 = -1000.;
    }

    outtree->Fill(); //we only fill the hist if the if statement(the pickup coincides with the signal).)
    neutron_energy_0 = 0; neutron_energy_1 = 0; neutron_energy_2 = 0; neutron_energy_3 = 0; neutron_energy_4 = 0; neutron_energy_5 = 0; neutron_energy_6 = 0; neutron_energy_7 = 0; neutron_energy_8 = 0; //reset energies for the next entry
                    }
    std::cout << "INFO: Writing output file with" << outtree->GetEntries() << "entries" << std::endl;

    outfile->Write();
    outfile->Close();
    inputFile->Close();
}
