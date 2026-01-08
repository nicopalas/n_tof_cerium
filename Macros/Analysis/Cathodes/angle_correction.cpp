#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <algorithm>
#include <cmath>
#include <string>


int main(){
    TFile *fin = TFile::Open("final_selection.root", "READ");
    if (!fin || fin->IsZombie()){
        std::cerr << "Error opening the file \n";
        return 1;
    }
    TTree *tin =(TTree*)fin->Get("coincidences");
    if (!tin){
        std::cerr << "Tree coincidences not found in file" << "\n";
        fin->Close();
        return 2;
    }
    Int_t nentries = tin->GetEntries();
    std::cout << "Number of entries = " << nentries << "\n";


    double tofA7[10], tofA8[10], tofA9[10];
    int multA7, multA8, multA9;

    // Positions detector
    double x7[10], y7[10], x8[10], y8[10], x9[10], y9[10];

    // Amplitudes Anodes
    float ampA7[10], ampA8[10], ampA9[10];

    // Amplitudes Cathodes
    float ampC71[10], ampC72[10], ampC73[10], ampC74[10];
    float ampC81[10], ampC82[10], ampC83[10], ampC84[10];
    float ampC91[10], ampC92[10], ampC93[10], ampC94[10];

    // Sums
    double sumX7[10], sumY7[10], sumX8[10], sumY8[10], sumX9[10], sumY9[10];

    // Differences
    double diffX7[10], diffY7[10], diffX8[10], diffY8[10], diffX9[10], diffY9[10];

    // Ratios
    float ratioX7[10], ratioY7[10], ratioX8[10], ratioY8[10], ratioX9[10], ratioY9[10];

    // Sum of amplitudes
    float sumX7_ampA7[10], sumY7_ampA7[10], sumX8_ampA8[10], sumY8_ampA8[10], sumX9_ampA9[10], sumY9_ampA9[10];

    tin->SetBranchAddress("tofA7", tofA7);
    tin->SetBranchAddress("tofA8", tofA8);
    tin->SetBranchAddress("tofA9", tofA9);

    tin->SetBranchAddress("x7", x7);
    tin->SetBranchAddress("y7", y7);
    tin->SetBranchAddress("x8", x8);
    tin->SetBranchAddress("y8", y8);
    tin->SetBranchAddress("x9", x9);
    tin->SetBranchAddress("y9", y9);
    tin->SetBranchAddress("ampA7", ampA7);
    tin->SetBranchAddress("ampA8", ampA8);
    tin->SetBranchAddress("ampA9", ampA9);

    tin->SetBranchAddress("ampC71", ampC71);
    tin->SetBranchAddress("ampC72", ampC72);
    tin->SetBranchAddress("ampC73", ampC73);
    tin->SetBranchAddress("ampC74", ampC74);

    tin->SetBranchAddress("ampC81", ampC81);
    tin->SetBranchAddress("ampC82", ampC82);
    tin->SetBranchAddress("ampC83", ampC83);
    tin->SetBranchAddress("ampC84", ampC84);

    tin->SetBranchAddress("ampC91", ampC91);
    tin->SetBranchAddress("ampC92", ampC92);
    tin->SetBranchAddress("ampC93", ampC93);
    tin->SetBranchAddress("ampC94", ampC94);

    tin->SetBranchAddress("sumX7", sumX7);
    tin->SetBranchAddress("sumY7", sumY7);
    tin->SetBranchAddress("sumX8", sumX8);
    tin->SetBranchAddress("sumY8", sumY8);
    tin->SetBranchAddress("sumX9", sumX9);
    tin->SetBranchAddress("sumY9", sumY9);

    tin->SetBranchAddress("diffX7", diffX7);
    tin->SetBranchAddress("diffY7", diffY7);
    tin->SetBranchAddress("diffX8", diffX8);
    tin->SetBranchAddress("diffY8", diffY8);
    tin->SetBranchAddress("diffX9", diffX9);
    tin->SetBranchAddress("diffY9", diffY9);

    tin->SetBranchAddress("ratioX7", ratioX7);
    tin->SetBranchAddress("ratioY7", ratioY7);
    tin->SetBranchAddress("ratioX8", ratioX8);
    tin->SetBranchAddress("ratioY8", ratioY8);
    tin->SetBranchAddress("ratioX9", ratioX9);
    tin->SetBranchAddress("ratioY9", ratioY9);
    tin->SetBranchAddress("sumX7_ampA7", sumX7_ampA7);
    tin->SetBranchAddress("sumY7_ampA7", sumY7_ampA7);
    tin->SetBranchAddress("sumX8_ampA8", sumX8_ampA8);
    tin->SetBranchAddress("sumY8_ampA8", sumY8_ampA8);
    tin->SetBranchAddress("sumX9_ampA9", sumX9_ampA9);
    tin->SetBranchAddress("sumY9_ampA9", sumY9_ampA9);

    tin->SetBranchAddress("mult7", &multA7);
    tin->SetBranchAddress("mult8", &multA8);
    tin->SetBranchAddress("mult9", &multA9);


    TFile *file_uranium = TFile::Open("events_uranium.root", "RECREATE");
    TTree *tout_uranium = new TTree("events_uranium", "Selected events for uranium target");
    TFile *file_gold = TFile::Open("events_gold.root", "RECREATE");
    TTree *tout_gold = new TTree("events_gold", "Selected events for gold target");
    // Variables to be stored in the output trees
    double neutron_energy_gold;
    double neutron_energy_uranium;
    double cos_theta_gold;
    double cos_theta_uranium;
    double phi_gold;
    double phi_uranium;
    double tof7_gold, tof8_gold;
    double tof8_uranium, tof9_uranium;
    float ampA7_gold, ampA8_gold;
    float ampA8_uranium, ampA9_uranium;
    double x7_gold, y7_gold, x8_gold, y8_gold;
    double x8_uranium, y8_uranium, x9_uranium, y9_uranium;
    float ampC71_gold, ampC72_gold, ampC73_gold, ampC74_gold;
    float ampC81_gold, ampC82_gold, ampC83_gold, ampC84_gold;
    float ampC81_uranium, ampC82_uranium, ampC83_uranium, ampC84_uranium;
    float ampC91_uranium, ampC92_uranium, ampC93_uranium, ampC94_uranium;
    double sumX7_gold, sumY7_gold, sumX8_gold, sumY8_gold;
    double sumX8_uranium, sumY8_uranium, sumX9_uranium, sumY9_uranium;
    double diffX7_gold, diffY7_gold, diffX8_gold, diffY8_gold;
    double diffX8_uranium, diffY8_uranium, diffX9_uranium, diffY9_uranium;
    float ratioX7_gold, ratioY7_gold, ratioX8_gold, ratioY8_gold;
    float ratioX8_uranium, ratioY8_uranium, ratioX9_uranium, ratioY9_uranium;
    float sumX7_ampA7_gold, sumY7_ampA7_gold, sumX8_ampA8_gold, sumY8_ampA8_gold;
    float sumX8_ampA8_uranium, sumY8_ampA8_uranium, sumX9_ampA9_uranium, sumY9_ampA9_uranium;
    float multA7_gold, multA8_gold;
    float multA8_uranium, multA9_uranium;

    tout_gold->Branch("neutron_energy", &neutron_energy_gold, "neutron_energy/D");
    tout_gold->Branch("tof7", &tof7_gold, "tof7/D");
    tout_gold->Branch("tof8", &tof8_gold, "tof8/D");
    tout_gold->Branch("ampA7", &ampA7_gold, "ampA7/F");
    tout_gold->Branch("ampA8", &ampA8_gold, "ampA8/F");
    tout_gold->Branch("x7", &x7_gold, "x7/D");
    tout_gold->Branch("y7", &y7_gold, "y7/D");
    tout_gold->Branch("x8", &x8_gold, "x8/D");
    tout_gold->Branch("y8", &y8_gold, "y8/D");
    tout_gold->Branch("ampC71", &ampC71_gold, "ampC71/F");
    tout_gold->Branch("ampC72", &ampC72_gold, "ampC72/F");
    tout_gold->Branch("ampC73", &ampC73_gold, "ampC73/F");
    tout_gold->Branch("ampC74", &ampC74_gold, "ampC74/F");
    tout_gold->Branch("ampC81", &ampC81_gold, "ampC81/F");
    tout_gold->Branch("ampC82", &ampC82_gold, "ampC82/F");
    tout_gold->Branch("ampC83", &ampC83_gold, "ampC83/F");
    tout_gold->Branch("ampC84", &ampC84_gold, "ampC84/F");
    tout_gold->Branch("sumX7", &sumX7_gold, "sumX7/D");
    tout_gold->Branch("sumY7", &sumY7_gold, "sumY7/D");
    tout_gold->Branch("sumX8", &sumX8_gold, "sumX8/D");
    tout_gold->Branch("sumY8", &sumY8_gold, "sumY8/D");
    tout_gold->Branch("diffX7", &diffX7_gold, "diffX7/D");
    tout_gold->Branch("diffY7", &diffY7_gold, "diffY7/D");
    tout_gold->Branch("diffX8", &diffX8_gold, "diffX8/D");
    tout_gold->Branch("diffY8", &diffY8_gold, "diffY8/D");
    tout_gold->Branch("ratioX7", &ratioX7_gold, "ratioX7/F");
    tout_gold->Branch("ratioY7", &ratioY7_gold, "ratioY7/F");
    tout_gold->Branch("ratioX8", &ratioX8_gold, "ratioX8/F");
    tout_gold->Branch("ratioY8", &ratioY8_gold, "ratioY8/F");
    tout_gold->Branch("sumX7_ampA7", &sumX7_ampA7_gold, "sumX7_ampA7/F");
    tout_gold->Branch("sumY7_ampA7", &sumY7_ampA7_gold, "sumY7_ampA7/F");
    tout_gold->Branch("sumX8_ampA8", &sumX8_ampA8_gold, "sumX8_ampA8/F");
    tout_gold->Branch("sumY8_ampA8", &sumY8_ampA8_gold, "sumY8_ampA8/F");

    tout_uranium->Branch("neutron_energy", &neutron_energy_uranium, "neutron_energy/D");
    tout_uranium->Branch("tof8", &tof8_uranium, "tof8/D");
    tout_uranium->Branch("tof9", &tof9_uranium, "tof9/D");
    tout_uranium->Branch("ampA8", &ampA8_uranium, "ampA8/F");
    tout_uranium->Branch("ampA9", &ampA9_uranium, "ampA9/F");
    tout_uranium->Branch("x8", &x8_uranium, "x8/D");
    tout_uranium->Branch("y8", &y8_uranium, "y8/D");
    tout_uranium->Branch("x9", &x9_uranium, "x9/D");
    tout_uranium->Branch("y9", &y9_uranium, "y9/D");
    tout_uranium->Branch("ampC81", &ampC81_uranium, "ampC81/F");
    tout_uranium->Branch("ampC82", &ampC82_uranium, "ampC82/F");
    tout_uranium->Branch("ampC83", &ampC83_uranium, "ampC83/F");
    tout_uranium->Branch("ampC84", &ampC84_uranium, "ampC84/F");
    tout_uranium->Branch("ampC91", &ampC91_uranium, "ampC91/F");
    tout_uranium->Branch("ampC92", &ampC92_uranium, "ampC92/F");
    tout_uranium->Branch("ampC93", &ampC93_uranium, "ampC93/F");
    tout_uranium->Branch("ampC94", &ampC94_uranium, "ampC94/F");
    tout_uranium->Branch("sumX8", &sumX8_uranium, "sumX8/D");
    tout_uranium->Branch("sumY8", &sumY8_uranium, "sumY8/D");
    tout_uranium->Branch("sumX9", &sumX9_uranium, "sumX9/D");
    tout_uranium->Branch("sumY9", &sumY9_uranium, "sumY9/D");
    tout_uranium->Branch("diffX8", &diffX8_uranium, "diffX8/D");
    tout_uranium->Branch("diffY8", &diffY8_uranium, "diffY8/D");
    tout_uranium->Branch("diffX9", &diffX9_uranium, "diffX9/D");
    tout_uranium->Branch("diffY9", &diffY9_uranium, "diffY9/D");
    tout_uranium->Branch("ratioX8", &ratioX8_uranium, "ratioX8/F");
    tout_uranium->Branch("ratioY8", &ratioY8_uranium, "ratioY8/F");
    tout_uranium->Branch("ratioX9", &ratioX9_uranium, "ratioX9/F");
    tout_uranium->Branch("ratioY9", &ratioY9_uranium, "ratioY9/F");
    tout_uranium->Branch("sumX8_ampA8", &sumX8_ampA8_uranium, "sumX8_ampA8/F");
    tout_uranium->Branch("sumY8_ampA8", &sumY8_ampA8_uranium, "sumY8_ampA8/F");
    tout_uranium->Branch("sumX9_ampA9", &sumX9_ampA9_uranium, "sumX9_ampA9/F");
    tout_uranium->Branch("sumY9_ampA9", &sumY9_ampA9_uranium, "sumY9_ampA9/F");
    Int_t counts_uranium = 0;
    Int_t counts_gold = 0;
    for (Int_t i= 0 ; i<nentries; i++){
        tin->GetEntry(i);
        if (i % 100000 == 1) {
            std::cout << "Counts gold: " << counts_gold << " Counts uranium: " << counts_uranium << std::endl;
        }
        for (int j = 0; j<multA7; j++){
            for (int k = 0; k<multA8;k++){
                neutron_energy_gold = 939.56542 *(1/TMath::Power(1 - TMath::Power(184.5/(299792458.0*(tofA7[j] + 726.616 + 185/0.3)*1E-9),2),0.5) - 1);
                if (multA7>0 && multA8>0 && ampA7[j]>10000 && ampA8[k]>2000 && tofA7[j]-tofA8[k]<3  && abs(tofA7[j]-tofA8[k])<10 && ampA7[j]+ampA8[k]>16000 && tofA7[j]>-650){ // Apply amplitude cut  
                    // Fill variables for the output tree
                    tof7_gold = tofA7[j];
                    tof8_gold = tofA8[k];
                    ampA7_gold = ampA7[j];
                    ampA8_gold = ampA8[k];
                    x7_gold = x7[j];
                    y7_gold = y7[j];
                    x8_gold = x8[k];
                    y8_gold = y8[k];
                    ampC71_gold = ampC71[j];
                    ampC72_gold = ampC72[j];
                    ampC73_gold = ampC73[j];
                    ampC74_gold = ampC74[j];
                    ampC81_gold = ampC81[k];
                    ampC82_gold = ampC82[k];
                    ampC83_gold = ampC83[k];
                    ampC84_gold = ampC84[k];
                    sumX7_gold = sumX7[j];
                    sumY7_gold = sumY7[j];
                    sumX8_gold = sumX8[k];
                    sumY8_gold = sumY8[k];
                    diffX7_gold = diffX7[j];
                    diffY7_gold = diffY7[j];
                    diffX8_gold = diffX8[k];
                    diffY8_gold = diffY8[k];
                    ratioX7_gold = ratioX7[j];
                    ratioY7_gold = ratioY7[j];
                    ratioX8_gold = ratioX8[k];
                    ratioY8_gold = ratioY8[k];
                    sumX7_ampA7_gold = sumX7_ampA7[j];
                    sumY7_ampA7_gold = sumY7_ampA7[j];
                    sumX8_ampA8_gold = sumX8_ampA8[k];
                    sumY8_ampA8_gold = sumY8_ampA8[k];
                    tout_gold->Fill();
                    counts_gold++;
                }
            }
        }

        for (int j = 0 ; j<multA8 ; j++){
            for (int k = 0; k<multA9; k++){
                neutron_energy_uranium = 939.56542 *(1/TMath::Power(1 - TMath::Power(184.5/(299792458.0*(tofA8[j] + 726.417 + 185/0.3)*1E-9),2),0.5) - 1);
                if (multA8>0 && multA9>0 && abs(tofA9[k]-tofA8[j])<10 && ampA9[k]>2000 && ampA8[j]>4000 && ampA9[k]+ampA8[j]>9500 && tofA8[j]>-650){
                    tof8_uranium = tofA8[j];
                    tof9_uranium = tofA9[k];
                    ampA8_uranium = ampA8[j];
                    ampA9_uranium = ampA9[k];
                    x8_uranium = x8[j];
                    y8_uranium = y8[j];
                    x9_uranium = x9[k];
                    y9_uranium = y9[k];
                    ampC81_uranium = ampC81[j];
                    ampC82_uranium = ampC82[j];
                    ampC83_uranium = ampC83[j];
                    ampC84_uranium = ampC84[j];
                    ampC91_uranium = ampC91[k];
                    ampC92_uranium = ampC92[k];
                    ampC93_uranium = ampC93[k];
                    ampC94_uranium = ampC94[k];
                    sumX8_uranium = sumX8[j];
                    sumY8_uranium = sumY8[j];
                    sumX9_uranium = sumX9[k];
                    sumY9_uranium = sumY9[k];
                    diffX8_uranium = diffX8[j];
                    diffY8_uranium = diffY8[j];
                    diffX9_uranium = diffX9[k];
                    diffY9_uranium = diffY9[k];
                    ratioX8_uranium = ratioX8[j];
                    ratioY8_uranium = ratioY8[j];
                    ratioX9_uranium = ratioX9[k];
                    ratioY9_uranium = ratioY9[k];
                    sumX8_ampA8_uranium = sumX8_ampA8[j];
                    sumY8_ampA8_uranium = sumY8_ampA8[j];
                    sumX9_ampA9_uranium = sumX9_ampA9[k];
                    sumY9_ampA9_uranium = sumY9_ampA9[k];
                    tout_uranium->Fill();
                    counts_uranium++;
            }

        }
    }
}
    file_uranium->cd();
    tout_uranium->Write();
    file_uranium->Close();
    file_gold->cd();
    tout_gold->Write();
    file_gold->Close();
    file_uranium->Close();
    fin->Close();
    std::cout << "Processing completed." << std::endl;
    std::cout << "Final Counts gold: " << counts_gold << " Final Counts uranium: " << counts_uranium << std::endl;
    return 0;        
}