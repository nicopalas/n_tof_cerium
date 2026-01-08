#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>

int main() {

    // =======================
    // ===== GOLD ============
    // =======================

    TFile *fin_gold = TFile::Open("events_gold.root", "UPDATE");
    if (!fin_gold || fin_gold->IsZombie()) {
        std::cerr << "Error opening events_gold.root\n";
        return 1;
    }

    TTree *tin_gold = (TTree*)fin_gold->Get("events_gold");
    if (!tin_gold) {
        std::cerr << "Tree events_gold not found\n";
        return 2;
    }

    double x7, y7, x8, y8;
    tin_gold->SetBranchAddress("x7", &x7);
    tin_gold->SetBranchAddress("y7", &y7);
    tin_gold->SetBranchAddress("x8", &x8);
    tin_gold->SetBranchAddress("y8", &y8);

    Long64_t n_gold = tin_gold->GetEntries();

    double sum_x7 = 0, sum_y7 = 0, sum_x8 = 0, sum_y8 = 0;

    for (Long64_t i = 0; i < n_gold; i++) {
        tin_gold->GetEntry(i);
        sum_x7 += x7;
        sum_y7 += y7;
        sum_x8 += x8;
        sum_y8 += y8;
    }

    double mean_x7 = sum_x7 / n_gold;
    double mean_y7 = sum_y7 / n_gold;
    double mean_x8 = sum_x8 / n_gold;
    double mean_y8 = sum_y8 / n_gold;

    TTree *tout_gold = tin_gold->CloneTree(0);

    double x7_corr, y7_corr, x8_corr, y8_corr;
    double cos_theta_gold, phi_gold;

    tout_gold->Branch("x7_corr", &x7_corr, "x7_corr/D");
    tout_gold->Branch("y7_corr", &y7_corr, "y7_corr/D");
    tout_gold->Branch("x8_corr", &x8_corr, "x8_corr/D");
    tout_gold->Branch("y8_corr", &y8_corr, "y8_corr/D");
    tout_gold->Branch("cos_theta", &cos_theta_gold, "cos_theta/D");
    tout_gold->Branch("phi", &phi_gold, "phi/D");

    for (Long64_t i = 0; i < n_gold; i++) {
        tin_gold->GetEntry(i);

        x7_corr = x7 - (mean_x7 + 17.0);
        y7_corr = y7 - mean_y7;

        x8_corr = x8 - (mean_x8 - 17.0);
        y8_corr = y8 - mean_y8;

        cos_theta_gold =
            ((-x8_corr + x7_corr + 34.0) * 0.7071 + 34.0 * 0.7071) /
            TMath::Sqrt(
                TMath::Power(y8_corr - y7_corr, 2) +
                TMath::Power(-x8_corr + x7_corr + 34.0, 2) +
                TMath::Power(34.0, 2)
            );

        phi_gold = TMath::ATan2(
            (y8_corr - y7_corr),
            -(-x8_corr + x7_corr + 34.0) * 0.7071 + 34.0 * 0.7071
        );

        tout_gold->Fill();
    }

    fin_gold->cd();
    tout_gold->Write("events_gold", TObject::kOverwrite);
    fin_gold->Close();

    // =======================
    // ===== URANIUM =========
    // =======================

    TFile *fin_u = TFile::Open("events_uranium.root", "UPDATE");
    if (!fin_u || fin_u->IsZombie()) {
        std::cerr << "Error opening events_uranium.root\n";
        return 3;
    }

    TTree *tin_u = (TTree*)fin_u->Get("events_uranium");
    if (!tin_u) {
        std::cerr << "Tree events_uranium not found\n";
        return 4;
    }

    double x8_u, y8_u, x9_u, y9_u;
    tin_u->SetBranchAddress("x8", &x8_u);
    tin_u->SetBranchAddress("y8", &y8_u);
    tin_u->SetBranchAddress("x9", &x9_u);
    tin_u->SetBranchAddress("y9", &y9_u);

    Long64_t n_u = tin_u->GetEntries();

    double sum_x8u = 0, sum_y8u = 0, sum_x9u = 0, sum_y9u = 0;

    for (Long64_t i = 0; i < n_u; i++) {
        tin_u->GetEntry(i);
        sum_x8u += x8_u;
        sum_y8u += y8_u;
        sum_x9u += x9_u;
        sum_y9u += y9_u;
    }

    double mean_x8_u = sum_x8u / n_u;
    double mean_y8_u = sum_y8u / n_u;
    double mean_x9_u = sum_x9u / n_u;
    double mean_y9_u = sum_y9u / n_u;

    TTree *tout_u = tin_u->CloneTree(0);

    double x8_corr_u, y8_corr_u, x9_corr_u, y9_corr_u;
    double cos_theta_u, phi_u;

    tout_u->Branch("x8_corr", &x8_corr_u, "x8_corr/D");
    tout_u->Branch("y8_corr", &y8_corr_u, "y8_corr/D");
    tout_u->Branch("x9_corr", &x9_corr_u, "x9_corr/D");
    tout_u->Branch("y9_corr", &y9_corr_u, "y9_corr/D");
    tout_u->Branch("cos_theta", &cos_theta_u, "cos_theta/D");
    tout_u->Branch("phi", &phi_u, "phi/D");

    for (Long64_t i = 0; i < n_u; i++) {
        tin_u->GetEntry(i);
        if (i==0){
            std::cout << "Mean x8 uranium: " << mean_x8_u << "\n";
            std::cout << "Mean y8 uranium: " << mean_y8_u << "\n";
            std::cout << "Mean x9 uranium: " << mean_x9_u << "\n";
            std::cout << "Mean y9 uranium: " << mean_y9_u << "\n";
        }

        x8_corr_u = x8_u-(mean_x8_u+17.0);
        y8_corr_u = y8_u - mean_y8_u;

        x9_corr_u = x9_u - (mean_x9_u - 17.0);
        y9_corr_u = y9_u - mean_y9_u;

        cos_theta_u =
            ((-x9_corr_u + x8_corr_u + 34.0) * 0.7071 + 34.0 * 0.7071) /
            TMath::Sqrt(
                TMath::Power(y9_corr_u - y8_corr_u, 2) +
                TMath::Power(-x9_corr_u + x8_corr_u + 34.0, 2) +
                TMath::Power(34.0, 2)
            );

        phi_u = TMath::ATan2(
            (y9_corr_u - y8_corr_u),
            -(-x9_corr_u + x8_corr_u + 34.0) * 0.7071 + 34.0 * 0.7071
        );

        tout_u->Fill();
    }

    fin_u->cd();
    tout_u->Write("events_uranium", TObject::kOverwrite);
    fin_u->Close();

    std::cout << "Correcciones aplicadas correctamente usando medias exactas del TTree.\n";
    return 0;
}


