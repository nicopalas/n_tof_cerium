#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

void plot() {

    // --- Open file and get tree ---
    TFile *f = TFile::Open("coincidences.root");
    if (!f || f->IsZombie()) { 
        std::cout << "Cannot open coincidences.root\n"; 
        return; 
    }

    TTree *tree = (TTree*)f->Get("training");
    if (!tree) { 
        std::cout << "Cannot find TTree 'training'\n"; 
        return; 
    }

    // --- Create alias for neutron energy ---
    tree->SetAlias(
        "neutron_energy_8",
        "939.565420 * (1.0 / TMath::Sqrt(1 - TMath::Power(184.5 / (299792458.0 * (tofA8 + 726.616 + 185/0.3) * 1e-9), 2)) - 1.0)"
    );

    // --- Activate only needed branches ---
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("tofA8", 1);
    tree->SetBranchStatus("dtA", 1);

    // --- Create log-spaced bin edges for neutron energy ---
    int nbinsX = 800;
    double xmin = 0.6;
    double xmax = 1000.0;

    double logmin = TMath::Log10(xmin);
    double logmax = TMath::Log10(xmax);

    std::vector<double> edgesX(nbinsX + 1);
    for (int i = 0; i <= nbinsX; i++)
        edgesX[i] = TMath::Power(10, logmin + (logmax - logmin) * i / nbinsX); // Logarithmic spacing between bins

    // --- Create 2D histogram ---
    TH2D *h2 = new TH2D("h2", "",
                        nbinsX, edgesX.data(),
                        500, -10, 10);  // adjust dtA range if needed

    // --- Draw / fill histogram ---
    TCanvas *c1 = new TCanvas("c1", "Plot", 900, 700);
    tree->Draw("dtA:neutron_energy_8 >> h2", "neutron_energy_8<1000", "COLZ");
    gPad->SetLogx();

    // --- Axis labels ---
    h2->GetXaxis()->SetTitle("E_{n} (MeV)");
    h2->GetYaxis()->SetTitle("t_{1} - t_{0} (ns)");

    c1->Update();
}
