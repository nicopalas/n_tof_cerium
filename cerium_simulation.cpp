#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TNamed.h>
#include <sstream>
#include <cmath>
#include <iostream>

// --------- User constants and options -------------
const int Z_CN = 58;                 // Cerium Z
const int A_CN = 140;                // Assumed mass number (change if desired)
const double beta_def = 0.6;         // deformation
const int Nevents = 500000;          // number of events per tree
const unsigned int RNG_seed = 12345; // reproducibility
// --------------------------------------------------

// Simple Viola systematics (common simple parameterization).
double violaTKE(int Z, int A) {
    // TKE in MeV for the compound nucleus (approximate Viola)
    // TKE = 0.1189 * Z^2 / A^{1/3} + 7.3
    return 0.1189 * (double)Z * (double)Z / std::pow((double)A, 1.0/3.0) + 7.3;
}

// Solve for d such that Coulomb TKE == Viola TKE for symmetric split
double solve_d_for_symmetric_viola(int Zcn, int Acn, double beta) {
    int Z1 = Zcn / 2;
    int Z2 = Zcn - Z1;
    double Afactor = (double)Acn / (double)Zcn;
    int A1 = (int)std::round(Afactor * (double)Z1);
    int A2 = Acn - A1;

    double T_viola = violaTKE(Zcn, Acn); // MeV
    double numerator = 1.44 * (double)Z1 * (double)Z2; // MeV*fm
    double R_needed = numerator / T_viola; // fm

    double r1 = std::pow((double)A1, 1.0/3.0) * (1.0 + 2.0*beta/3.0);
    double r2 = std::pow((double)A2, 1.0/3.0) * (1.0 + 2.0*beta/3.0);
    double d = R_needed/1.2 - (r1 + r2);
    return d;
}

// sample Z from gaussian around mean with sigma
int sampleZFromGaussian(TRandom3 &rnd, double mean, double sigma) {
    double z = rnd.Gaus(mean, sigma);
    if (z < 1.0) z = 1.0;
    if (z > (double)(Z_CN - 1)) z = (double)(Z_CN - 1);
    int zi = (int)std::round(z);
    if (zi < 1) zi = 1;
    if (zi > Z_CN - 1) zi = Z_CN - 1;
    return zi;
}

int main() {
    std::cout << "Fission simulation for Ce (Z=" << Z_CN << ", A=" << A_CN << ")\n";
    std::cout << "Events per tree: " << Nevents << "\n";

    TRandom3 rnd(RNG_seed);

    // compute d to match Viola for symmetric splitting
    std::cout << "TKE using viola: " << violaTKE(Z_CN, A_CN) << " MeV\n";
    double d_sym = solve_d_for_symmetric_viola(Z_CN, A_CN, beta_def);
    int Z1_sym = Z_CN / 2;
    int Z2_sym = Z_CN - Z1_sym;
    double Afactor_sym = (double)A_CN / (double)Z_CN;
    int A1_sym = (int)std::round(Afactor_sym * (double)Z1_sym);
    int A2_sym = A_CN - A1_sym;

    double r1_sym = std::pow((double)A1_sym, 1.0/3.0) * (1.0 + 2.0*beta_def/3.0);
    double r2_sym = std::pow((double)A2_sym, 1.0/3.0) * (1.0 + 2.0*beta_def/3.0);
    double R_sym = 1.2 * (r1_sym + r2_sym + d_sym);
    double TKE_sym = 1.44 * (double)Z1_sym * (double)Z2_sym / R_sym;

    std::cout << "Symmetric split TKE (using d_sym, beta): " << TKE_sym << " MeV\n";
    std::cout << "Computed neck distance d (to match Viola for symmetric split): " << d_sym << " fm\n";

    // Prepare ROOT file and trees
    TFile *fout = new TFile("fission_Ce140.root","RECREATE");

    // Trees for symmetric and asymmetric
    TTree *tsym = new TTree("symmetricTree","Symmetric-split sample");
    TTree *tasym = new TTree("asymmetricTree","Asymmetric (Z=34 shell) sample");

    // Event branches - EXACT names as requested
    int event;
    int Z1, Z2, A1post, A2post;
    double tke_tot, tke1, tke2;

    tsym->Branch("event", &event, "event/I");
    tsym->Branch("Z1", &Z1, "Z1/I");
    tsym->Branch("Z2", &Z2, "Z2/I");
    tsym->Branch("A1", &A1post, "A1/I");
    tsym->Branch("A2", &A2post, "A2/I");
    tsym->Branch("TKE", &tke_tot, "TKE/D");
    tsym->Branch("KE1", &tke1, "KE1/D");
    tsym->Branch("KE2", &tke2, "KE2/D");

    // Asymmetric tree branches - EXACT names as requested
    tasym->Branch("event", &event, "event/I");
    tasym->Branch("Z1", &Z1, "Z1/I");
    tasym->Branch("Z2", &Z2, "Z2/I");
    tasym->Branch("A1", &A1post, "A1/I");
    tasym->Branch("A2", &A2post, "A2/I");
    tasym->Branch("TKE", &tke_tot, "TKE/D");
    tasym->Branch("KE1", &tke1, "KE1/D");
    tasym->Branch("KE2", &tke2, "KE2/D");

    // Histograms
    TH1D *hZ_sym = new TH1D("hZ_sym","Z distribution symmetric;Z;counts", Z_CN, 0.5, Z_CN+0.5);
    TH1D *hA_sym = new TH1D("hA_sym","A distribution symmetric;A;counts", A_CN, 0.5, A_CN+0.5);
    TH2D *hTKE_vs_Z_sym = new TH2D("hTKE_vs_Z_sym","TKE vs Z1 symmetric;Z1;TKE (MeV)", Z_CN, 0.5, Z_CN+0.5, 200, 0, 150);
    TH2D *hTKE_vs_A_sym = new TH2D("hTKE_vs_A_sym","TKE vs A1 symmetric;A1;TKE (MeV)", A_CN, 0.5, A_CN+0.5, 200, 0, 150);

    TH1D *hZ_asym = new TH1D("hZ_asym","Z distribution asymmetric;Z;counts", Z_CN, 0.5, Z_CN+0.5);
    TH1D *hA_asym = new TH1D("hA_asym","A distribution asymmetric;A;counts", A_CN, 0.5, A_CN+0.5);
    TH2D *hTKE_vs_Z_asym = new TH2D("hTKE_vs_Z_asym","TKE vs Z1 asymmetric;Z1;TKE (MeV)", Z_CN, 0.5, Z_CN+0.5, 200, 0, 150);
    TH2D *hTKE_vs_A_asym = new TH2D("hTKE_vs_A_asym","TKE vs A1 asymmetric;A1;TKE (MeV)", A_CN, 0.5, A_CN+0.5, 200, 0, 150);

    double Afactor = (double)A_CN / (double)Z_CN;

    // Temporary variables for calculations
    double r1, r2, R;

    // ----- Symmetric sampling -----
    double mean_sym_Z = (double)Z_CN / 2.0;
    double sigma_sym_Z = 4.0;

    for (int i = 0; i < Nevents; ++i) {
        event = i;
        Z1 = sampleZFromGaussian(rnd, mean_sym_Z, sigma_sym_Z);
        Z2 = Z_CN - Z1;
        A1post = (int)std::round(Afactor * (double)Z1);
        if (A1post < 1) A1post = 1;
        if (A1post > A_CN-1) A1post = A_CN-1;
        A2post = A_CN - A1post;

        r1 = std::pow((double)A1post, 1.0/3.0) * (1.0 + 2.0*beta_def/3.0);
        r2 = std::pow((double)A2post, 1.0/3.0) * (1.0 + 2.0*beta_def/3.0);
        R = 1.2 * (r1 + r2 + d_sym);
        tke_tot = 1.44 * (double)Z1 * (double)Z2 / R;

        // per-fragment kinetic energies
        tke1 = (double)A2post / (double)(A1post + A2post) * tke_tot;
        tke2 = (double)A1post / (double)(A1post + A2post) * tke_tot;

        tsym->Fill();
        hZ_sym->Fill(Z1);
        hA_sym->Fill(A1post);
        hTKE_vs_Z_sym->Fill(Z1, tke_tot);
        hTKE_vs_A_sym->Fill(A1post, tke_tot);
    }

    // ----- Asymmetric sampling -----
    double Zshell = 34.0;
    double sigma_asym = 2.5;

    for (int i = 0; i < Nevents; ++i) {
        event = i;
        Z1 = sampleZFromGaussian(rnd, Zshell, sigma_asym);
        Z2 = Z_CN - Z1;

        A1post = (int)std::round(Afactor * (double)Z1);
        if (A1post < 1) A1post = 1;
        if (A1post > A_CN-1) A1post = A_CN-1;
        A2post = A_CN - A1post;

        r1 = std::pow((double)A1post, 1.0/3.0) * (1.0 + 2.0*beta_def/3.0);
        r2 = std::pow((double)A2post, 1.0/3.0) * (1.0 + 2.0*beta_def/3.0);
        R = 1.2 * (r1 + r2 + d_sym);
        tke_tot = 1.44 * (double)Z1 * (double)Z2 / R;

        // per-fragment kinetic energies
        tke1 = (double)A2post / (double)(A1post + A2post) * tke_tot;
        tke2 = (double)A1post / (double)(A1post + A2post) * tke_tot;

        tasym->Fill();
        hZ_asym->Fill(Z1);
        hA_asym->Fill(A1post);
        hTKE_vs_Z_asym->Fill(Z1, tke_tot);
        hTKE_vs_A_asym->Fill(A1post, tke_tot);
    }

    // Save everything
    fout->cd();
    tsym->Write();
    tasym->Write();

    hZ_sym->Write();
    hA_sym->Write();
    hTKE_vs_Z_sym->Write();
    hTKE_vs_A_sym->Write();

    hZ_asym->Write();
    hA_asym->Write();
    hTKE_vs_Z_asym->Write();
    hTKE_vs_A_asym->Write();

    double T_viola = violaTKE(Z_CN, A_CN);
    {
        std::ostringstream oss;
        oss << "d_sym=" << d_sym;
        TNamed dname("d_sym", oss.str().c_str());
        dname.Write();

        oss.str(""); 
        oss.clear();
        oss << "T_viola=" << T_viola;
        TNamed tname("T_viola", oss.str().c_str());
        tname.Write();
    }

    std::cout << "Wrote output file: fission_Ce140.root\n";

    fout->Close();
    delete fout;
    return 0;
}