#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

const int NDET = 10;      // detectors 0–9
const int MAXHITS = 5;   // must be >= maximum multiplicity per detector

void calibrate_amplitudes(const char* infile) {
    // Open file + tree
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Error opening file\n"; return; }
    TTree *t = (TTree*)fin->Get("coincidences");
    if (!t) { std::cerr << "No TTree 'coincidences' found\n"; fin->Close(); return; }

    // Branch setup
    Int_t mult = 0;
    Int_t mults[NDET] = {0};
    Float_t amps[NDET][MAXHITS];
    t->SetBranchAddress("mult", &mult);
    for (int d=0; d<NDET; d++) {
        TString brmult = Form("mult%d", d);
        TString bramp  = Form("amp%d", d);
        t->SetBranchAddress(brmult, &mults[d]);
        t->SetBranchAddress(bramp, amps[d]); // amps[d] decays to Float_t*
    }

    // Histograms for γ flash
    TH1F *hflash[NDET];
    for (int d=0; d<NDET; d++) {
        hflash[d] = new TH1F(Form("hflash%d",d), Form("Gamma flash det %d",d),
                             400, 2000, 30000);
        hflash[d]->Sumw2();
    }

    // Fill γ flash events
    Long64_t nev = t->GetEntries();
    for (Long64_t i=0;i<nev;i++) {
        t->GetEntry(i);
        if (mult <= 8) continue; // gamma flash condition
        for (int d=0; d<NDET; d++) {
            for (int j=0;j<mults[d] && j<MAXHITS; j++) {
                hflash[d]->Fill(amps[d][j]);
            }
        }
    }

    // Fit γ flash peaks
    TSpectrum spectrum(8);
    double peakpos1[NDET]; for (int d=0; d<NDET; d++) peakpos1[d] = 0;

    TCanvas *cAll = new TCanvas("cAll", "Gamma flash all detectors", 1200, 800);
    cAll->Divide(5, 2);

    for (int d=0; d<NDET; d++) {
        cAll->cd(d+1);
        hflash[d]->Draw();

        Int_t nfound = spectrum.Search(hflash[d], 2, "", 0.4);
        double *xpeaks = spectrum.GetPositionX();
        std::vector<double> peaks(xpeaks, xpeaks + nfound);
        std::sort(peaks.begin(), peaks.end());

        double p1 = 0;
        for (double p : peaks) if (TMath::Abs(p - 6500) < 1500) { p1 = p; break; }
        if (p1 == 0) {
            int binMax = hflash[d]->GetMaximumBin();
            p1 = hflash[d]->GetBinCenter(binMax);
        }

        TString fname = Form("f1_det%d", d);
        TF1 *f1 = new TF1(fname, "gaus", p1-2500, p1+2500);
        f1->SetParameters(hflash[d]->GetMaximum(), p1, 500.0);
        hflash[d]->Fit(f1, "RQ");
        peakpos1[d] = f1->GetParameter(1);

        std::cout << "Detector " << d << " peak at " << peakpos1[d] << "\n";
    }

    cAll->Update();
    cAll->SaveAs("gamma_flash_peaks.pdf");
    delete cAll;

    // Reference peak
    double refPeak = 0; int count=0;
    for (int d=0; d<NDET; d++) if (peakpos1[d]>0) { refPeak += peakpos1[d]; count++; }
    refPeak /= count;
    std::cout << "Reference peak = " << refPeak << "\n";

    // Calibration factors
    double calibFactor[NDET];
    for (int d=0; d<NDET; d++) {
        calibFactor[d] = (peakpos1[d]>0) ? refPeak/peakpos1[d] : 1.0;
        std::cout << "Detector " << d << " calib = " << calibFactor[d] << "\n";
    }

    // Histograms for calibrated spectra (all events)
    TH1F *hcal[NDET];
    for (int d=0; d<NDET; d++) {
        hcal[d] = new TH1F(Form("hcal%d",d), "", 400, 0, 30000);
        hcal[d]->Sumw2();
    }

    // Histograms for background (mult < 8)
    TH1F *hbg[NDET];
    for (int d=0; d<NDET; d++) {
        hbg[d] = new TH1F(Form("hbg%d",d), "", 400, 0, 10000);
        hbg[d]->Sumw2();
    }

    // Fill calibrated + background
    for (Long64_t i=0;i<nev;i++) {
        t->GetEntry(i);
        for (int d=0; d<NDET; d++) {
            for (int j=0;j<mults[d] && j<MAXHITS; j++) {
                double ampCal = amps[d][j] * calibFactor[d];
                if (mult>8){
                hcal[d]->Fill(ampCal);
                }
                else{
                    hbg[d]->Fill(ampCal);
                }
            }
        }
    }

    // Overlay calibrated spectra (linear scale, full y range)
    TCanvas *cComp = new TCanvas("cComp", "Calibrated spectra comparison", 1000, 800);
    hcal[0]->SetLineColor(kBlack);
    hcal[0]->Draw("HIST");
    for (int d=1; d<NDET; d++) {
        hcal[d]->SetLineColor(TColor::GetColorPalette((d * 50) % TColor::GetNumberOfColors()));
        hcal[d]->Draw("HIST SAME");
    }
    TLegend *leg1 = new TLegend(0.7, 0.6, 0.92, 0.92);
    for (int d=0; d<NDET; d++) leg1->AddEntry(hcal[d], Form("Det %d", d), "l");
    hcal[0]->SetStats(0);
    leg1->Draw();
    hcal[0]->SetMaximum(12000);
    cComp->SaveAs("calibrated_spectra.pdf");
    delete cComp;

    // Overlay background spectra (log y)
    TCanvas *cBg = new TCanvas("cBg", "Background spectra", 1000, 800);
    cBg->SetLogy();
    hbg[0]->SetLineColor(kBlack);
    hbg[0]->Draw("HIST");
    for (int d=1; d<NDET; d++) {
        hbg[d]->SetLineColor(TColor::GetColorPalette((d * 45) % TColor::GetNumberOfColors()));
        hbg[d]->Draw("HIST SAME");
    }
    TLegend *leg2 = new TLegend(0.7, 0.6, 0.92, 0.92);
    for (int d=0; d<NDET; d++) leg2->AddEntry(hbg[d], Form("Det %d", d), "l");
    leg2->Draw();
    hbg[0]->SetStats(0);
    cBg->SaveAs("background_spectra.pdf");
    delete cBg;

    fin->Close();
}


