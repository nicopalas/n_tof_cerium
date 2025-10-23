#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TVector3.h>
#include <TMath.h>
#include <TLine.h>
#include <TCanvas.h>

#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

const int NDET   = 10;
const int MAXA   = 5;    // max multiplicity anodes
const int MAXC   = 20;   // max multiplicity cathodes
const double L_mm   = 200.0; // detector size (mm)
const double detector_spacing = 50.0; // spacing along z between consecutive detectors (mm)
const double d_cathode_anode = 3.2; // distance between cathode and anode (mm)

// rotation around Y function
void rotY(double x, double y, double z, double ang,
          double &xp, double &yp, double &zp)
{
    double ca = cos(ang);
    double sa = sin(ang);
    xp =  x*ca + z*sa;
    yp =  y;
    zp = -x*sa + z*ca;
}

void positions(const char *infile="input.root",
               const char *treename="T",
               const char *outfile="output.root")
{
    // --- Input ---
    TFile *fin = TFile::Open(infile,"READ");
    if(!fin || fin->IsZombie()){ cerr<<"No se pudo abrir "<<infile<<"\n"; return; }
    TTree *tin = (TTree*)fin->Get(treename);
    if(!tin){ cerr<<"No se encontró el árbol "<<treename<<"\n"; return; }
    Long64_t nentries = tin->GetEntries();

    // branches arrays - FIXED: Added proper array declarations
    Double_t tofA[NDET][MAXA];
    Double_t tofC1[NDET][MAXC], tofC2[NDET][MAXC], tofC3[NDET][MAXC], tofC4[NDET][MAXC];
    Int_t mult[NDET];
    Int_t multC1[NDET], multC2[NDET], multC3[NDET], multC4[NDET];

    // set branch addresses (assumes branches named tof0, tof01, mult0, mult01, etc.)
    for(int k=0;k<NDET;k++){
        tin->SetBranchAddress(Form("tof%d",k),  tofA[k]);
        tin->SetBranchAddress(Form("tof%d1",k), tofC1[k]);
        tin->SetBranchAddress(Form("tof%d2",k), tofC2[k]);
        tin->SetBranchAddress(Form("tof%d3",k), tofC3[k]);
        tin->SetBranchAddress(Form("tof%d4",k), tofC4[k]);
        tin->SetBranchAddress(Form("mult%d",k), &mult[k]);
        tin->SetBranchAddress(Form("mult%d1",k), &multC1[k]);
        tin->SetBranchAddress(Form("mult%d2",k), &multC2[k]);
        tin->SetBranchAddress(Form("mult%d3",k), &multC3[k]);
        tin->SetBranchAddress(Form("mult%d4",k), &multC4[k]);
    }

    // histograms to find reflection peaks
    vector<TH1D*> hX(NDET), hY(NDET);
    for(int k=0;k<NDET;k++){
        hX[k] = new TH1D(Form("hX%d",k), Form("hX%d",k), 2000, -500, 500);
        hY[k] = new TH1D(Form("hY%d",k), Form("hY%d",k), 2000, -500, 500);
    }

    for(Long64_t ev=0; ev<nentries; ev++){
        tin->GetEntry(ev);
        for(int k=0;k<NDET;k++){
            if (mult[k]<1) continue; // no hits on anode
            for(int i=0;i<mult[k];i++){
                // Y cathodes (C1,C2)
                for(int j=0;j<multC1[k] && j<MAXC;j++){
                    for(int l=0;l<multC2[k] && l<MAXC;l++){
                        double sumDiff = tofC2[k][l] + tofC1[k][j] - 2*tofA[k][i];
                        double y = (tofC2[k][l] - tofC1[k][j]);
                        if(k>=6){
                            if(sumDiff>50 && sumDiff<150) hY[k]->Fill(y);
                        } else {
                            if(sumDiff>40 && sumDiff<130) hY[k]->Fill(y);
                        }
                    }
                }
                // X cathodes (C3,C4)
                for(int j=0;j<multC3[k] && j<MAXC;j++){
                    for(int l=0;l<multC4[k] && l<MAXC;l++){
                        double sumDiff = tofC4[k][l] + tofC3[k][j] - 2*tofA[k][i];
                        double x = (tofC4[k][l] - tofC3[k][j]);
                        if(k>=6){
                            if(sumDiff>50 && sumDiff<150) hX[k]->Fill(x);
                        } else {
                            if(sumDiff>40 && sumDiff<130) hX[k]->Fill(x);
                        }
                    }
                }
            }
        }
    }

    // speeds of delay lines (mm per time-unit)
    vector<double> vX(NDET,0), vY(NDET,0);

    TCanvas *cCalib = new TCanvas("cCalib","Calibracion Delay Lines",1200,800);
    cCalib->Divide(2, NDET);

    for(int k=0;k<NDET;k++){
        // peak finding
        TSpectrum sx(40), sy(40);
        sx.Search(hX[k],2,"nodraw",0.01);
        sy.Search(hY[k],2,"nodraw",0.01);

        // collect peaks around reflections (expected near +/-110)
        vector<double> peaksX, peaksY;
        for(int j=0;j<sx.GetNPeaks();j++){
            double p = sx.GetPositionX()[j];
            if(fabs(p) > 90 && fabs(p) < 120) peaksX.push_back(p);
        }
        for(int j=0;j<sy.GetNPeaks();j++){
            double p = sy.GetPositionX()[j];
            if(fabs(p) > 90 && fabs(p) < 120) peaksY.push_back(p);
        }

        sort(peaksX.begin(), peaksX.end());
        sort(peaksY.begin(), peaksY.end());

        double px=0, py=0;
        if(peaksX.size()>=2) {
            // pick one near +110 and one near -110 (if present) and average abs positions
            auto it_pos = min_element(peaksX.begin(), peaksX.end(), [](double a, double b){
                return fabs(a - 110) < fabs(b - 110);
            });
            auto it_neg = min_element(peaksX.begin(), peaksX.end(), [](double a, double b){
                return fabs(a + 110) < fabs(b + 110);
            });
            px = 0.5 * (fabs(*it_pos) + fabs(*it_neg));
        } else if(peaksX.size()==1) {
            px = fabs(peaksX[0]);
        }

        if(peaksY.size()>=2) {
            auto it_pos = min_element(peaksY.begin(), peaksY.end(), [](double a, double b){
                return fabs(a - 110) < fabs(b - 110);
            });
            auto it_neg = min_element(peaksY.begin(), peaksY.end(), [](double a, double b){
                return fabs(a + 110) < fabs(b + 110);
            });
            py = 0.5 * (fabs(*it_pos) + fabs(*it_neg));
        } else if(peaksY.size()==1) {
            py = fabs(peaksY[0]);
        }

        if(px>1e-6) vX[k] = L_mm/px;
        if(py>1e-6) vY[k] = L_mm/py;

        // plotting with lines marking the chosen peaks
        cCalib->cd(2*k+1);
        gPad->SetLogy();
        hX[k]->Draw();
        if(peaksX.size() >= 1) {
            for(double p : peaksX){
                TLine *l = new TLine(p, 1, p, hX[k]->GetMaximum());
                l->SetLineColor(kRed);
                l->SetLineStyle(1);
                l->Draw("same");
            }
        }

        cCalib->cd(2*k+2);
        gPad->SetLogy();
        hY[k]->Draw();
        if(peaksY.size() >= 1) {
            for(double p : peaksY){
                TLine *l = new TLine(p, 1, p, hY[k]->GetMaximum());
                l->SetLineColor(kBlue);
                l->SetLineStyle(2);
                l->Draw("same");
            }
        }

        cout<<"det "<<k<<" px="<<px<<" vX="<<vX[k]
            <<"   py="<<py<<" vY="<<vY[k]<<"\n";
    } // end per-detector calibration

    cCalib->SaveAs("calibracion_peaks.pdf");

    // --- Output ---
    TFile *fout = TFile::Open(outfile,"RECREATE");

    // new branches for positions
    

