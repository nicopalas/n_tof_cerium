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
    TTree *tout = tin->CloneTree(0); // clone structure only

    // FIXED: Properly declare vectors for output branches
    vector<double> x_detector[NDET];
    vector<double> y_detector[NDET]; 
    vector<double> x_beam[NDET];
    vector<double> y_beam[NDET];
    vector<double> costheta[NDET-1];
    vector<double> phi[NDET-1];
    vector<double> x_target[NDET-1];
    vector<double> y_target[NDET-1];
    
    double vx_det[NDET];
    double vy_det[NDET];

    // FIXED: Initialize arrays
    for(int k=0;k<NDET;k++){
        vx_det[k] = 0;
        vy_det[k] = 0;
    }

    for(int k=0;k<NDET;k++){
        tout->Branch(Form("x_detector%d",k), &x_detector[k]);
        tout->Branch(Form("y_detector%d",k), &y_detector[k]);
        tout->Branch(Form("x_beam%d",k), &x_beam[k]);
        tout->Branch(Form("y_beam%d",k), &y_beam[k]);
        
        tout->Branch(Form("vX%d",k), &vx_det[k], Form("vX%d/D",k));
        tout->Branch(Form("vY%d",k), &vy_det[k], Form("vY%d/D",k));
        
        if(k < NDET-1) {
            tout->Branch(Form("x_target%d",k), &x_target[k]);
            tout->Branch(Form("y_target%d",k), &y_target[k]);
            tout->Branch(Form("costheta%d",k), &costheta[k]);
            tout->Branch(Form("phi%d",k), &phi[k]);
        }
    }

    // --- Loop events again to fill output tree with reconstructed positions & angles ---
    for(Long64_t ev=0; ev<nentries; ev++){
        tin->GetEntry(ev);

        // clear vectors for this event
        for(int k=0;k<NDET;k++){
            x_detector[k].clear(); 
            y_detector[k].clear(); 
            x_beam[k].clear(); 
            y_beam[k].clear();
            vx_det[k] = vX[k];  // Set velocities for this event
            vy_det[k] = vY[k];
        }
        for(int k=0;k<NDET-1;k++){
            costheta[k].clear(); 
            phi[k].clear();
            x_target[k].clear();
            y_target[k].clear();
        }

        // position reconstruction using calibrated speeds (vX, vY)
        for(int k=0; k<NDET; k++){
            if (mult[k] < 1 || multC1[k] < 1 || multC2[k] < 1 || multC3[k] < 1 || multC4[k] < 1) continue;

            for(int i=0; i<mult[k]; i++){
                // y from C1/C2
                for (int j=0; j<multC1[k] && j<MAXC; j++){
                    for (int l=0; l<multC2[k] && l<MAXC; l++){
                        double sumDiff = tofC2[k][l] + tofC1[k][j] - 2*tofA[k][i];
                        if(k>=6){
                            if(sumDiff>60 && sumDiff<140){
                                double yi = (tofC2[k][l] - tofC1[k][j]) * (vY[k]/2.0);
                                if (std::isfinite(yi) && fabs(yi) > 1e-6 && fabs(yi) < L_mm){
                                    y_detector[k].push_back(yi);
                                }
                            }
                        } else {
                            if(sumDiff>40 && sumDiff<130){
                                double yi = (tofC2[k][l] - tofC1[k][j]) * (vY[k]/2.0);
                                if (std::isfinite(yi) && fabs(yi) > 1e-6 && fabs(yi) < L_mm){
                                    y_detector[k].push_back(yi);
                                }
                            }
                        }
                    }
                }
                // x from C3/C4
                for (int j=0; j<multC3[k] && j<MAXC; j++){
                    for (int l=0; l<multC4[k] && l<MAXC; l++){
                        double sumDiff = tofC4[k][l] + tofC3[k][j] - 2*tofA[k][i];
                        if(k>=6){
                            if(sumDiff>60 && sumDiff<140){
                                double xi = (tofC4[k][l] - tofC3[k][j]) * (vX[k]/2.0);
                                if (std::isfinite(xi) && fabs(xi) > 1e-6){
                                    // FIXED: Removed incorrect offset calculation
                                    x_detector[k].push_back(xi);
                                }
                            }
                        } else {
                            if(sumDiff>40 && sumDiff<130){
                                double xi = (tofC4[k][l] - tofC3[k][j]) * (vX[k]/2.0);
                                if (std::isfinite(xi) && fabs(xi) > 1e-6){
                                    // FIXED: Removed incorrect offset calculation
                                    x_detector[k].push_back(xi);
                                }
                            }
                        }
                    }
                }
            } // end multiplicities for this detector
        } // end detectors loop for positions

        // FIXED: Angle calculation - properly handle the loops
        for(int k=0; k<NDET-1; k++){
            for (size_t i0 = 0; i0 < x_detector[k].size(); ++i0) {
                for (size_t i1 = 0; i1 < x_detector[k+1].size(); ++i1) {
                    for (size_t j0 = 0; j0 < y_detector[k].size(); ++j0) {
                        for (size_t j1 = 0; j1 < y_detector[k+1].size(); ++j1) {
                            double x0 = x_detector[k][i0];
                            double y0 = y_detector[k][j0];
                            double z0 = k * detector_spacing;

                            double x1 = x_detector[k+1][i1];
                            double y1 = y_detector[k+1][j1];
                            double z1 = (k+1) * detector_spacing;
                            
                            // FIXED: Proper target position calculation
                            double x0targ = x0 + (x1 - x0) / detector_spacing * (detector_spacing/2.0 + d_cathode_anode);
                            double y0targ = y0 + (y1 - y0) / detector_spacing * (detector_spacing/2.0 - d_cathode_anode);
                            
                            double x0p, y0p, z0p, x1p, y1p, z1p;
                            rotY(x0, y0, z0, TMath::Pi()/4.0, x0p, y0p, z0p);
                            rotY(x1, y1, z1, TMath::Pi()/4.0, x1p, y1p, z1p);

                            x_target[k].push_back(x0targ);
                            y_target[k].push_back(y0targ);

                            double Vx = x1p - x0p;
                            double Vy = y1p - y0p;
                            double Vz = z1p - z0p;

                            double Vmod = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
                            if(Vmod > 1e-9){
                                double cosTh = Vz / Vmod;
                                costheta[k].push_back(cosTh);
                            } else {
                                costheta[k].push_back(-9999.0);
                            }
                            phi[k].push_back(atan2(Vy, Vx));
                        }
                    }
                }
            }
        }

        tout->Fill();
    }
    cout<<"Processed "<<nentries<<" events.\n";

    fout->Write();
    fout->Close();
    fin->Close();

    cout<<"File "<<outfile<<" written \n";
}

