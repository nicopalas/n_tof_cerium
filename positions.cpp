

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TVector3.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TCanvas.h>
using namespace std;

const int NDET   = 10;
const int MAXA   = 5;    // max multiplicity anodes
const int MAXC   = 20;   // max multiplicity cathodes
const double L_mm   = 200.0; // delay line length
const double z_dist = 25.0;  // target–detector distance

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

void anode_cathode_positions(const char *infile="input.root",
                             const char *treename="T",
                             const char *outfile="output.root")
{
    // --- Input ---
    TFile *fin = TFile::Open(infile,"READ");
    if(!fin || fin->IsZombie()){ cerr<<"No se pudo abrir "<<infile<<"\n"; return; }
    TTree *tin = (TTree*)fin->Get(treename);
    if(!tin){ cerr<<"No se encontró el árbol "<<treename<<"\n"; return; }
    Long64_t nentries = tin->GetEntries();

    Double_t tofA[NDET][MAXA];
    Double_t tofC1[NDET][MAXC], tofC2[NDET][MAXC], tofC3[NDET][MAXC], tofC4[NDET][MAXC];
    Int_t mult[NDET];

    for(int k=0;k<NDET;k++){
        tin->SetBranchAddress(Form("tof%d",k),  tofA[k]);
        tin->SetBranchAddress(Form("tof%d1",k), tofC1[k]);
        tin->SetBranchAddress(Form("tof%d2",k), tofC2[k]);
        tin->SetBranchAddress(Form("tof%d3",k), tofC3[k]);
        tin->SetBranchAddress(Form("tof%d4",k), tofC4[k]);
        tin->SetBranchAddress(Form("mult%d",k), &mult[k]);
    }
vector<TH1D*> hX(NDET), hY(NDET);
for(int k=0;k<NDET;k++){
    hX[k] = new TH1D(Form("hX%d",k),"",2000,-500,500); // delay lines have 200 mm, so these limits should contain the important events
    hY[k] = new TH1D(Form("hY%d",k),"",2000,-500,500);
}
for(Long64_t ev=0; ev<nentries; ev++){
    tin->GetEntry(ev);
    for(int k=0;k<NDET;k++){
        for(int i=0;i<mult[k];i++){
            hX[k]->Fill(tofC4[k][i] - tofC3[k][i]); // fill histograms
            hY[k]->Fill(tofC2[k][i] - tofC1[k][i]);
        }
    }
}


// speeds of delay lines
vector<double> vX(NDET,0), vY(NDET,0);
TCanvas *cCalib = new TCanvas("cCalib","Calibracion Delay Lines",1200,800);
cCalib->Divide(2,NDET);

for(int k=0;k<NDET;k++){
    // peak finding
    TSpectrum sx(20); // maximum of 20 peaks...
    TSpectrum sy(20);
    sx.Search(hX[k],2,"nodraw",0.05); // we do not expect big peaks so low discrimination... we will filter later by position
    sy.Search(hY[k],2,"nodraw",0.05);

    // peak searching
    vector<double> peaksX, peaksY;
    for(int j=0;j<sx.GetNPeaks();j++){
        double p = sx.GetPositionX()[j];
        if(fabs(p) > 90 && fabs(p) < 120) peaksX.push_back(p); // avoid peaks in the center (not related to reflections)
    }
    for(int j=0;j<sy.GetNPeaks();j++){
        double p = sy.GetPositionX()[j];
        if(fabs(p) > 90 && fabs(p) < 120) peaksY.push_back(p);
    }

    sort(peaksX.begin(), peaksX.end()); //sorting by positions
    sort(peaksY.begin(), peaksY.end());

    double px=0, py=0;
    if(peaksX.size()>=2) {
        // Find peaks closest to +110 and -110
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

    // plotting
    cCalib->cd(2*k+1);
    gPad->SetLogy();
    hX[k]->Draw();
    // Draw a line for the peaks used for speed calculation
    if(peaksX.size() >= 2) {
        auto it_pos = min_element(peaksX.begin(), peaksX.end(), [](double a, double b){
            return fabs(a - 110) < fabs(b - 110);
        });
        auto it_neg = min_element(peaksX.begin(), peaksX.end(), [](double a, double b){
            return fabs(a + 110) < fabs(b + 110);
        });
        TLine *l_pos = new TLine(*it_pos, 0, *it_pos, hX[k]->GetMaximum());
        l_pos->SetLineColor(kRed); l_pos->SetLineStyle(1); l_pos->Draw("same");
        TLine *l_neg = new TLine(*it_neg, 0, *it_neg, hX[k]->GetMaximum());
        l_neg->SetLineColor(kRed); l_neg->SetLineStyle(1); l_neg->Draw("same");
    } else if(peaksX.size() == 1) {
        TLine *l = new TLine(peaksX[0], 0, peaksX[0], hX[k]->GetMaximum());
        l->SetLineColor(kRed); l->SetLineStyle(1); l->Draw("same");
    }

    cCalib->cd(2*k+2);
    gPad->SetLogy(); // logarithmic scale
    hY[k]->Draw();

    if(peaksY.size() >= 2) {
        auto it_pos = min_element(peaksY.begin(), peaksY.end(), [](double a, double b){
            return fabs(a - 110) < fabs(b - 110);
        });
        auto it_neg = min_element(peaksY.begin(), peaksY.end(), [](double a, double b){
            return fabs(a + 110) < fabs(b + 110);
        });
        TLine *l_pos = new TLine(*it_pos, 0, *it_pos, hY[k]->GetMaximum());
        l_pos->SetLineColor(kBlue); l_pos->SetLineStyle(2); l_pos->Draw("same");
        TLine *l_neg = new TLine(*it_neg, 0, *it_neg, hY[k]->GetMaximum());
        l_neg->SetLineColor(kBlue); l_neg->SetLineStyle(2); l_neg->Draw("same");
    } else if(peaksY.size() == 1) {
        TLine *l = new TLine(peaksY[0], 0, peaksY[0], hY[k]->GetMaximum());
        l->SetLineColor(kBlue); l->SetLineStyle(2); l->Draw("same");
    }

    cout<<"det "<<k<<" px="<<px<<" vX="<<vX[k]
        <<"   py="<<py<<" vY="<<vY[k]<<"\n";
}

// save the figure
cCalib->SaveAs("calibracion_peaks.pdf");


    // --- Output ---
    TFile *fout = TFile::Open(outfile,"RECREATE");
    TTree *tout = tin->CloneTree(0); // clone the initial tree structure

    // Nuevos branches
    vector<double> x[NDET], y[NDET], costheta[NDET-1], phi[NDET-1];
    double vx_det[NDET], vy_det[NDET];

    for(int k=0;k<NDET;k++){
        tout->Branch(Form("x%d",k), &x[k]);
        tout->Branch(Form("y%d",k), &y[k]);
        tout->Branch(Form("vX%d",k), &vx_det[k], Form("vX%d/D",k));
        tout->Branch(Form("vY%d",k), &vy_det[k], Form("vY%d/D",k));
        if(k < NDET-1) {
            tout->Branch(Form("costheta%d",k), &costheta[k]);
            tout->Branch(Form("phi%d",k), &phi[k]);
        }
    }

    // angles and positions calculations
    for(Long64_t ev=0; ev<nentries; ev++){
        tin->GetEntry(ev);

        // clear vectors
        for(int k=0;k<NDET;k++){ x[k].clear(); y[k].clear(); }
        for(int k=0;k<NDET-1;k++){ costheta[k].clear(); phi[k].clear(); }

        // position reconstruction
        for(int k=0; k<NDET; k++){
            for(int i=0; i<mult[k]; i++){
                double yi = (tofC2[k][i] - tofC1[k][i]) * (vY[k]/2.0);
                double xi = (tofC4[k][i] - tofC3[k][i]) * (vX[k]/2.0);
                x[k].push_back(xi);
                y[k].push_back(yi);
            }
            vx_det[k] = vX[k];
            vy_det[k] = vY[k];
        }

        // calculate angles between consecutive pairs
        for(int k=0; k<NDET-1; k++){
            int npairs = min(x[k].size(), x[k+1].size());

            for(int i=0; i<npairs; i++){
                
                // positions
                double x0 = x[k][i] - z_dist*tan(5*TMath::Pi()/4); // rotation of 225º 
                double y0 = y[k][i];
                double z0 = +z_dist;
                double x1 = x[k+1][i] + z_dist*tan(5*TMath::Pi()/4);
                double y1 = y[k+1][i];
                double z1 = -z_dist;

                double Vx = x1 - x0;
                double Vy = y1 - y0;
                double Vz = z1 - z0;

                double Vmod = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

                // cos(theta) with W=(1,0,-1) vector along the beam axis
                if(Vmod > 1e-9){
                    double cosTh = (Vx - Vz) / (Vmod * sqrt(2.0));
                    costheta[k].push_back(cosTh);
                } else {
                    costheta[k].push_back(-9999);
                }

                // phi angle calculation
                double x0p,y0p,z0p, x1p,y1p,z1p;
                rotY(x0,y0,z0, 5*TMath::Pi()/4, x0p,y0p,z0p);
                rotY(x1,y1,z1, 5*TMath::Pi()/4, x1p,y1p,z1p);
                double dXp = x1p - x0p;
                double dYp = y1p - y0p;
                phi[k].push_back( atan2(dYp, dXp) );
            }
        }

        tout->Fill();
    }

    fout->Write();
    fout->Close();
    fin->Close();

    cout<<"File "<<outfile<<" written \n";
}


