// anode_cathode_positions.C
// Uso: root -l -q 'anode_cathode_positions.C("input.root","T","output.root")'

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TVector3.h>
#include <TMath.h>
#include <vector>
#include <iostream>
using namespace std;

const int NDET   = 10;
const int MAXA   = 5;    // max multiplicidad para anodos
const int MAXC   = 20;   // max multiplicidad para cátodos
const double L_mm   = 200.0; // longitud delay line
const double z_dist = 25.0;  // distancia target–detector

// Rotación alrededor del eje Y (ángulo en rad)
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

    // --- Bind inputs ---
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

    // --- Histogramas para calibración ---
    vector<TH1D*> hX(NDET), hY(NDET);
    for(int k=0;k<NDET;k++){
        hX[k] = new TH1D(Form("hX%d",k),"",2000,0,1000);
        hY[k] = new TH1D(Form("hY%d",k),"",2000,0,1000);
    }
    for(Long64_t ev=0; ev<nentries; ev++){
        tin->GetEntry(ev);
        for(int k=0;k<NDET;k++){
            for(int i=0;i<mult[k];i++){
                hX[k]->Fill(tofC4[k][i] + tofC3[k][i] - 2*tofA[k][i]);
                hY[k]->Fill(tofC2[k][i] + tofC1[k][i] - 2*tofA[k][i]);
            }
        }
    }

    // --- Calibración de velocidades ---
    vector<double> vX(NDET,0), vY(NDET,0);
    for(int k=0;k<NDET;k++){
        TSpectrum sx(5); sx.Search(hX[k],2,"",0.1);
        TSpectrum sy(5); sy.Search(hY[k],2,"",0.1);
        double px = (sx.GetNPeaks()>0 ? sx.GetPositionX()[0] : hX[k]->GetMean());
        double py = (sy.GetNPeaks()>0 ? sy.GetPositionX()[0] : hY[k]->GetMean());
        if(px>1e-6) vX[k] = L_mm/px;
        if(py>1e-6) vY[k] = L_mm/py;
        cout<<"det "<<k<<" vX="<<vX[k]<<" vY="<<vY[k]<<"\n";
    }

    // --- Output ---
    TFile *fout = TFile::Open(outfile,"RECREATE");
    TTree *tout = tin->CloneTree(0); // copia TODOS los branches originales

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

    // --- Loop de eventos ---
    for(Long64_t ev=0; ev<nentries; ev++){
        tin->GetEntry(ev);

        // limpiar vectores al inicio de cada evento
        for(int k=0;k<NDET;k++){ x[k].clear(); y[k].clear(); }
        for(int k=0;k<NDET-1;k++){ costheta[k].clear(); phi[k].clear(); }

        // reconstruir posiciones x,y
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

        // calcular ángulos entre pares consecutivos
        for(int k=0; k<NDET-1; k++){
            int npairs = min(x[k].size(), x[k+1].size());

            for(int i=0; i<npairs; i++){
                
                // posiciones en 3D (z0 y z1)
                double x0 = x[k][i] - z_dist*tan(5*TMath::Pi()/4); // rotación +45° alrededor de Y
                double y0 = y[k][i];
                double z0 = +z_dist;
                double x1 = x[k+1][i] + z_dist*tan(5*TMath::Pi()/4); // rotación +45° alrededor de Y
                double y1 = y[k+1][i];
                double z1 = -z_dist;

                double Vx = x1 - x0;
                double Vy = y1 - y0;
                double Vz = z1 - z0;

                double Vmod = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

                // cos(theta) con W=(1,0,-1)
                if(Vmod > 1e-9){
                    double cosTh = (Vx - Vz) / (Vmod * TMath::Sqrt2());
                    costheta[k].push_back(cosTh);
                } else {
                    costheta[k].push_back(-9999);
                }

                // phi tras rotación +45° alrededor de Y
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

    cout<<"Archivo "<<outfile<<" escrito con TODOS los branches originales + posiciones, ángulos y velocidades\n";
}


