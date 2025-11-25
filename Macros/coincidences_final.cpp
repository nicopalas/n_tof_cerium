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

const int NDET = 10;     // total de detectores en los datos
const int MAXA = 50;     // máximo número de anodos por detector (entrada)
const int MAXC = 100;    // máximo número de catodos por tipo (entrada)
const int MAXOUT = 400;  // máximo número de catodos almacenados en el árbol de salida
const double L_mm = 100.0;
const double SUM_MIN = 80.0;
const double SUM_MAX = 120.0;

Int_t repeated_c; // number of repeated combinations for the same anode (change in cathodes)
Int_t repeated_a; // number of repeated combinations for the same cathode (change in anodes)

Int_t event = -1;
Int_t mult = -1;
// Time-of-flight
double tofA0[10], tofA1[10], tofA2[10], tofA3[10], tofA4[10];
double tofA5[10], tofA6[10], tofA7[10], tofA8[10], tofA9[10];

double dtA;

// Positions detector
double x0[10], y0v[10], x1[10], y1v[10], x2[10], y2[10], x3[10], y3[10];
double x4[10], y4[10], x5[10], y5[10], x6[10], y6[10], x7[10], y7[10];
double x8[10], y8[10], x9[10], y9[10];

// Amplitudes Anodes
float ampA0[10], ampA1[10], ampA2[10], ampA3[10], ampA4[10];
float ampA5[10], ampA6[10], ampA7[10], ampA8[10], ampA9[10];

// Amplitudes Cathodes
float ampC01[10], ampC02[10], ampC03[10], ampC04[10];
float ampC11[10], ampC12[10], ampC13[10], ampC14[10];
float ampC21[10], ampC22[10], ampC23[10], ampC24[10];
float ampC31[10], ampC32[10], ampC33[10], ampC34[10];
float ampC41[10], ampC42[10], ampC43[10], ampC44[10];
float ampC51[10], ampC52[10], ampC53[10], ampC54[10];
float ampC61[10], ampC62[10], ampC63[10], ampC64[10];
float ampC71[10], ampC72[10], ampC73[10], ampC74[10];
float ampC81[10], ampC82[10], ampC83[10], ampC84[10];
float ampC91[10], ampC92[10], ampC93[10], ampC94[10];

// Sums
double sumX0[10], sumY0[10], sumX1[10], sumY1[10], sumX2[10], sumY2[10];
double sumX3[10], sumY3[10], sumX4[10], sumY4[10], sumX5[10], sumY5[10];
double sumX6[10], sumY6[10], sumX7[10], sumY7[10], sumX8[10], sumY8[10];
double sumX9[10], sumY9[10];

// Differences
double diffX0[10], diffY0[10], diffX1[10], diffY1[10], diffX2[10], diffY2[10];
double diffX3[10], diffY3[10], diffX4[10], diffY4[10], diffX5[10], diffY5[10];
double diffX6[10], diffY6[10], diffX7[10], diffY7[10], diffX8[10], diffY8[10];
double diffX9[10], diffY9[10];

// Ratios
float ratioX0[10], ratioY0[10], ratioX1[10], ratioY1[10], ratioX2[10], ratioY2[10];
float ratioX3[10], ratioY3[10], ratioX4[10], ratioY4[10], ratioX5[10], ratioY5[10];
float ratioX6[10], ratioY6[10], ratioX7[10], ratioY7[10], ratioX8[10], ratioY8[10];
float ratioX9[10], ratioY9[10];

// Sum of amplitudes
float sumX0_ampA0[10], sumY0_ampA0[10], sumX1_ampA1[10], sumY1_ampA1[10];
float sumX2_ampA2[10], sumY2_ampA2[10], sumX3_ampA3[10], sumY3_ampA3[10];
float sumX4_ampA4[10], sumY4_ampA4[10], sumX5_ampA5[10], sumY5_ampA5[10];
float sumX6_ampA6[10], sumY6_ampA6[10], sumX7_ampA7[10], sumY7_ampA7[10];
float sumX8_ampA8[10], sumY8_ampA8[10], sumX9_ampA9[10], sumY9_ampA9[10];



struct hit { 
    int det; double tof_A0; float amp_A0; double tofC1_A0, tofC2_A0, tofC3_A0, tofC4_A0; 
    float ampC1_A0, ampC2_A0, ampC3_A0, ampC4_A0; 
    double sumX_0, sumY_0, diffX_0, diffY_0, ratioX_0, ratioY_0,x_0, y_0; 
    float sumX_amp_0, sumY_amp_0; 
    float ratioX_amp_0, ratioY_amp_0; };



// ------- printHit() -------
void printHit(const hit& h)
{
    std::cout << "Hit -- det = " << h.det << "\n";
    std::cout << "  tof_A0        = " << h.tof_A0 << "\n";
    std::cout << "  amp_A0        = " << h.amp_A0 << "\n";

    std::cout << "  tofC1_A0      = " << h.tofC1_A0 << "\n";
    std::cout << "  tofC2_A0      = " << h.tofC2_A0 << "\n";
    std::cout << "  tofC3_A0      = " << h.tofC3_A0 << "\n";
    std::cout << "  tofC4_A0      = " << h.tofC4_A0 << "\n";

    std::cout << "  ampC1_A0      = " << h.ampC1_A0 << "\n";
    std::cout << "  ampC2_A0      = " << h.ampC2_A0 << "\n";
    std::cout << "  ampC3_A0      = " << h.ampC3_A0 << "\n";
    std::cout << "  ampC4_A0      = " << h.ampC4_A0 << "\n";

    std::cout << "  sumX_0        = " << h.sumX_0 << "\n";
    std::cout << "  sumY_0        = " << h.sumY_0 << "\n";
    std::cout << "  diffX_0       = " << h.diffX_0 << "\n";
    std::cout << "  diffY_0       = " << h.diffY_0 << "\n";
    std::cout << "  ratioX_0      = " << h.ratioX_0 << "\n";
    std::cout << "  ratioY_0      = " << h.ratioY_0 << "\n";

    std::cout << "  x_0           = " << h.x_0 << "\n";
    std::cout << "  y_0           = " << h.y_0 << "\n";

    std::cout << "  sumX_amp_0    = " << h.sumX_amp_0 << "\n";
    std::cout << "  sumY_amp_0    = " << h.sumY_amp_0 << "\n";
    std::cout << "  ratioX_amp_0  = " << h.ratioX_amp_0 << "\n";
    std::cout << "  ratioY_amp_0  = " << h.ratioY_amp_0 << "\n";

    std::cout << "---------------------------------------------\n";
}






void positions(const char *infile = "in.root",
               const char *treename = "T",
               const char *outfile = "positions.root")
{

    TFile *fin = TFile::Open(infile, "READ");
    if(!fin || fin->IsZombie()){ std::cerr<<"Error abriendo "<<infile<<"\n"; return; }
    TTree *tin = (TTree*)fin->Get(treename);
    if(!tin){ std::cerr<<"Tree "<<treename<<" no encontrado en "<<infile<<"\n"; fin->Close(); return; }
    Long64_t nentries = tin->GetEntries();
    std::cout<<"Input: "<<infile<<"  entries="<<nentries<<"\n";

    //detectors used
    const int NDET_USED = 10;
    int detList[NDET_USED] = {0,1,2,3,4,5,6,7,8,9};

    //set variables
    static Double_t tofA[NDET][MAXA];
    static Float_t  ampA[NDET][MAXA];
    static Double_t tofC1[NDET][MAXC], tofC2[NDET][MAXC], tofC3[NDET][MAXC], tofC4[NDET][MAXC];
    static Float_t  ampC1[NDET][MAXC], ampC2[NDET][MAXC], ampC3[NDET][MAXC], ampC4[NDET][MAXC];
    static Int_t multA[NDET], multC1[NDET], multC2[NDET], multC3[NDET], multC4[NDET];

    // branches
    for(int detIdx=0; detIdx<NDET_USED; ++detIdx){
        int k = detList[detIdx];
        if(!tin->GetBranch(Form("tof%d",k))) std::cerr<<"Warning: branch tof"<<k<<" no encontrada\n";
        else tin->SetBranchAddress(Form("tof%d",k), tofA[k]);

        if(!tin->GetBranch(Form("amp%d",k))) std::cerr<<"Warning: branch amp"<<k<<" no encontrada\n";
        else tin->SetBranchAddress(Form("amp%d",k), ampA[k]);

        if(!tin->GetBranch(Form("mult%d",k))) std::cerr<<"Warning: branch mult"<<k<<" no encontrada\n";
        else tin->SetBranchAddress(Form("mult%d",k), &multA[k]);

        if(!tin->GetBranch(Form("tof%d1",k))) std::cerr<<"Warning: branch tof"<<k<<"1 no encontrada\n";
        else tin->SetBranchAddress(Form("tof%d1",k), tofC1[k]);

        if(!tin->GetBranch(Form("tof%d2",k))) std::cerr<<"Warning: branch tof"<<k<<"2 no encontrada\n";
        else tin->SetBranchAddress(Form("tof%d2",k), tofC2[k]);

        if(!tin->GetBranch(Form("tof%d3",k))) std::cerr<<"Warning: branch tof"<<k<<"3 no encontrada\n";
        else tin->SetBranchAddress(Form("tof%d3",k), tofC3[k]);

        if(!tin->GetBranch(Form("tof%d4",k))) std::cerr<<"Warning: branch tof"<<k<<"4 no encontrada\n";
        else tin->SetBranchAddress(Form("tof%d4",k), tofC4[k]);

        if(!tin->GetBranch(Form("mult%d1",k))) std::cerr<<"Warning: branch mult"<<k<<"1 no encontrada\n";
        else tin->SetBranchAddress(Form("mult%d1",k), &multC1[k]);

        if(!tin->GetBranch(Form("mult%d2",k))) std::cerr<<"Warning: branch mult"<<k<<"2 no encontrada\n";
        else tin->SetBranchAddress(Form("mult%d2",k), &multC2[k]);

        if(!tin->GetBranch(Form("mult%d3",k))) std::cerr<<"Warning: branch mult"<<k<<"3 no encontrada\n";
        else tin->SetBranchAddress(Form("mult%d3",k), &multC3[k]);

        if(!tin->GetBranch(Form("mult%d4",k))) std::cerr<<"Warning: branch mult"<<k<<"4 no encontrada\n";
        else tin->SetBranchAddress(Form("mult%d4",k), &multC4[k]);

        if(!tin->GetBranch(Form("amp%d1",k))) std::cerr<<"Warning: branch amp"<<k<<"1 no encontrada\n";
        else tin->SetBranchAddress(Form("amp%d1",k), ampC1[k]);

        if(!tin->GetBranch(Form("amp%d2",k))) std::cerr<<"Warning: branch amp"<<k<<"2 no encontrada\n";
        else tin->SetBranchAddress(Form("amp%d2",k), ampC2[k]);

        if(!tin->GetBranch(Form("amp%d3",k))) std::cerr<<"Warning: branch amp"<<k<<"3 no encontrada\n";
        else tin->SetBranchAddress(Form("amp%d3",k), ampC3[k]);

        if(!tin->GetBranch(Form("amp%d4",k))) std::cerr<<"Warning: branch amp"<<k<<"4 no encontrada\n";
        else tin->SetBranchAddress(Form("amp%d4",k), ampC4[k]);
    }

    // Histograms per detector in detList
    std::vector<TH1D*> hSumX(NDET_USED), hSumY(NDET_USED), hDiffX(NDET_USED), hDiffY(NDET_USED);
    for(int idx=0; idx<NDET_USED; idx++){
        int k = detList[idx];
        hSumX[idx]  = new TH1D(Form("hSumX_%d",k), Form("Det %d sumX",k), 200, 0, 200);
        hSumY[idx]  = new TH1D(Form("hSumY_%d",k), Form("Det %d sumY",k), 200, 0, 200);
        hDiffX[idx] = new TH1D(Form("hDiffX_%d",k), Form("Det %d diffX",k), 400, -200, 200);
        hDiffY[idx] = new TH1D(Form("hDiffY_%d",k), Form("Det %d diffY",k), 400, -200, 200);
    }

    Long64_t printLimit = std::min<Long64_t>(10, nentries);
    for(Long64_t ev=0; ev<1000000; ++ev){
        tin->GetEntry(ev);

        // error handling if multiplicities are negative or too large
        for(int id=0; id<NDET_USED; ++id){
            int k = detList[id];
            if(multA[k] < 0) multA[k] = 0;
            if(multC1[k] < 0) multC1[k] = 0;
            if(multC2[k] < 0) multC2[k] = 0;
            if(multC3[k] < 0) multC3[k] = 0;
            if(multC4[k] < 0) multC4[k] = 0;

            if(multA[k] > MAXA){ std::cerr<<"Warning: multA["<<k<<"]="<<multA[k]<<" > MAXA -> truncated\n"; multA[k] = MAXA; }
            if(multC1[k] > MAXC){ std::cerr<<"Warning: multC1["<<k<<"]="<<multC1[k]<<" > MAXC -> truncated\n"; multC1[k] = MAXC; }
            if(multC2[k] > MAXC){std::cerr<<"Warning: multC2["<<k<<"]="<<multC2[k]<<" > MAXC -> truncated\n"; multC2[k] = MAXC; }
            if(multC3[k] > MAXC){ std::cerr<<"Warning: multC3["<<k<<"]="<<multC3[k]<<" > MAXC -> truncated\n"; multC3[k] = MAXC; }
            if(multC4[k] > MAXC){  std::cerr<<"Warning: multC4["<<k<<"]="<<multC4[k]<<" > MAXC -> truncated\n"; multC4[k] = MAXC; }
        }

        for(int idx=0; idx<NDET_USED; ++idx){
            int k = detList[idx];
            if(multA[k] < 1) continue;
            for(int i=0;i<multA[k]; ++i){
                for(int j1=0;j1<multC1[k]; ++j1){
                    for(int j2=0;j2<multC2[k]; ++j2){
                        double sumX = tofC2[k][j2] + tofC1[k][j1] - 2*tofA[k][i];
                        double diffX = tofC2[k][j2] - tofC1[k][j1];
                        hSumX[idx]->Fill(sumX);
                        if(sumX > SUM_MIN && sumX < SUM_MAX && fabs(diffX) < (L_mm / 0.8)){ 
                        //we will only fill diff if within physical limits (ad hoc)
                        hDiffX[idx]->Fill(diffX);
                        }
                    }
                }
                for(int j3=0;j3<multC3[k]; ++j3){
                    for(int j4=0;j4<multC4[k]; ++j4){
                        double sumY = tofC4[k][j4] + tofC3[k][j3] - 2*tofA[k][i];
                        double diffY = tofC4[k][j4] - tofC3[k][j3];
                        hSumY[idx]->Fill(sumY);
                        if(sumY > SUM_MIN && sumY < SUM_MAX && fabs(diffY) < (L_mm / 0.8)){
                            hDiffY[idx]->Fill(diffY);
                        }
                    }
                }
            }
        }

        if(ev < printLimit){
            // print multiplicities of the first events for debug
            int k0 = detList[0];
            std::cout << "Event " << ev
                      << " multA["<<k0<<"]="<<multA[k0]
                      << " multC1="<<multC1[k0]<<" multC2="<<multC2[k0]
                      << " multC3="<<multC3[k0]<<" multC4="<<multC4[k0]
                      << std::endl;
        }
    }

    // Calibration: find vX, vY and means/sigmas of sums
    std::vector<double> vX(NDET_USED,0), vY(NDET_USED,0);
    std::vector<double> meanSumX(NDET_USED,100), meanSumY(NDET_USED,100);
    std::vector<double> meanWidthX(NDET_USED,2), meanWidthY(NDET_USED,2);

    for(int idx=0; idx<NDET_USED; ++idx){
        int k = detList[idx];
        TSpectrum sx(20), sy(20);
        sx.Search(hDiffX[idx], 2, "nodraw", 0.05);
        sy.Search(hDiffY[idx], 2, "nodraw", 0.05);

        std::vector<double> px, py;
        for(int i=0;i<sx.GetNPeaks(); ++i){
            double p = sx.GetPositionX()[i];
            if(fabs(p-110)<5 || fabs(p+110)<5) px.push_back(p);
        }
        for(int i=0;i<sy.GetNPeaks(); ++i){
            double p = sy.GetPositionX()[i];
            if(fabs(p-110)<5 || fabs(p+110)<5) py.push_back(p);
        }
        std::sort(px.begin(),px.end());
        std::sort(py.begin(),py.end());
        if(px.size()>=2) vX[idx] = L_mm / (0.5*(fabs(px.front())+fabs(px.back())));
        if(py.size()>=2) vY[idx] = L_mm / (0.5*(fabs(py.front())+fabs(py.back())));
        std::cout << "Detector " << k
                  << "  vX = " << vX[idx]
                  << "  vY = " << vY[idx]
                  << std::endl;

        // fitting
        if (hSumX[idx]->GetEntries() > 50) {
            double meanHist = hSumX[idx]->GetMean();
            double peak = hSumX[idx]->GetMaximum();
            // gaus(0)+pol1(3): params [0]=A, [1]=mean, [2]=sigma, [3]=a0, [4]=a1
            TF1 *fX = new TF1(Form("fX_%d", k), "gaus(0) + pol1(3)", SUM_MIN, SUM_MAX);
            fX->SetParNames("A", "mean", "sigma", "a0", "a1");
            fX->SetParameter(0, peak);
            fX->SetParameter(1, 100.0);
            fX->SetParameter(2, 2.0);
            // background initial guess. the mean content of the histogram as plateau and 0 as the slope
            double meanContent = (hSumX[idx]->GetNbinsX()>0) ? (hSumX[idx]->Integral() / hSumX[idx]->GetNbinsX()) : 0.0;
            fX->SetParameter(3, meanContent);
            fX->SetParameter(4, 0.0);

            // 
            fX->SetParLimits(1, 90.0, 110.0);  // center around 100
            fX->SetParLimits(2, 0.5, 5.0);     // expected width (above its resolution)
            fX->SetParLimits(3, 0.0, 10.0 * std::max(1.0, meanContent)); // background
            fX->SetParLimits(4, -3.1, 3.1);    // slope

            hSumX[idx]->Fit(fX, "RQ"); // R: use range, Q: quiet

            meanSumX[idx]   = fX->GetParameter(1);
            meanWidthX[idx] = fX->GetParameter(2);

            std::cout << "Detector " << k
                      << "  meanSumX = " << meanSumX[idx]
                      << "  sigma = " << meanWidthX[idx]
                      << "  (chi2/ndf=" << (fX->GetNDF()? fX->GetChisquare()/fX->GetNDF() : -1) << ")"
                      << std::endl;
        }

        // fitting-
        if (hSumY[idx]->GetEntries() > 50) {
            double meanHist = hSumY[idx]->GetMean();
            double peak = hSumY[idx]->GetMaximum();
            TF1 *fY = new TF1(Form("fY_%d", k), "gaus(0) + pol1(3)", SUM_MIN, SUM_MAX);
            fY->SetParNames("A", "mean", "sigma", "a0", "a1");
            fY->SetParameter(0, peak);
            fY->SetParameter(1, 100.0);
            fY->SetParameter(2, 2.0);
            double meanContent = (hSumY[idx]->GetNbinsX()>0) ? (hSumY[idx]->Integral() / hSumY[idx]->GetNbinsX()) : 0.0;
            fY->SetParameter(3, meanContent);
            fY->SetParameter(4, 0.0);

            fY->SetParLimits(1, 90.0, 110.0);
            fY->SetParLimits(2, 0.5, 5.0);
            fY->SetParLimits(3, 0.0, 10.0 * std::max(1.0, meanContent));
            fY->SetParLimits(4, -3.1, 3.1);

            hSumY[idx]->Fit(fY, "RQ");

            meanSumY[idx]   = fY->GetParameter(1);
            meanWidthY[idx] = fY->GetParameter(2);

            std::cout << "Detector " << k
                      << "  meanSumY = " << meanSumY[idx]
                      << "  sigma = " << meanWidthY[idx]
                      << "  (chi2/ndf=" << (fY->GetNDF()? fY->GetChisquare()/fY->GetNDF() : -1) << ")"
                      << std::endl;
        }
    } // end calibration loop
    TFile *fout = new TFile(outfile,"RECREATE");
    TTree *tbdt = new TTree("coincidences","coincidences_tree");
    int mult0, mult1, mult2, mult3, mult4, mult5, mult6, mult7, mult8, mult9;

    // mult branches
tbdt->Branch("mult0", &mult0, "mult0/I");
tbdt->Branch("mult1", &mult1, "mult1/I");
tbdt->Branch("mult2", &mult2, "mult2/I");
tbdt->Branch("mult3", &mult3, "mult3/I");
tbdt->Branch("mult4", &mult4, "mult4/I");
tbdt->Branch("mult5", &mult5, "mult5/I");
tbdt->Branch("mult6", &mult6, "mult6/I");
tbdt->Branch("mult7", &mult7, "mult7/I");
tbdt->Branch("mult8", &mult8, "mult8/I");
tbdt->Branch("mult9", &mult9, "mult9/I");
tbdt->Branch("mult", &mult, "mult/I");
tbdt->Branch("event", &event, "event/I");

// TOF
tbdt->Branch("tofA0", tofA0, "tofA0[mult0]/D");
tbdt->Branch("tofA1", tofA1, "tofA1[mult1]/D");
tbdt->Branch("tofA2", tofA2, "tofA2[mult2]/D");
tbdt->Branch("tofA3", tofA3, "tofA3[mult3]/D");
tbdt->Branch("tofA4", tofA4, "tofA4[mult4]/D");
tbdt->Branch("tofA5", tofA5, "tofA5[mult5]/D");
tbdt->Branch("tofA6", tofA6, "tofA6[mult6]/D");
tbdt->Branch("tofA7", tofA7, "tofA7[mult7]/D");
tbdt->Branch("tofA8", tofA8, "tofA8[mult8]/D");
tbdt->Branch("tofA9", tofA9, "tofA9[mult9]/D");

// detector positions
tbdt->Branch("x0", x0, "x0[mult0]/D"); tbdt->Branch("y0", y0v, "y0[mult0]/D");
tbdt->Branch("x1", x1, "x1[mult1]/D"); tbdt->Branch("y1", y1v, "y1[mult1]/D");
tbdt->Branch("x2", x2, "x2[mult2]/D"); tbdt->Branch("y2", y2, "y2[mult2]/D");
tbdt->Branch("x3", x3, "x3[mult3]/D"); tbdt->Branch("y3", y3, "y3[mult3]/D");
tbdt->Branch("x4", x4, "x4[mult4]/D"); tbdt->Branch("y4", y4, "y4[mult4]/D");
tbdt->Branch("x5", x5, "x5[mult5]/D"); tbdt->Branch("y5", y5, "y5[mult5]/D");
tbdt->Branch("x6", x6, "x6[mult6]/D"); tbdt->Branch("y6", y6, "y6[mult6]/D");
tbdt->Branch("x7", x7, "x7[mult7]/D"); tbdt->Branch("y7", y7, "y7[mult7]/D");
tbdt->Branch("x8", x8, "x8[mult8]/D"); tbdt->Branch("y8", y8, "y8[mult8]/D");
tbdt->Branch("x9", x9, "x9[mult9]/D"); tbdt->Branch("y9", y9, "y9[mult9]/D");

// Amplitudes
tbdt->Branch("ampA0", ampA0, "ampA0[mult0]/F");
tbdt->Branch("ampA1", ampA1, "ampA1[mult1]/F");
tbdt->Branch("ampA2", ampA2, "ampA2[mult2]/F");
tbdt->Branch("ampA3", ampA3, "ampA3[mult3]/F");
tbdt->Branch("ampA4", ampA4, "ampA4[mult4]/F");
tbdt->Branch("ampA5", ampA5, "ampA5[mult5]/F");
tbdt->Branch("ampA6", ampA6, "ampA6[mult6]/F");
tbdt->Branch("ampA7", ampA7, "ampA7[mult7]/F");
tbdt->Branch("ampA8", ampA8, "ampA8[mult8]/F");
tbdt->Branch("ampA9", ampA9, "ampA9[mult9]/F");

// Amplitudes cathode
tbdt->Branch("ampC01", ampC01, "ampC01[mult0]/F"); tbdt->Branch("ampC02", ampC02, "ampC02[mult0]/F");
tbdt->Branch("ampC03", ampC03, "ampC03[mult0]/F"); tbdt->Branch("ampC04", ampC04, "ampC04[mult0]/F");

tbdt->Branch("ampC11", ampC11, "ampC11[mult1]/F"); tbdt->Branch("ampC12", ampC12, "ampC12[mult1]/F");
tbdt->Branch("ampC13", ampC13, "ampC13[mult1]/F"); tbdt->Branch("ampC14", ampC14, "ampC14[mult1]/F");

tbdt->Branch("ampC21", ampC21, "ampC21[mult2]/F"); tbdt->Branch("ampC22", ampC22, "ampC22[mult2]/F");
tbdt->Branch("ampC23", ampC23, "ampC23[mult2]/F"); tbdt->Branch("ampC24", ampC24, "ampC24[mult2]/F");

tbdt->Branch("ampC31", ampC31, "ampC31[mult3]/F"); tbdt->Branch("ampC32", ampC32, "ampC32[mult3]/F");
tbdt->Branch("ampC33", ampC33, "ampC33[mult3]/F"); tbdt->Branch("ampC34", ampC34, "ampC34[mult3]/F");

tbdt->Branch("ampC41", ampC41, "ampC41[mult4]/F"); tbdt->Branch("ampC42", ampC42, "ampC42[mult4]/F");
tbdt->Branch("ampC43", ampC43, "ampC43[mult4]/F"); tbdt->Branch("ampC44", ampC44, "ampC44[mult4]/F");

tbdt->Branch("ampC51", ampC51, "ampC51[mult5]/F"); tbdt->Branch("ampC52", ampC52, "ampC52[mult5]/F");
tbdt->Branch("ampC53", ampC53, "ampC53[mult5]/F"); tbdt->Branch("ampC54", ampC54, "ampC54[mult5]/F");

tbdt->Branch("ampC61", ampC61, "ampC61[mult6]/F"); tbdt->Branch("ampC62", ampC62, "ampC62[mult6]/F");
tbdt->Branch("ampC63", ampC63, "ampC63[mult6]/F"); tbdt->Branch("ampC64", ampC64, "ampC64[mult6]/F");

tbdt->Branch("ampC71", ampC71, "ampC71[mult7]/F"); tbdt->Branch("ampC72", ampC72, "ampC72[mult7]/F");
tbdt->Branch("ampC73", ampC73, "ampC73[mult7]/F"); tbdt->Branch("ampC74", ampC74, "ampC74[mult7]/F");

tbdt->Branch("ampC81", ampC81, "ampC81[mult8]/F"); tbdt->Branch("ampC82", ampC82, "ampC82[mult8]/F");
tbdt->Branch("ampC83", ampC83, "ampC83[mult8]/F"); tbdt->Branch("ampC84", ampC84, "ampC84[mult8]/F");

tbdt->Branch("ampC91", ampC91, "ampC91[mult9]/F"); tbdt->Branch("ampC92", ampC92, "ampC92[mult9]/F");
tbdt->Branch("ampC93", ampC93, "ampC93[mult9]/F"); tbdt->Branch("ampC94", ampC94, "ampC94[mult9]/F");

// sumX, sumY
tbdt->Branch("sumX0", sumX0, "sumX0[mult0]/D"); tbdt->Branch("sumY0", sumY0, "sumY0[mult0]/D");
tbdt->Branch("sumX1", sumX1, "sumX1[mult1]/D"); tbdt->Branch("sumY1", sumY1, "sumY1[mult1]/D");
tbdt->Branch("sumX2", sumX2, "sumX2[mult2]/D"); tbdt->Branch("sumY2", sumY2, "sumY2[mult2]/D");
tbdt->Branch("sumX3", sumX3, "sumX3[mult3]/D"); tbdt->Branch("sumY3", sumY3, "sumY3[mult3]/D");
tbdt->Branch("sumX4", sumX4, "sumX4[mult4]/D"); tbdt->Branch("sumY4", sumY4, "sumY4[mult4]/D");
tbdt->Branch("sumX5", sumX5, "sumX5[mult5]/D"); tbdt->Branch("sumY5", sumY5, "sumY5[mult5]/D");
tbdt->Branch("sumX6", sumX6, "sumX6[mult6]/D"); tbdt->Branch("sumY6", sumY6, "sumY6[mult6]/D");
tbdt->Branch("sumX7", sumX7, "sumX7[mult7]/D"); tbdt->Branch("sumY7", sumY7, "sumY7[mult7]/D");
tbdt->Branch("sumX8", sumX8, "sumX8[mult8]/D"); tbdt->Branch("sumY8", sumY8, "sumY8[mult8]/D");
tbdt->Branch("sumX9", sumX9, "sumX9[mult9]/D"); tbdt->Branch("sumY9", sumY9, "sumY9[mult9]/D");

// diffX, diffY
tbdt->Branch("diffX0", diffX0, "diffX0[mult0]/D"); tbdt->Branch("diffY0", diffY0, "diffY0[mult0]/D");
tbdt->Branch("diffX1", diffX1, "diffX1[mult1]/D"); tbdt->Branch("diffY1", diffY1, "diffY1[mult1]/D");
tbdt->Branch("diffX2", diffX2, "diffX2[mult2]/D"); tbdt->Branch("diffY2", diffY2, "diffY2[mult2]/D");
tbdt->Branch("diffX3", diffX3, "diffX3[mult3]/D"); tbdt->Branch("diffY3", diffY3, "diffY3[mult3]/D");
tbdt->Branch("diffX4", diffX4, "diffX4[mult4]/D"); tbdt->Branch("diffY4", diffY4, "diffY4[mult4]/D");
tbdt->Branch("diffX5", diffX5, "diffX5[mult5]/D"); tbdt->Branch("diffY5", diffY5, "diffY5[mult5]/D");
tbdt->Branch("diffX6", diffX6, "diffX6[mult6]/D"); tbdt->Branch("diffY6", diffY6, "diffY6[mult6]/D");
tbdt->Branch("diffX7", diffX7, "diffX7[mult7]/D"); tbdt->Branch("diffY7", diffY7, "diffY7[mult7]/D");
tbdt->Branch("diffX8", diffX8, "diffX8[mult8]/D"); tbdt->Branch("diffY8", diffY8, "diffY8[mult8]/D");
tbdt->Branch("diffX9", diffX9, "diffX9[mult9]/D"); tbdt->Branch("diffY9", diffY9, "diffY9[mult9]/D");

// ratioX, ratioY
tbdt->Branch("ratioX0", ratioX0, "ratioX0[mult0]/F"); tbdt->Branch("ratioY0", ratioY0, "ratioY0[mult0]/F");
tbdt->Branch("ratioX1", ratioX1, "ratioX1[mult1]/F"); tbdt->Branch("ratioY1", ratioY1, "ratioY1[mult1]/F");
tbdt->Branch("ratioX2", ratioX2, "ratioX2[mult2]/F"); tbdt->Branch("ratioY2", ratioY2, "ratioY2[mult2]/F");
tbdt->Branch("ratioX3", ratioX3, "ratioX3[mult3]/F"); tbdt->Branch("ratioY3", ratioY3, "ratioY3[mult3]/F");
tbdt->Branch("ratioX4", ratioX4, "ratioX4[mult4]/F"); tbdt->Branch("ratioY4", ratioY4, "ratioY4[mult4]/F");
tbdt->Branch("ratioX5", ratioX5, "ratioX5[mult5]/F"); tbdt->Branch("ratioY5", ratioY5, "ratioY5[mult5]/F");
tbdt->Branch("ratioX6", ratioX6, "ratioX6[mult6]/F"); tbdt->Branch("ratioY6", ratioY6, "ratioY6[mult6]/F");
tbdt->Branch("ratioX7", ratioX7, "ratioX7[mult7]/F"); tbdt->Branch("ratioY7", ratioY7, "ratioY7[mult7]/F");
tbdt->Branch("ratioX8", ratioX8, "ratioX8[mult8]/F"); tbdt->Branch("ratioY8", ratioY8, "ratioY8[mult8]/F");
tbdt->Branch("ratioX9", ratioX9, "ratioX9[mult9]/F"); tbdt->Branch("ratioY9", ratioY9, "ratioY9[mult9]/F");

// sumX_ampA
tbdt->Branch("sumX0_ampA0", sumX0_ampA0, "sumX0_ampA0[mult0]/F"); tbdt->Branch("sumY0_ampA0", sumY0_ampA0, "sumY0_ampA0[mult0]/F");
tbdt->Branch("sumX1_ampA1", sumX1_ampA1, "sumX1_ampA1[mult1]/F"); tbdt->Branch("sumY1_ampA1", sumY1_ampA1, "sumY1_ampA1[mult1]/F");
tbdt->Branch("sumX2_ampA2", sumX2_ampA2, "sumX2_ampA2[mult2]/F"); tbdt->Branch("sumY2_ampA2", sumY2_ampA2, "sumY2_ampA2[mult2]/F");
tbdt->Branch("sumX3_ampA3", sumX3_ampA3, "sumX3_ampA3[mult3]/F"); tbdt->Branch("sumY3_ampA3", sumY3_ampA3, "sumY3_ampA3[mult3]/F");
tbdt->Branch("sumX4_ampA4", sumX4_ampA4, "sumX4_ampA4[mult4]/F"); tbdt->Branch("sumY4_ampA4", sumY4_ampA4, "sumY4_ampA4[mult4]/F");
tbdt->Branch("sumX5_ampA5", sumX5_ampA5, "sumX5_ampA5[mult5]/F"); tbdt->Branch("sumY5_ampA5", sumY5_ampA5, "sumY5_ampA5[mult5]/F");
tbdt->Branch("sumX6_ampA6", sumX6_ampA6, "sumX6_ampA6[mult6]/F"); tbdt->Branch("sumY6_ampA6", sumY6_ampA6, "sumY6_ampA6[mult6]/F");
tbdt->Branch("sumX7_ampA7", sumX7_ampA7, "sumX7_ampA7[mult7]/F"); tbdt->Branch("sumY7_ampA7", sumY7_ampA7, "sumY7_ampA7[mult7]/F");
tbdt->Branch("sumX8_ampA8", sumX8_ampA8, "sumX8_ampA8[mult8]/F"); tbdt->Branch("sumY8_ampA8", sumY8_ampA8, "sumY8_ampA8[mult8]/F");
tbdt->Branch("sumX9_ampA9", sumX9_ampA9, "sumX9_ampA9[mult9]/F"); tbdt->Branch("sumY9_ampA9", sumY9_ampA9, "sumY9_ampA9[mult9]/F");



std::vector<hit> hits;  // vector to store the hits for a certain event

Int_t events_good=0; // check how many events we are keeping
using HitKey = std::tuple<
    int,      // det
    double, float,       // tof_A0, amp_A0
    double, double,       // sumX_0, sumY_0
    double, double,       // diffX_0, diffY_0
    float, float,       // ratioX_0, ratioY_0
    float, float,       // sumX_amp_0, sumY_amp_0
    double, double        // x_0, y_0
>;

std::set<HitKey> usedHits;   // set to store unique hits. why are we using set? because it automatically returns true if the element is already present



for (Long64_t ev = 0; ev < nentries; ++ev) {
    if (ev%100000 == 0){
        std::cout << "ratio event_good: " << 100.0 * events_good / (ev + 1) << std::endl;
        std::cout << "% of events= " << 100.0*ev/(nentries) << std::endl;



    }

    tin->GetEntry(ev);

    mult0 = mult1 = mult2 = mult3 = mult4 =
    mult5 = mult6 = mult7 = mult8 = mult9 = 0;

    mult = 0;
    hits.clear();               // clear the vector and other variables for each event


    event = ev;

    for (int idx = 0; idx < NDET; ++idx) { // loop over all detectors
        if (multA[idx] < 1) continue;

        for (int i = 0; i < multA[idx]; i++) { // loop over all anode hits for a certain detector

            if (multC1[idx]<1 || multC2[idx]<1 || multC3[idx]<1 || multC4[idx]<1 ){
                continue;
            }
            // loop over all combinations of C1, C2, C3, C4 to search for valid cathodes for a certain anode
            for (int j1 = 0; j1 < multC1[idx]; j1++)
            for (int j2 = 0; j2 < multC2[idx]; j2++)
            for (int j3 = 0; j3 < multC3[idx]; j3++)
            for (int j4 = 0; j4 < multC4[idx]; j4++){
                    float ampC1_val = ampC1[idx][j1];
                    float ampC2_val = ampC2[idx][j2];
                    float ampC3_val = ampC3[idx][j3];
                    float ampC4_val = ampC4[idx][j4];
                    double sumX = tofC2[idx][j2] + tofC1[idx][j1] - 2 * tofA[idx][i];
                    double diffX = tofC2[idx][j2] - tofC1[idx][j1];
                    double ratioX = ampC2[idx][j2] / ampC1[idx][j1];
                    double ratio_amp_X = (ampC2[idx][j2] + ampC1[idx][j1]) / ampA[idx][i];

                    double sumY = tofC4[idx][j4] + tofC3[idx][j3] - 2 * tofA[idx][i];
                    double diffY = tofC4[idx][j4] - tofC3[idx][j3];
                    double ratioY = ampC4[idx][j4] / ampC3[idx][j3];
                    double ratio_amp_Y = (ampC4[idx][j4] + ampC3[idx][j3]) / ampA[idx][i];
                    double x_0 = vX[idx] * diffX / 2.0;
                    double y_0 = vY[idx] * diffY / 2.0;

                    // apply cuts
                   if (fabs(sumY - meanSumY[idx]) > 2 * meanWidthY[idx]
                    || fabs(diffY) > (L_mm / vY[idx]) || ratioY > 2.0 || ratioY < 0.5
                    || ratio_amp_Y < 0.2 || ratio_amp_Y > 1.8 || fabs(sumX - meanSumX[idx]) > 2 * meanWidthX[idx] 
                    || fabs(diffX) > (L_mm / vX[idx])  || ratioX > 2.0 || ratioX < 0.5  || ratio_amp_X < 0.2 || ratio_amp_X > 1.8) {
                                continue;
                    }
                    // if passed all cuts, create the hit
                    hit h;
                    h.det = idx;
                    h.tof_A0 = tofA[idx][i];
                    h.amp_A0 = ampA[idx][i];
                    h.tofC1_A0 = tofC1[idx][j1];
                    h.tofC2_A0 = tofC2[idx][j2];
                    h.tofC3_A0 = tofC3[idx][j3];
                    h.tofC4_A0 = tofC4[idx][j4];

                    h.ampC1_A0 = ampC1_val;
                    h.ampC2_A0 = ampC2_val;
                    h.ampC3_A0 = ampC3_val;
                    h.ampC4_A0 = ampC4_val;

                    h.sumX_0 = sumX;
                    h.sumY_0 = sumY;
                    h.diffX_0 = diffX;
                    h.diffY_0 = diffY;
                    h.ratioX_0 = ratioX;
                    h.ratioY_0 = ratioY;

                    h.x_0 = x_0;
                    h.y_0 = y_0;

                    h.sumX_amp_0 = ratio_amp_X;;
                    h.sumY_amp_0 = ratio_amp_Y;
                    HitKey key = std::make_tuple(h.det,h.tof_A0, h.amp_A0, h.sumX_0, h.sumY_0,
                    h.diffX_0, h.diffY_0, h.ratioX_0, h.ratioY_0,h.sumX_amp_0, h.sumY_amp_0, h.x_0, h.y_0);
                            if (usedHits.count(key) != 0) {
                                    continue;   // already exists the hit then skip
                    } 

                    usedHits.insert(key);     // mark as used
                    hits.push_back(h);        // adding
            }
        }
    }

    for (const auto &hit : hits) { // loop over all valid hits found in the event and fill the tree branches
        mult++;

        switch (hit.det) {
    case 0:
        tofA0[mult0] = hit.tof_A0;
        ampA0[mult0] = hit.amp_A0;
        sumX0[mult0] = hit.sumX_0;
        sumY0[mult0] = hit.sumY_0;
        diffX0[mult0] = hit.diffX_0;
        diffY0[mult0] = hit.diffY_0;
        ratioX0[mult0] = hit.ratioX_0;
        ratioY0[mult0] = hit.ratioY_0;
        sumX0_ampA0[mult0] = hit.sumX_amp_0;
        sumY0_ampA0[mult0] = hit.sumY_amp_0;
        x0[mult0] = hit.x_0;
        y0v[mult0] = hit.y_0;
        ampC01[mult0] = hit.ampC1_A0;
        ampC02[mult0] = hit.ampC2_A0;
        ampC03[mult0] = hit.ampC3_A0;
        ampC04[mult0] = hit.ampC4_A0;
        mult0++;
        break;

    case 1:
        tofA1[mult1] = hit.tof_A0;
        ampA1[mult1] = hit.amp_A0;
        sumX1[mult1] = hit.sumX_0;
        sumY1[mult1] = hit.sumY_0;
        diffX1[mult1] = hit.diffX_0;
        diffY1[mult1] = hit.diffY_0;
        ratioX1[mult1] = hit.ratioX_0;
        ratioY1[mult1] = hit.ratioY_0;
        sumX1_ampA1[mult1] = hit.sumX_amp_0;
        sumY1_ampA1[mult1] = hit.sumY_amp_0;
        ampC11[mult1] = hit.ampC1_A0;
        ampC12[mult1] = hit.ampC2_A0;
        ampC13[mult1] = hit.ampC3_A0;
        ampC14[mult1] = hit.ampC4_A0;
        x1[mult1] = hit.x_0;
        y1v[mult1] = hit.y_0;
        mult1++;
        break;

    case 2:
        tofA2[mult2] = hit.tof_A0;
        ampA2[mult2] = hit.amp_A0;
        sumX2[mult2] = hit.sumX_0;
        sumY2[mult2] = hit.sumY_0;
        diffX2[mult2] = hit.diffX_0;
        diffY2[mult2] = hit.diffY_0;
        ratioX2[mult2] = hit.ratioX_0;
        ratioY2[mult2] = hit.ratioY_0;
        sumX2_ampA2[mult2] = hit.sumX_amp_0;
        sumY2_ampA2[mult2] = hit.sumY_amp_0;
        ampC21[mult2] = hit.ampC1_A0;
        ampC22[mult2] = hit.ampC2_A0;
        ampC23[mult2] = hit.ampC3_A0;
        ampC24[mult2] = hit.ampC4_A0;
        x2[mult2] = hit.x_0;
        y2[mult2] = hit.y_0;
        mult2++;
        break;

    case 3:
        tofA3[mult3] = hit.tof_A0;
        ampA3[mult3] = hit.amp_A0;
        sumX3[mult3] = hit.sumX_0;
        sumY3[mult3] = hit.sumY_0;
        diffX3[mult3] = hit.diffX_0;
        diffY3[mult3] = hit.diffY_0;
        ratioX3[mult3] = hit.ratioX_0;
        ratioY3[mult3] = hit.ratioY_0;
        sumX3_ampA3[mult3] = hit.sumX_amp_0;
        sumY3_ampA3[mult3] = hit.sumY_amp_0;
        ampC31[mult3] = hit.ampC1_A0;
        ampC32[mult3] = hit.ampC2_A0;
        ampC33[mult3] = hit.ampC3_A0;
        ampC34[mult3] = hit.ampC4_A0;
        x3[mult3] = hit.x_0;
        y3[mult3] = hit.y_0;
        mult3++;
        break;

    case 4:
        tofA4[mult4] = hit.tof_A0;
        ampA4[mult4] = hit.amp_A0;
        sumX4[mult4] = hit.sumX_0;
        sumY4[mult4] = hit.sumY_0;
        diffX4[mult4] = hit.diffX_0;
        diffY4[mult4] = hit.diffY_0;
        ratioX4[mult4] = hit.ratioX_0;
        ratioY4[mult4] = hit.ratioY_0;
        sumX4_ampA4[mult4] = hit.sumX_amp_0;
        sumY4_ampA4[mult4] = hit.sumY_amp_0;
        ampC41[mult4] = hit.ampC1_A0;
        ampC42[mult4] = hit.ampC2_A0;
        ampC43[mult4] = hit.ampC3_A0;
        ampC44[mult4] = hit.ampC4_A0;
        x4[mult4] = hit.x_0;
        y4[mult4] = hit.y_0;
        mult4++;
        break;

    case 5:
        tofA5[mult5] = hit.tof_A0;
        ampA5[mult5] = hit.amp_A0;
        sumX5[mult5] = hit.sumX_0;
        sumY5[mult5] = hit.sumY_0;
        diffX5[mult5] = hit.diffX_0;
        diffY5[mult5] = hit.diffY_0;
        ratioX5[mult5] = hit.ratioX_0;
        ratioY5[mult5] = hit.ratioY_0;
        sumX5_ampA5[mult5] = hit.sumX_amp_0;
        sumY5_ampA5[mult5] = hit.sumY_amp_0;
        ampC51[mult5] = hit.ampC1_A0;
        ampC52[mult5] = hit.ampC2_A0;
        ampC53[mult5] = hit.ampC3_A0;
        ampC54[mult5] = hit.ampC4_A0;
        x5[mult5] = hit.x_0;
        y5[mult5] = hit.y_0;    
        mult5++;
        break;

    case 6:
        tofA6[mult6] = hit.tof_A0;
        ampA6[mult6] = hit.amp_A0;
        sumX6[mult6] = hit.sumX_0;
        sumY6[mult6] = hit.sumY_0;
        diffX6[mult6] = hit.diffX_0;
        diffY6[mult6] = hit.diffY_0;
        ratioX6[mult6] = hit.ratioX_0;
        ratioY6[mult6] = hit.ratioY_0;
        sumX6_ampA6[mult6] = hit.sumX_amp_0;
        sumY6_ampA6[mult6] = hit.sumY_amp_0;
        ampC61[mult6] = hit.ampC1_A0;
        ampC62[mult6] = hit.ampC2_A0;
        ampC63[mult6] = hit.ampC3_A0;
        ampC64[mult6] = hit.ampC4_A0;
        x6[mult6] = hit.x_0;
        y6[mult6] = hit.y_0;
        mult6++;
        break;

    case 7:
        tofA7[mult7] = hit.tof_A0;
        ampA7[mult7] = hit.amp_A0;
        sumX7[mult7] = hit.sumX_0;
        sumY7[mult7] = hit.sumY_0;
        diffX7[mult7] = hit.diffX_0;
        diffY7[mult7] = hit.diffY_0;
        ratioX7[mult7] = hit.ratioX_0;
        ratioY7[mult7] = hit.ratioY_0;
        sumX7_ampA7[mult7] = hit.sumX_amp_0;
        sumY7_ampA7[mult7] = hit.sumY_amp_0;
        ampC71[mult7] = hit.ampC1_A0;
        ampC72[mult7] = hit.ampC2_A0;
        ampC73[mult7] = hit.ampC3_A0;
        ampC74[mult7] = hit.ampC4_A0;
        x7[mult7] = hit.x_0;
        y7[mult7] = hit.y_0;
        mult7++;
        break;

    case 8:
        tofA8[mult8] = hit.tof_A0;
        ampA8[mult8] = hit.amp_A0;
        sumX8[mult8] = hit.sumX_0;
        sumY8[mult8] = hit.sumY_0;
        diffX8[mult8] = hit.diffX_0;
        diffY8[mult8] = hit.diffY_0;
        ratioX8[mult8] = hit.ratioX_0;
        ratioY8[mult8] = hit.ratioY_0;
        sumX8_ampA8[mult8] = hit.sumX_amp_0;
        sumY8_ampA8[mult8] = hit.sumY_amp_0;
        ampC81[mult8] = hit.ampC1_A0;
        ampC82[mult8] = hit.ampC2_A0;
        ampC83[mult8] = hit.ampC3_A0;
        ampC84[mult8] = hit.ampC4_A0;
        x8[mult8] = hit.x_0;
        y8[mult8] = hit.y_0;
        mult8++;
        break;

    case 9:
        tofA9[mult9] = hit.tof_A0;
        ampA9[mult9] = hit.amp_A0;
        sumX9[mult9] = hit.sumX_0;
        sumY9[mult9] = hit.sumY_0;
        diffX9[mult9] = hit.diffX_0;
        diffY9[mult9] = hit.diffY_0;
        ratioX9[mult9] = hit.ratioX_0;
        ratioY9[mult9] = hit.ratioY_0;
        sumX9_ampA9[mult9] = hit.sumX_amp_0;
        sumY9_ampA9[mult9] = hit.sumY_amp_0;
        ampC91[mult9] = hit.ampC1_A0;
        ampC92[mult9] = hit.ampC2_A0;
        ampC93[mult9] = hit.ampC3_A0;
        ampC94[mult9] = hit.ampC4_A0;
        x9[mult9] = hit.x_0;
        y9[mult9] = hit.y_0;
        mult9++;
        break;
}
    }
    if (!hits.empty()) { // only store valid events with at least one valid hit
        tbdt->Fill();
        events_good++;
    }
}
fout->cd();
tbdt->Write();
fout->Close(); // close the file

std::cout << "File written and closed successfully.\n";
}

        