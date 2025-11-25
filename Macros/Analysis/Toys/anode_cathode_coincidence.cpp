#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <tuple>

const int NDET = 10;     // total de detectores en los datos
const int MAXA = 50;     // máximo número de anodos por detector (entrada)
const int MAXC = 100;    // máximo número de catodos por tipo (entrada)
const int MAXOUT = 400;  // máximo número de catodos almacenados en el árbol de salida
const double L_mm = 100.0;
const double SUM_MIN = 80.0;
const double SUM_MAX = 120.0;


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
    const int NDET_USED = 2;
    int detList[NDET_USED] = {7,8};

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
                        //we will only fill diff if within physical limits
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

    // Calibration: find scale factors vX, vY and means/sigmas of sums
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
            // background initial guess
            double meanContent = (hSumX[idx]->GetNbinsX()>0) ? (hSumX[idx]->Integral() / hSumX[idx]->GetNbinsX()) : 0.0;
            fX->SetParameter(3, meanContent);
            fX->SetParameter(4, 0.0);

            // límites razonables
            fX->SetParLimits(1, 90.0, 110.0);  // center around 100
            fX->SetParLimits(2, 0.5, 5.0);     // expected physical width
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
    TTree *tbdt = new TTree("training","Training tree for BDT");

    // variables
    Int_t event;
    Double_t tofA8, tofA9;
    Double_t dtA, x8, y8, x9, y9;
    Float_t ampA8, ampA9;
    Double_t sumX8, sumY8, diffX8, diffY8, sumX9, sumY9, diffX9, diffY9;
    Float_t ratioX8, ratioY8, ratioX9, ratioY9;
    Float_t sumX8_ampA8, sumY8_ampA8, sumX9_ampA9, sumY9_ampA9;
    Int_t isGood;
    Int_t nCombPerAnode; 

    tbdt->Branch("comb", &nCombPerAnode, "comb/I");   
    tbdt->Branch("event",&event,"event/I");
    tbdt->Branch("dtA",&dtA,"dtA/D");
    tbdt->Branch("tofA8",&tofA8,"tofA8/D");
    tbdt->Branch("tofA9",&tofA9,"tofA9/D");
    tbdt->Branch("ampA8",&ampA8,"ampA8/F");
    tbdt->Branch("ampA9",&ampA9,"ampA9/F");
    tbdt->Branch("sumX8",&sumX8,"sumX8/D");
    tbdt->Branch("sumY8",&sumY8,"sumY8/D");
    tbdt->Branch("x8", &x8, "x8/D");
    tbdt->Branch("y8", &y8, "y8/D");
    tbdt->Branch("diffX8",&diffX8,"diffX8/D");
    tbdt->Branch("diffY8",&diffY8,"diffY8/D");
    tbdt->Branch("x9", &x9, "x9/D");
    tbdt->Branch("y9", &y9, "y9/D");
    tbdt->Branch("ratioX8",&ratioX8,"ratioX8/F");
    tbdt->Branch("ratioY8",&ratioY8,"ratioY8/F");
    tbdt->Branch("sumX8_ampA8",&sumX8_ampA8,"sumX8_ampA8/F");
    tbdt->Branch("sumY8_ampA8",&sumY8_ampA8,"sumY8_ampA8/F");

    tbdt->Branch("sumX9",&sumX9,"sumX9/D");
    tbdt->Branch("sumY9",&sumY9,"sumY9/D");
    tbdt->Branch("diffX9",&diffX9,"diffX9/D");
    tbdt->Branch("diffY9",&diffY9,"diffY9/D");
    tbdt->Branch("ratioX9",&ratioX9,"ratioX9/F");
    tbdt->Branch("ratioY9",&ratioY9,"ratioY9/F");
    tbdt->Branch("sumX9_ampA9",&sumX9_ampA9,"sumX9_ampA9/F");
    tbdt->Branch("sumY9_ampA9",&sumY9_ampA9,"sumY9_ampA9/F");

    tbdt->Branch("isGood",&isGood,"isGood/I");
    int good_combinations = 0;
    typedef std::tuple<
    double, double, float, float,        // tofA8, tofA9, ampA8, ampA9
    double, double, double, double,      // sumX8, sumY8, diffX8, diffY8
    float, float, float, float,          // ratioX8, ratioY8, sumX8_ampA8, sumY8_ampA8
    double, double, double, double,      // sumX9, sumY9, diffX9, diffY9
    float, float, float, float, double, double, double,double           // ratioX9, ratioY9, sumX9_ampA9, sumY9_ampA9
> ComboKey;

std::set<ComboKey> allComb; 

for(Long64_t ev = 0; ev < nentries; ++ev) {
    tin->GetEntry(ev);
    event = ev;

    if(ev % 1000000 == 0) std::cout << "Processing event " << ev << "\n";

    int d8 = detList[0];
    int d9 = detList[1];

    if(multA[d8]<1 || multA[d9]<1) continue;

    for(int i8 = 0; i8 < multA[d8]; ++i8) {

            tofA8 = tofA[d8][i8]; ampA8 = ampA[d8][i8];
            tofA9 = tofA[d9][i9]; ampA9 = ampA[d9][i9];
            dtA = tofA9 - tofA8;

            if(fabs(dtA) >= 8.0) continue;

            // ----------------------
            // search cathode combinations detector 8
            // ----------------------
            struct Combo8 { int k1,k2,k3,k4; double sumX,sumY,diffX,diffY,ratioX,ratioY,sumX_amp,sumY_amp,x,y; };
            std::vector<Combo8> valid8;
            if (multC1[d8]<1 || multC2[d8]<1 || multC3[d8]<1 || multC4[d8]<1) continue;

            for(int k1=0; k1<multC1[d8]; ++k1)
            for(int k2=0; k2<multC2[d8]; ++k2)
            for(int k3=0; k3<multC3[d8]; ++k3)
            for(int k4=0; k4<multC4[d8]; ++k4){
                double sX = tofC1[d8][k1] + tofC2[d8][k2] - 2*tofA8;
                double dX = tofC2[d8][k2] - tofC1[d8][k1];
                double x_8 = (dX * vX[0]) / 2.0;
                double rX = (ampC2[d8][k2]!=0) ? ampC1[d8][k1]/ampC2[d8][k2] : 0;
                double sX_amp = (ampC1[d8][k1]+ampC2[d8][k2])/ampA8;

                double sY = tofC3[d8][k3] + tofC4[d8][k4] - 2*tofA8;
                double dY = tofC4[d8][k4] - tofC3[d8][k3];
                double y_8 = (dY * vY[0]) / 2.0;
                double rY = (ampC4[d8][k4]!=0) ? ampC3[d8][k3]/ampC4[d8][k4] : 0;
                double sY_amp = (ampC3[d8][k3]+ampC4[d8][k4])/ampA8;
                if(rY>=0.5 && rY<=2.0 && sY_amp>=0.2 && sY_amp<=1.6 &&
                   fabs(sY-meanSumY[0])<=2*meanWidthY[0] && fabs(dY)<=L_mm/vY[0] &&
                   rX>=0.5 && rX<=2.0 && sX_amp>=0.2 && sX_amp<=1.6 &&
                   fabs(sX-meanSumX[0])<=2*meanWidthX[0] && fabs(dX)<=L_mm/vX[0]) {
                    valid8.push_back({k1,k2,k3,k4,sX,sY,dX,dY,rX,rY,sX_amp,sY_amp, x_8,y_8});
                }
            }

            if(valid8.empty()) continue; //pass
            // ----------------------
            // search cathode combinations detector 9
            // ----------------------
            if (multC1[d9]<1 || multC2[d9]<1 || multC3[d9]<1 || multC4[d9]<1) continue;
            struct Combo9 { int l1,l2,l3,l4; double sumX,sumY,diffX,diffY,ratioX,ratioY,sumX_amp,sumY_amp, x,y; };
            std::vector<Combo9> valid9;

            for(int l1=0; l1<multC1[d9]; ++l1)
            for(int l2=0; l2<multC2[d9]; ++l2)
            for(int l3=0; l3<multC3[d9]; ++l3)
            for(int l4=0; l4<multC4[d9]; ++l4){
                double sX = tofC1[d9][l1] + tofC2[d9][l2] - 2*tofA9;
                double dX = tofC2[d9][l2] - tofC1[d9][l1];
                double x_9 = (dX * vX[1]) / 2.0;
                double rX = (ampC2[d9][l2]!=0) ? ampC1[d9][l1]/ampC2[d9][l2] : 0;
                double sX_amp = (ampC1[d9][l1]+ampC2[d9][l2])/ampA9;

                double sY = tofC3[d9][l3] + tofC4[d9][l4] - 2*tofA9;
                double dY = tofC4[d9][l4] - tofC3[d9][l3];
                double y_9 = (dY * vY[1]) / 2.0;
                double rY = (ampC4[d9][l4]!=0) ? ampC3[d9][l3]/ampC4[d9][l4] : 0;
                double sY_amp = (ampC3[d9][l3]+ampC4[d9][l4])/ampA9;

                if(rY>=0.5 && rY<=2.0 && sY_amp>=0.2 && sY_amp<=1.6 &&
                   fabs(sY-meanSumY[1])<=2*meanWidthY[1] && fabs(dY)<=L_mm/vY[1] &&
                   rX>=0.5 && rX<=2.0 && sX_amp>=0.2 && sX_amp<=1.6 &&
                   fabs(sX-meanSumX[1])<=2*meanWidthX[1] && fabs(dX)<=L_mm/vX[1]) {
                    valid9.push_back({l1,l2,l3,l4,sX,sY,dX,dY,rX,rY,sX_amp,sY_amp, x_9,y_9});
                }
            }

            if(valid9.empty()) continue; // pass to next anode pair

            // ----------------------
            // combine valid combos from det 8 and 9
            // ----------------------
            nCombPerAnode = 0;
            for(auto &c8 : valid8){
                for(auto &c9 : valid9){

                    sumX8 = c8.sumX; sumY8 = c8.sumY;
                    diffX8 = c8.diffX; diffY8 = c8.diffY;
                    ratioX8 = c8.ratioX; ratioY8 = c8.ratioY;
                    sumX8_ampA8 = c8.sumX_amp; sumY8_ampA8 = c8.sumY_amp;
                    x8 = c8.x; y8 = c8.y;

                    sumX9 = c9.sumX; sumY9 = c9.sumY;
                    diffX9 = c9.diffX; diffY9 = c9.diffY;
                    ratioX9 = c9.ratioX; ratioY9 = c9.ratioY;
                    sumX9_ampA9 = c9.sumX_amp; sumY9_ampA9 = c9.sumY_amp;
                    x9 = c9.x; y9 = c9.y;

                    tofA8 = tofA[d8][i8]; ampA8 = ampA[d8][i8];
                    tofA9 = tofA[d9][i9]; ampA9 = ampA[d9][i9];

                    // Crear key para filtrar duplicados globales
                    ComboKey key = std::make_tuple(
                        tofA8,tofA9,ampA8,ampA9,
                        sumX8,sumY8,diffX8,diffY8,
                        ratioX8,ratioY8,sumX8_ampA8,sumY8_ampA8,
                        sumX9,sumY9,diffX9,diffY9,
                        ratioX9,ratioY9,sumX9_ampA9,sumY9_ampA9, x8, y8, x9, y9
                    );

                    if(allComb.count(key) == 0){
                        isGood = 1;
                        tbdt->Fill();
                        good_combinations++;
                        nCombPerAnode++;
                        allComb.insert(key); // registrar combinación global
                    }
                }
            }

        }
    }
}

std::cout << "Total good combinations: " << good_combinations << std::endl;



fout->cd();
tbdt->Write();
fout->Close();
fin->Close();
std::cout << "Training tree written to " << outfile << "\n";
}