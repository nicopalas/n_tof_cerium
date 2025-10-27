// attach_cathodes_to_anodes.cpp
// Final version per user's specification.
// Compile with: g++ -O2 `root-config --cflags --libs` -o attach_cathodes_to_anodes attach_cathodes_to_anodes.cpp
// Run: ./attach_cathodes_to_anodes <run_number> [amp_threshold=0.0] [time_for_coincidence=50.0]

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>

#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iomanip>

struct Signal {
    Int_t RunNumber=0;
    Int_t time=0;
    Double_t psTime=0;
    Int_t BunchNumber=0;
    Int_t PSpulse=0;
    Float_t PulseIntensity=0;
    Int_t detn=0;
    Double_t tof=0;
    Float_t amp=0;
    bool used=false;
};

static bool is_valid_cathode(int detn) {
    if (detn <= 0 || detn > 99) return false;
    int ones = detn % 10;
    int tens = detn / 10;
    if (tens < 0 || tens > 9) return false;
    return (ones >= 1 && ones <= 4);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run_number> [amp_threshold=0.0] [time_for_coincidence=50.0]\n";
        return 1;
    }
    int run_number = std::atoi(argv[1]);
    float amp_threshold = 0.0;
    double time_for_coincidence = 50.0;
    if (argc >= 3) amp_threshold = std::atof(argv[2]);
    if (argc >= 4) time_for_coincidence = std::atof(argv[3]);

    std::cout << "[INFO] attach_cathodes_to_anodes starting for run " << run_number
              << " (amp_threshold=" << amp_threshold
              << ", time_for_coincidence=" << time_for_coincidence << ")\n";

    // Paths (adjust if needed)
    std::string data_file = "/nucl_lustre/n_tof_INTC_P_665/DATA/run" + std::to_string(run_number) + ".root";
    std::string anode_file = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/coincidences_raw/toy/out_run" + std::to_string(run_number) + ".root";
    std::string out_file = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Cathodes/toy/out_cathodes_" + std::to_string(run_number) + ".root";

    // Open data file
    TFile *fdata = TFile::Open(data_file.c_str());
    if (!fdata || fdata->IsZombie()) {
        std::cerr << "[ERROR] cannot open data file: " << data_file << std::endl;
        return 2;
    }

    // PKUP tree -> build map beam params -> pkup tof
    TTree *pkupTree = nullptr;
    fdata->GetObject("PKUP", pkupTree);
    if (!pkupTree) {
        std::cerr << "[ERROR] PKUP tree not found in " << data_file << std::endl;
        fdata->Close();
        return 3;
    }

    Int_t pk_RunNumber, pk_time, pk_BunchNumber, pk_PSpulse;
    Double_t pk_psTime, pk_tof;
    Float_t pk_PulseIntensity;
    pkupTree->SetBranchStatus("*",0);
    pkupTree->SetBranchStatus("RunNumber",1); pkupTree->SetBranchAddress("RunNumber",&pk_RunNumber);
    pkupTree->SetBranchStatus("time",1); pkupTree->SetBranchAddress("time",&pk_time);
    pkupTree->SetBranchStatus("psTime",1); pkupTree->SetBranchAddress("psTime",&pk_psTime);
    pkupTree->SetBranchStatus("BunchNumber",1); pkupTree->SetBranchAddress("BunchNumber",&pk_BunchNumber);
    pkupTree->SetBranchStatus("PSpulse",1); pkupTree->SetBranchAddress("PSpulse",&pk_PSpulse);
    pkupTree->SetBranchStatus("PulseIntensity",1); pkupTree->SetBranchAddress("PulseIntensity",&pk_PulseIntensity);
    pkupTree->SetBranchStatus("tof",1); pkupTree->SetBranchAddress("tof",&pk_tof);

    std::map<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t>, Double_t> pkup_map;
    Long64_t n_pk = pkupTree->GetEntriesFast();
    for (Long64_t i=0;i<n_pk;++i) {
        pkupTree->GetEntry(i);
        auto key = std::make_tuple(pk_RunNumber, pk_time, pk_psTime, pk_BunchNumber, pk_PSpulse, pk_PulseIntensity);
        pkup_map[key] = pk_tof;
    }
    std::cout << "[INFO] PKUP map built with " << pkup_map.size() << " entries\n";

    // PPAC tree -> read signals and correct TOF
    TTree *ppacTree = nullptr;
    fdata->GetObject("PPAC", ppacTree);
    if (!ppacTree) {
        std::cerr << "[ERROR] PPAC tree not found in " << data_file << std::endl;
        fdata->Close();
        return 4;
    }

    Signal tmp;
    ppacTree->SetBranchStatus("*",0);
    ppacTree->SetBranchStatus("RunNumber",1); ppacTree->SetBranchAddress("RunNumber",&tmp.RunNumber);
    ppacTree->SetBranchStatus("time",1); ppacTree->SetBranchAddress("time",&tmp.time);
    ppacTree->SetBranchStatus("psTime",1); ppacTree->SetBranchAddress("psTime",&tmp.psTime);
    ppacTree->SetBranchStatus("BunchNumber",1); ppacTree->SetBranchAddress("BunchNumber",&tmp.BunchNumber);
    ppacTree->SetBranchStatus("PSpulse",1); ppacTree->SetBranchAddress("PSpulse",&tmp.PSpulse);
    ppacTree->SetBranchStatus("PulseIntensity",1); ppacTree->SetBranchAddress("PulseIntensity",&tmp.PulseIntensity);
    ppacTree->SetBranchStatus("detn",1); ppacTree->SetBranchAddress("detn",&tmp.detn);
    ppacTree->SetBranchStatus("tof",1); ppacTree->SetBranchAddress("tof",&tmp.tof);
    ppacTree->SetBranchStatus("amp",1); ppacTree->SetBranchAddress("amp",&tmp.amp);

    std::vector<Signal> ppac_signals;
    Long64_t n_ppac = ppacTree->GetEntriesFast();
    ppac_signals.reserve(n_ppac);

    for (Long64_t i=0;i<n_ppac;++i) {
        ppacTree->GetEntry(i);
        auto key = std::make_tuple(tmp.RunNumber, tmp.time, tmp.psTime, tmp.BunchNumber, tmp.PSpulse, tmp.PulseIntensity);
        auto it = pkup_map.find(key);
        if (it == pkup_map.end()) {
            // skip if no pkup match
            continue;
        }
        tmp.tof -= it->second;          // PKUP correction
        if (tmp.detn < 63) tmp.tof += 15; // local correction (as in toy.cpp)
        if (std::fabs(tmp.amp) > amp_threshold && tmp.tof > -900) {
            tmp.used = false;
            ppac_signals.push_back(tmp);
        }
    }
    std::cout << "[INFO] PPAC signals loaded after PKUP correction: " << ppac_signals.size() << std::endl;

    // Sort PPAC signals and build index: beam params -> list of indices
    std::sort(ppac_signals.begin(), ppac_signals.end(), [](const Signal &a, const Signal &b){
        if (a.RunNumber != b.RunNumber) return a.RunNumber < b.RunNumber;
        if (a.time != b.time) return a.time < b.time;
        if (a.psTime != b.psTime) return a.psTime < b.psTime;
        if (a.BunchNumber != b.BunchNumber) return a.BunchNumber < b.BunchNumber;
        if (a.PSpulse != b.PSpulse) return a.PSpulse < b.PSpulse;
        if (a.PulseIntensity != b.PulseIntensity) return a.PulseIntensity < b.PulseIntensity;
        return a.tof < b.tof;
    });

    std::map<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t>, std::vector<size_t>> ppac_index;
    for (size_t i=0;i<ppac_signals.size(); ++i) {
        auto &s = ppac_signals[i];
        auto key = std::make_tuple(s.RunNumber, s.time, s.psTime, s.BunchNumber, s.PSpulse, s.PulseIntensity);
        ppac_index[key].push_back(i);
    }

    // Open anode coincidences file and tree
    TFile *fan = TFile::Open(anode_file.c_str());
    if (!fan || fan->IsZombie()) {
        std::cerr << "[ERROR] cannot open anode file: " << anode_file << std::endl;
        fdata->Close();
        return 5;
    }
    TTree *anodeTree = nullptr;
    fan->GetObject("nTOF_coincidences", anodeTree);
    if (!anodeTree) {
        std::cerr << "[ERROR] nTOF_coincidences not found in " << anode_file << std::endl;
        fan->Close(); fdata->Close();
        return 6;
    }

    // Read anode branches
    Int_t A_RunNumber=0, A_time=0, A_BunchNumber=0, A_PSpulse=0;
    Double_t A_psTime=0;
    Float_t A_PulseIntensity=0;
    const int MAXAN = 40;
    Int_t A_mult0=0,A_mult1=0,A_mult2=0,A_mult3=0,A_mult4=0,A_mult5=0,A_mult6=0,A_mult7=0,A_mult8=0,A_mult9=0;
    Double_t A_tof0[MAXAN], A_tof1[MAXAN], A_tof2[MAXAN], A_tof3[MAXAN], A_tof4[MAXAN], A_tof5[MAXAN], A_tof6[MAXAN], A_tof7[MAXAN], A_tof8[MAXAN], A_tof9[MAXAN];
    Float_t  A_amp0[MAXAN], A_amp1[MAXAN], A_amp2[MAXAN], A_amp3[MAXAN], A_amp4[MAXAN], A_amp5[MAXAN], A_amp6[MAXAN], A_amp7[MAXAN], A_amp8[MAXAN], A_amp9[MAXAN];

    anodeTree->SetBranchStatus("*",1);
    anodeTree->SetBranchAddress("RunNumber",&A_RunNumber);
    anodeTree->SetBranchAddress("time",&A_time);
    anodeTree->SetBranchAddress("psTime",&A_psTime);
    anodeTree->SetBranchAddress("BunchNumber",&A_BunchNumber);
    anodeTree->SetBranchAddress("PSpulse",&A_PSpulse);
    anodeTree->SetBranchAddress("PulseIntensity",&A_PulseIntensity);

    anodeTree->SetBranchAddress("mult0",&A_mult0);
    anodeTree->SetBranchAddress("mult1",&A_mult1);
    anodeTree->SetBranchAddress("mult2",&A_mult2);
    anodeTree->SetBranchAddress("mult3",&A_mult3);
    anodeTree->SetBranchAddress("mult4",&A_mult4);
    anodeTree->SetBranchAddress("mult5",&A_mult5);
    anodeTree->SetBranchAddress("mult6",&A_mult6);
    anodeTree->SetBranchAddress("mult7",&A_mult7);
    anodeTree->SetBranchAddress("mult8",&A_mult8);
    anodeTree->SetBranchAddress("mult9",&A_mult9);

    anodeTree->SetBranchAddress("tof0",A_tof0);
    anodeTree->SetBranchAddress("tof1",A_tof1);
    anodeTree->SetBranchAddress("tof2",A_tof2);
    anodeTree->SetBranchAddress("tof3",A_tof3);
    anodeTree->SetBranchAddress("tof4",A_tof4);
    anodeTree->SetBranchAddress("tof5",A_tof5);
    anodeTree->SetBranchAddress("tof6",A_tof6);
    anodeTree->SetBranchAddress("tof7",A_tof7);
    anodeTree->SetBranchAddress("tof8",A_tof8);
    anodeTree->SetBranchAddress("tof9",A_tof9);

    anodeTree->SetBranchAddress("amp0",A_amp0);
    anodeTree->SetBranchAddress("amp1",A_amp1);
    anodeTree->SetBranchAddress("amp2",A_amp2);
    anodeTree->SetBranchAddress("amp3",A_amp3);
    anodeTree->SetBranchAddress("amp4",A_amp4);
    anodeTree->SetBranchAddress("amp5",A_amp5);
    anodeTree->SetBranchAddress("amp6",A_amp6);
    anodeTree->SetBranchAddress("amp7",A_amp7);
    anodeTree->SetBranchAddress("amp8",A_amp8);
    anodeTree->SetBranchAddress("amp9",A_amp9);

    Long64_t n_anodes = anodeTree->GetEntriesFast();
    std::cout << "[INFO] anode coincidences count: " << n_anodes << std::endl;

    // Prepare output file and tree
    TFile *fout = TFile::Open(out_file.c_str(),"RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] cannot create output file: " << out_file << std::endl;
        fan->Close(); fdata->Close();
        return 7;
    }

    TTree *out = new TTree("coincidences","anodes + cathodes coincidences");

    // Output anode branches (we will write only kept anode tof entries)
    Int_t out_RunNumber, out_time, out_BunchNumber, out_PSpulse;
    Double_t out_psTime;
    Float_t out_PulseIntensity;
    Int_t out_mult0,out_mult1,out_mult2,out_mult3,out_mult4,out_mult5,out_mult6,out_mult7,out_mult8,out_mult9;
    const int MAXOUT = 80; // max entries per anode bucket (adjustable)
    Double_t out_tof0[MAXOUT], out_tof1[MAXOUT], out_tof2[MAXOUT], out_tof3[MAXOUT], out_tof4[MAXOUT],
             out_tof5[MAXOUT], out_tof6[MAXOUT], out_tof7[MAXOUT], out_tof8[MAXOUT], out_tof9[MAXOUT];
    Float_t  out_amp0[MAXOUT], out_amp1[MAXOUT], out_amp2[MAXOUT], out_amp3[MAXOUT], out_amp4[MAXOUT],
             out_amp5[MAXOUT], out_amp6[MAXOUT], out_amp7[MAXOUT], out_amp8[MAXOUT], out_amp9[MAXOUT];

    out->Branch("RunNumber",&out_RunNumber,"RunNumber/I");
    out->Branch("time",&out_time,"time/I");
    out->Branch("psTime",&out_psTime,"psTime/D");
    out->Branch("BunchNumber",&out_BunchNumber,"BunchNumber/I");
    out->Branch("PSpulse",&out_PSpulse,"PSpulse/I");
    out->Branch("PulseIntensity",&out_PulseIntensity,"PulseIntensity/F");

    out->Branch("mult0",&out_mult0,"mult0/I");
    out->Branch("mult1",&out_mult1,"mult1/I");
    out->Branch("mult2",&out_mult2,"mult2/I");
    out->Branch("mult3",&out_mult3,"mult3/I");
    out->Branch("mult4",&out_mult4,"mult4/I");
    out->Branch("mult5",&out_mult5,"mult5/I");
    out->Branch("mult6",&out_mult6,"mult6/I");
    out->Branch("mult7",&out_mult7,"mult7/I");
    out->Branch("mult8",&out_mult8,"mult8/I");
    out->Branch("mult9",&out_mult9,"mult9/I");

    out->Branch("tof0",out_tof0,"tof0[out_mult0]/D");
    out->Branch("tof1",out_tof1,"tof1[out_mult1]/D");
    out->Branch("tof2",out_tof2,"tof2[out_mult2]/D");
    out->Branch("tof3",out_tof3,"tof3[out_mult3]/D");
    out->Branch("tof4",out_tof4,"tof4[out_mult4]/D");
    out->Branch("tof5",out_tof5,"tof5[out_mult5]/D");
    out->Branch("tof6",out_tof6,"tof6[out_mult6]/D");
    out->Branch("tof7",out_tof7,"tof7[out_mult7]/D");
    out->Branch("tof8",out_tof8,"tof8[out_mult8]/D");
    out->Branch("tof9",out_tof9,"tof9[out_mult9]/D");

    out->Branch("amp0",out_amp0,"amp0[out_mult0]/F");
    out->Branch("amp1",out_amp1,"amp1[out_mult1]/F");
    out->Branch("amp2",out_amp2,"amp2[out_mult2]/F");
    out->Branch("amp3",out_amp3,"amp3[out_mult3]/F");
    out->Branch("amp4",out_amp4,"amp4[out_mult4]/F");
    out->Branch("amp5",out_amp5,"amp5[out_mult5]/F");
    out->Branch("amp6",out_amp6,"amp6[out_mult6]/F");
    out->Branch("amp7",out_amp7,"amp7[out_mult7]/F");
    out->Branch("amp8",out_amp8,"amp8[out_mult8]/F");
    out->Branch("amp9",out_amp9,"amp9[out_mult9]/F");

    // Cathode branches: create tofXX, ampXX, multXX for detn in 01..94 (valid cathodes only)
    const int MAX_DETN = 95; // index 0 unused
    const int MAXCAT = 100;  // max entries per cathode (adjustable)

    // allocate arrays
    std::vector<int> mult_cat(MAX_DETN,0);
    std::vector<std::vector<double>> tof_cat(MAX_DETN, std::vector<double>(MAXCAT, 0.0));
    std::vector<std::vector<float>> amp_cat(MAX_DETN, std::vector<float>(MAXCAT, 0.0f));

    // keep branch descriptor strings alive
    std::vector<std::string> branch_descriptors;

    for (int detn = 1; detn <= 94; ++detn) {
        if (!is_valid_cathode(detn)) continue;
        // branch names: e.g. "tof01", "amp01", "mult01"
        std::ostringstream sdet; sdet << std::setw(2) << std::setfill('0') << detn;
        std::string detstr = sdet.str();

        std::string bmult = "mult" + detstr;
        std::string btof  = "tof"  + detstr;
        std::string bamp  = "amp"  + detstr;

        // descriptor strings must live until file write - store in vector
        std::string desc_mult = bmult + "/I";
        std::string desc_tof  = btof  + "[ " + bmult + " ]/D";
        std::string desc_amp  = bamp  + "[ " + bmult + " ]/F";

        branch_descriptors.push_back(desc_mult);
        branch_descriptors.push_back(desc_tof);
        branch_descriptors.push_back(desc_amp);

        out->Branch(bmult.c_str(), &mult_cat[detn], desc_mult.c_str());
        out->Branch(btof.c_str(), tof_cat[detn].data(), desc_tof.c_str());
        out->Branch(bamp.c_str(), amp_cat[detn].data(), desc_amp.c_str());
    }

    // Main event loop: for each anode event, attach cathodes only to their corresponding anode
    for (Long64_t ev = 0; ev < n_anodes; ++ev) {
        anodeTree->GetEntry(ev);

        // copy beam info
        out_RunNumber = A_RunNumber;
        out_time = A_time;
        out_psTime = A_psTime;
        out_BunchNumber = A_BunchNumber;
        out_PSpulse = A_PSpulse;
        out_PulseIntensity = A_PulseIntensity;

        // reset outputs
        out_mult0=out_mult1=out_mult2=out_mult3=out_mult4=out_mult5=out_mult6=out_mult7=out_mult8=out_mult9=0;
        memset(out_tof0,0,sizeof(out_tof0)); memset(out_tof1,0,sizeof(out_tof1)); memset(out_tof2,0,sizeof(out_tof2));
        memset(out_tof3,0,sizeof(out_tof3)); memset(out_tof4,0,sizeof(out_tof4)); memset(out_tof5,0,sizeof(out_tof5));
        memset(out_tof6,0,sizeof(out_tof6)); memset(out_tof7,0,sizeof(out_tof7)); memset(out_tof8,0,sizeof(out_tof8));
        memset(out_tof9,0,sizeof(out_tof9));
        memset(out_amp0,0,sizeof(out_amp0)); memset(out_amp1,0,sizeof(out_amp1)); memset(out_amp2,0,sizeof(out_amp2));
        memset(out_amp3,0,sizeof(out_amp3)); memset(out_amp4,0,sizeof(out_amp4)); memset(out_amp5,0,sizeof(out_amp5));
        memset(out_amp6,0,sizeof(out_amp6)); memset(out_amp7,0,sizeof(out_amp7)); memset(out_amp8,0,sizeof(out_amp8));
        memset(out_amp9,0,sizeof(out_amp9));

        // reset cathode arrays for valid detn
        for (int d=1; d<=94; ++d) {
            if (!is_valid_cathode(d)) continue;
            mult_cat[d] = 0;
            std::fill(tof_cat[d].begin(), tof_cat[d].end(), 0.0);
            std::fill(amp_cat[d].begin(), amp_cat[d].end(), 0.0f);
        }

        // find ppac signals with same beam params
        auto key = std::make_tuple(out_RunNumber, out_time, out_psTime, out_BunchNumber, out_PSpulse, out_PulseIntensity);
        auto it = ppac_index.find(key);

        // We'll build per-anode vectors of anode tof/amp that DO have >=1 cathode match
        std::vector<double> keep_tof[10];
        std::vector<float>  keep_amp[10];

        if (it != ppac_index.end()) {
            // For each anode index and each anode tof entry, search only its corresponding cathodes
            for (int an = 0; an <= 9; ++an) {
                int base = an * 10;
                Double_t *A_tof_arr = nullptr;
                Float_t  *A_amp_arr = nullptr;
                Int_t A_mult = 0;
                switch(an) {
                    case 0: A_tof_arr = A_tof0; A_amp_arr = A_amp0; A_mult = A_mult0; break;
                    case 1: A_tof_arr = A_tof1; A_amp_arr = A_amp1; A_mult = A_mult1; break;
                    case 2: A_tof_arr = A_tof2; A_amp_arr = A_amp2; A_mult = A_mult2; break;
                    case 3: A_tof_arr = A_tof3; A_amp_arr = A_amp3; A_mult = A_mult3; break;
                    case 4: A_tof_arr = A_tof4; A_amp_arr = A_amp4; A_mult = A_mult4; break;
                    case 5: A_tof_arr = A_tof5; A_amp_arr = A_amp5; A_mult = A_mult5; break;
                    case 6: A_tof_arr = A_tof6; A_amp_arr = A_amp6; A_mult = A_mult6; break;
                    case 7: A_tof_arr = A_tof7; A_amp_arr = A_amp7; A_mult = A_mult7; break;
                    case 8: A_tof_arr = A_tof8; A_amp_arr = A_amp8; A_mult = A_mult8; break;
                    case 9: A_tof_arr = A_tof9; A_amp_arr = A_amp9; A_mult = A_mult9; break;
                }

                for (int k = 0; k < A_mult; ++k) {
                    double tof_an = A_tof_arr[k];
                    float amp_an = A_amp_arr[k];

                    bool anode_has_match = false;

                    // scan all ppac signals in this beam (only these indices)
                    for (size_t idx : it->second) {
                        Signal &s = ppac_signals[idx];
                        // We MUST mark as used AFTER checking; so first check, then set used regardless
                        if (!is_valid_cathode(s.detn)) {
                            // still mark as used after checking (per user's request)
                            s.used = true;
                            continue;
                        }
                        // check detn belongs to this anode's cathode set
                        if (s.detn < base+1 || s.detn > base+4) {
                            // not in this anode's cathode group -> check completes; mark used
                            s.used = true;
                            continue;
                        }

                        // Now s.detn is one of base+1..base+4. Check tof window.
                        if (s.tof >= tof_an && s.tof <= tof_an + 250.0) {
                            // match: append to cathode arrays for this detn
                            int d = s.detn;
                            if (mult_cat[d] < MAXCAT) {
                                tof_cat[d][mult_cat[d]] = s.tof;
                                amp_cat[d][mult_cat[d]] = s.amp;
                                mult_cat[d] += 1;
                            }
                            anode_has_match = true;
                        }
                        // After checking this signal for this anode entry, mark it as used (whether matched or not)
                        s.used = true;
                    } // end loop over ppac entries for this beam

                    // If this anode tof had >=1 matching cathode, keep the anode tof/amp in output for this anode
                    if (anode_has_match) {
                        keep_tof[an].push_back(tof_an);
                        keep_amp[an].push_back(amp_an);
                    }
                } // end loop anode tof entries
            } // end loop anodes 0..9
        } else {
            // No PPAC signals with same beam params: nothing to match.
            // Per rule, we will not keep any anode (all mults remain zero) and event will be skipped.
        }

        // Now copy kept anode buckets into out arrays (only anode entries that had >=1 cat match)
        auto copy_kept = [&](int anode_idx, Double_t *out_tofs, Float_t *out_amps, Int_t &out_mult) {
            out_mult = 0;
            for (size_t ii = 0; ii < keep_tof[anode_idx].size() && out_mult < MAXOUT; ++ii) {
                out_tofs[out_mult] = keep_tof[anode_idx][ii];
                out_amps[out_mult] = keep_amp[anode_idx][ii];
                out_mult++;
            }
        };

        copy_kept(0, out_tof0, out_amp0, out_mult0);
        copy_kept(1, out_tof1, out_amp1, out_mult1);
        copy_kept(2, out_tof2, out_amp2, out_mult2);
        copy_kept(3, out_tof3, out_amp3, out_mult3);
        copy_kept(4, out_tof4, out_amp4, out_mult4);
        copy_kept(5, out_tof5, out_amp5, out_mult5);
        copy_kept(6, out_tof6, out_amp6, out_mult6);
        copy_kept(7, out_tof7, out_amp7, out_mult7);
        copy_kept(8, out_tof8, out_amp8, out_mult8);
        copy_kept(9, out_tof9, out_amp9, out_mult9);

        // Only write the event if at least one anode has >=1 cathode match (i.e., out_multX > 0 for some X)
        bool any_anode_kept = (out_mult0>0)||(out_mult1>0)||(out_mult2>0)||(out_mult3>0)||(out_mult4>0) ||
                              (out_mult5>0)||(out_mult6>0)||(out_mult7>0)||(out_mult8>0)||(out_mult9>0);

        if (any_anode_kept) {
            out->Fill();
        } else {
            // not saved
        }
    } // end event loop

    // write and close
    out->Write();
    fout->Close();
    fan->Close();
    fdata->Close();

    std::cout << "[INFO] Finished. Output written to " << out_file << std::endl;
    return 0;
}
