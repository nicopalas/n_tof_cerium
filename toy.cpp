#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdlib> // for exit()

struct Signal {
    Int_t RunNumber, time, BunchNumber, PSpulse, detn;
    Double_t psTime, tof;
    Float_t PulseIntensity, amp;
    bool used = false;
};

void toy(int run_number, float threshold, double time_for_coincidence) {
    
    std::cout << "INFO: Processing run_number " << run_number << std::endl;

    // --- Open ROOT file ---
    std::string filename = "/nucl_lustre/n_tof_INTC_P_665/DATA/run" + std::to_string(run_number) + ".root";
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Cannot open file " << filename << std::endl;
        return;
    }
    std::cout << "INFO: File opened successfully." << std::endl;

    // --- Load PKUP tree ---
    TTree *pkupTree = nullptr;
    file->GetObject("PKUP", pkupTree);
    if (!pkupTree) {
        std::cerr << "ERROR: TTree 'PKUP' not found in " << filename << std::endl;
        file->Close();
        return;
    }

    Int_t pkup_RunNumber, pkup_time, pkup_BunchNumber, pkup_PSpulse;
    Double_t pkup_psTime, pkup_tof;
    Float_t pkup_PulseIntensity;

    pkupTree->SetBranchStatus("*", 0);
    pkupTree->SetBranchStatus("RunNumber", 1); pkupTree->SetBranchAddress("RunNumber", &pkup_RunNumber);
    pkupTree->SetBranchStatus("time", 1); pkupTree->SetBranchAddress("time", &pkup_time);
    pkupTree->SetBranchStatus("psTime", 1); pkupTree->SetBranchAddress("psTime", &pkup_psTime);
    pkupTree->SetBranchStatus("BunchNumber", 1); pkupTree->SetBranchAddress("BunchNumber", &pkup_BunchNumber);
    pkupTree->SetBranchStatus("PSpulse", 1); pkupTree->SetBranchAddress("PSpulse", &pkup_PSpulse);
    pkupTree->SetBranchStatus("PulseIntensity", 1); pkupTree->SetBranchAddress("PulseIntensity", &pkup_PulseIntensity);
    pkupTree->SetBranchStatus("tof", 1); pkupTree->SetBranchAddress("tof", &pkup_tof);

    Long64_t n_pkup = pkupTree->GetEntriesFast();
    std::cout << "INFO: PKUP entries: " << n_pkup << std::endl;

    // build a map for quick PKUP lookup. key: tuple of beam parameters, value: tof
    // for each PPAN entry we will search the matching PKUP entry using this map (i.e. same beam parameters between both). 
    // that way we can apply the TOF correction easily (tof=tof-pkup_tof).
    std::map<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t>, Double_t> pkup_map;
    for (Long64_t i = 0; i < n_pkup; ++i) {
        pkupTree->GetEntry(i);
        auto key = std::make_tuple(pkup_RunNumber, pkup_time, pkup_psTime,
                                   pkup_BunchNumber, pkup_PSpulse, pkup_PulseIntensity);
        pkup_map[key] = pkup_tof;
    }

    std::cout << "INFO: PKUP lookup map built with " << pkup_map.size() << " entries." << std::endl;

    // PPAN tree
    TTree *tree = nullptr;
    file->GetObject("PPAN", tree);
    if (!tree) {
        std::cerr << "ERROR: TTree 'PPAN' not found in " << filename << std::endl;
        file->Close();
        return;
    }

    Signal s;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("RunNumber", 1); tree->SetBranchAddress("RunNumber", &s.RunNumber);
    tree->SetBranchStatus("time", 1); tree->SetBranchAddress("time", &s.time);
    tree->SetBranchStatus("psTime", 1); tree->SetBranchAddress("psTime", &s.psTime);
    tree->SetBranchStatus("BunchNumber", 1); tree->SetBranchAddress("BunchNumber", &s.BunchNumber);
    tree->SetBranchStatus("PSpulse", 1); tree->SetBranchAddress("PSpulse", &s.PSpulse);
    tree->SetBranchStatus("PulseIntensity", 1); tree->SetBranchAddress("PulseIntensity", &s.PulseIntensity);
    tree->SetBranchStatus("detn", 1); tree->SetBranchAddress("detn", &s.detn);
    tree->SetBranchStatus("tof", 1); tree->SetBranchAddress("tof", &s.tof);
    tree->SetBranchStatus("amp", 1); tree->SetBranchAddress("amp", &s.amp);

    Long64_t nentries = tree->GetEntriesFast();
    if (nentries <= 0) {
        std::cerr << "ERROR: No entries in TTree 'PPAN'!" << std::endl;
        file->Close();
        return;
    }
    std::cout << "INFO: PPAN entries: " << nentries << std::endl;

    std::vector<Signal> signals;
    signals.reserve(nentries);

    // correct tof using PKUP and loading signals
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        //key as a tuple of the beam parameter of a certain signal.
        auto key = std::make_tuple(s.RunNumber, s.time, s.psTime, s.BunchNumber, s.PSpulse, s.PulseIntensity);
        // it represents an iterator to the map element with the matching key.
        // gives us the PKUP signal with same beam parameters as the current PPAN signal. its fir
        auto it = pkup_map.find(key);

        if (it == pkup_map.end()) { // error handling: no matching PKUP entry found
            std::cerr << "FATAL ERROR: No matching PKUP entry found for PPAN signal!"
                      << "\n       RunNumber=" << s.RunNumber
                      << " time=" << s.time
                      << " psTime=" << s.psTime
                      << " BunchNumber=" << s.BunchNumber
                      << " PSpulse=" << s.PSpulse
                      << " PulseIntensity=" << s.PulseIntensity
                      << std::endl;
            file->Close();
            std::exit(EXIT_FAILURE);
        }

        s.tof -= it->second; // tof correction by subtracting tof of the PKUP

        if (s.detn < 7) s.tof += 15; // local correction (due to different time references between DAQ systems). must be corrected strictly with the gamma flash
        if (std::abs(s.amp) > threshold && s.tof > -900) signals.push_back(s); // only keep signals above threshold to avoid unnecessary storage
    }
    

    std::cout << "INFO: Signals loaded and TOF corrected: " << signals.size() << std::endl;

    //sorting signals for coincidence search. this will help grouping them and speed up the process.
    // all the possible candidates for a certain coincidence will be near each other 
    // tof ordering within the same RunNumber, time, psTime, BunchNumber, PSpulse, PulseIntensity.
    std::sort(signals.begin(), signals.end(), [](const Signal &a, const Signal &b) {
        if (a.RunNumber != b.RunNumber) return a.RunNumber < b.RunNumber;
        if (a.time != b.time) return a.time < b.time;
        if (a.psTime != b.psTime) return a.psTime < b.psTime;
        if (a.BunchNumber != b.BunchNumber) return a.BunchNumber < b.BunchNumber;
        if (a.PSpulse != b.PSpulse) return a.PSpulse < b.PSpulse;
        if (a.PulseIntensity != b.PulseIntensity) return a.PulseIntensity < b.PulseIntensity;
        return a.tof < b.tof;
    });
    std::cout << "INFO: Signals sorted." << std::endl;

    // coincidence search
    std::vector<std::vector<Signal>> configurations;

    for (size_t i = 0; i < signals.size(); ++i) {
        if (signals[i].used || std::abs(signals[i].amp) <= threshold) continue;
        std::vector<Signal> coincidence;

        for (size_t j = i + 1; j < signals.size(); ++j) {
            if (std::abs(signals[j].amp) <= threshold) continue;
            // matching criteria in beam-related parameters
            if (signals[j].RunNumber != signals[i].RunNumber ||
                signals[j].time != signals[i].time ||
                signals[j].psTime != signals[i].psTime ||
                signals[j].BunchNumber != signals[i].BunchNumber ||
                signals[j].PSpulse != signals[i].PSpulse ||
                signals[j].PulseIntensity != signals[i].PulseIntensity)
                break;

            if (std::abs(signals[j].tof - signals[i].tof) >= time_for_coincidence) break;
            if (coincidence.size() < 1) {
                signals[i].used = true;
                coincidence.push_back(signals[i]);
            }
            signals[j].used = true;
            coincidence.push_back(signals[j]);
        }

        if (coincidence.size() > 1) configurations.push_back(coincidence);
    }

    std::cout << "INFO: Initial coincidences found: " << configurations.size() << std::endl;
    std::cout << "INFO: Signals used in coincidences: "
              << std::count_if(signals.begin(), signals.end(), [](const Signal &sig){ return sig.used; })
              << std::endl;

        // Lonely coincidences merging
    int merged_signals = 0;

    for (auto &sig : signals) {
        if (sig.used || std::abs(sig.amp) <= threshold) continue;

        bool merged = false;

        for (auto &conf : configurations) {
            if (conf.empty()) continue;

            const auto &first = conf.front();
            const auto &last  = conf.back();

            // Check same beam parameters
            if (sig.RunNumber == first.RunNumber &&
                sig.time == first.time &&
                sig.psTime == first.psTime &&
                sig.BunchNumber == first.BunchNumber &&
                sig.PSpulse == first.PSpulse &&
                sig.PulseIntensity == first.PulseIntensity) {

                // since conf is sorted by tof, check only if there is overlapping with the last or first signal.
                if (std::abs(sig.tof - first.tof) < time_for_coincidence ||
                    std::abs(sig.tof - last.tof) < time_for_coincidence) {

                    sig.used = true;
                    conf.push_back(sig);
                    std::sort(conf.begin(), conf.end(), [](const Signal &a, const Signal &b) {
                        return a.tof < b.tof;
                    });

                    std::cout << "Merging signal detn=" << sig.detn
                              << " tof=" << sig.tof
                              << " into coincidence with range ["
                              << first.tof << ", " << last.tof << "]" << std::endl;

                    merged = true;
                    merged_signals++;
                    break;
                }
            }
        }

        if (!merged) {
            sig.used = true;
        }
    }

    std::cout << "INFO: Lonely signals merged: " << merged_signals << std::endl;

    // Merging overlapping coincidences (only adjacent ones)
    int merge_operations = 0;

    for (size_t i = 0; i + 1 < configurations.size();) {
        auto &conf_i = configurations[i];
        auto &conf_j = configurations[i + 1];

        if (conf_i.empty() || conf_j.empty()) {
            ++i;
            continue;
        }

        const auto &a_first = conf_i.front();
        const auto &a_last  = conf_i.back();
        const auto &b_first = conf_j.front();
        const auto &b_last  = conf_j.back();

        // check same beam parameters
        if (a_first.RunNumber == b_first.RunNumber &&
            a_first.time == b_first.time &&
            a_first.psTime == b_first.psTime &&
            a_first.BunchNumber == b_first.BunchNumber &&
            a_first.PSpulse == b_first.PSpulse &&
            a_first.PulseIntensity == b_first.PulseIntensity) {

            // check for overlap in tof window between the two configurations (again, only need to check the edges)
            // since configurations is sorted by tof, we only need to check the last of the first with the first of the second and vice-versa
            if (std::abs(b_first.tof - a_last.tof) < time_for_coincidence ||
                std::abs(a_first.tof - b_last.tof) < time_for_coincidence) {

                conf_i.insert(conf_i.end(), conf_j.begin(), conf_j.end());
                std::sort(conf_i.begin(), conf_i.end(), [](const Signal &a, const Signal &b) {
                    return a.tof < b.tof;
                });
                configurations.erase(configurations.begin() + i + 1);
                merge_operations++;
                continue; // re-check current i with next configuration
            }
        }

        ++i;
    }

    std::cout << "INFO: Overlapping configurations merged: " << merge_operations << std::endl;
    std::cout << "Configurations:" << std::endl;
    for (size_t i = 0; i < configurations.size(); ++i) {
        std::cout << "Configuration " << i << " (size=" << configurations[i].size() << "):" << std::endl;
        for (const auto& sig : configurations[i]) {
            std::cout << "  detn=" << sig.detn
                      << ", tof=" << sig.tof
                      << ", amp=" << sig.amp
                      << ", used=" << sig.used
                      << std::endl;
        }
    }

    std::cout << "INFO: Merge operations: " << merge_operations << std::endl;
    std::cout << "INFO: Final number of configurations: " << configurations.size() << std::endl;
}
