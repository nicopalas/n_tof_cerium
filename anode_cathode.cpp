#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <set>
#include <TF1.h>
#include <TMath.h>
#include <vector>
Int_t RunNumber, eventTime, BunchNumber, PSpulse;
Double_t psTime;
Double_t tof;
Float_t PulseIntensity;
Int_t mult=1;
Int_t time_for_coincidence = 300; //time difference between signals CHECK IF TIME IS IN ns!!!!
Int_t max_mul = 5; //max number of signals in coincidence in a single detector
Int_t mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
Double_t tof0[5],tof1[5],tof2[5],tof3[5],tof4[5],tof5[5],tof6[5],tof7[5],tof8[5],tof9[5]; //must be hardcoded if arrays are defined at file scope
Float_t amp0[5],amp1[5],amp2[5],amp3[5],amp4[5],amp5[5],amp6[5],amp7[5],amp8[5],amp9[5];  //must be hardcoded if arrays are defined at file scope
Double_t gamma_flash0[5], gamma_flash1[5], gamma_flash2[5], gamma_flash3[5], gamma_flash4[5], gamma_flash5[5], gamma_flash6[5], gamma_flash7[5], gamma_flash8[5], gamma_flash9[5];
Int_t detn_all[10]; 
Int_t detn_cathode[40];
Double_t tof01[20], tof02[20], tof03[20], tof04[20], tof11[20], tof12[20], tof13[20], tof14[20], tof21[20], tof22[20], tof23[20], tof24[20], tof31[20], tof32[20], tof33[20], tof34[20], tof41[20], tof42[20], tof43[20], tof44[20], tof51[20], tof52[20], tof53[20], tof54[20], tof61[20], tof62[20], tof63[20], tof64[20], tof71[20], tof72[20], tof73[20], tof74[20], tof81[20], tof82[20], tof83[20], tof84[20], tof91[20], tof92[20], tof93[20], tof94[20];
Float_t amp01[20], amp02[20], amp03[20], amp04[20], amp11[20], amp12[20], amp13[20], amp14[20], amp21[20], amp22[20], amp23[20], amp24[20], amp31[20], amp32[20], amp33[20], amp34[20], amp41[20], amp42[20], amp43[20], amp44[20], amp51[20], amp52[20], amp53[20], amp54[20], amp61[20], amp62[20], amp63[20], amp64[20], amp71[20], amp72[20], amp73[20], amp74[20], amp81[20], amp82[20], amp83[20], amp84[20], amp91[20], amp92[20], amp93[20], amp94[20];
Int_t mult01 = 0, mult02 = 0, mult03 = 0, mult04 = 0, mult11 = 0, mult12 = 0, mult13 = 0, mult14 = 0, mult21 = 0, mult22 = 0, mult23 = 0, mult24 = 0, mult31 = 0, mult32 = 0, mult33 = 0, mult34 = 0, mult41 = 0, mult42 = 0, mult43 = 0, mult44 = 0, mult51 = 0, mult52 = 0, mult53 = 0, mult54 = 0, mult61 = 0, mult62 = 0, mult63 = 0, mult64 = 0, mult71 = 0, mult72 = 0, mult73 = 0, mult74 = 0, mult81 = 0, mult82 = 0, mult83 = 0, mult84 = 0, mult91 = 0, mult92 = 0, mult93 = 0, mult94 = 0;
Int_t cathode_RunNumber, cathode_time, cathode_BunchNumber, cathode_PSpulse;
Double_t cathode_psTime;
Float_t cathode_PulseIntensity;
Int_t first_detn_all[40]; 
Double_t first_tof01[100], first_tof02[100], first_tof03[100], first_tof04[100], first_tof11[100], first_tof12[100], first_tof13[100], first_tof14[100], first_tof21[100], first_tof22[100], first_tof23[100], first_tof24[100], first_tof31[100], first_tof32[100], first_tof33[100], first_tof34[100], first_tof41[100], first_tof42[100], first_tof43[100], first_tof44[100], first_tof51[100], first_tof52[100], first_tof53[100], first_tof54[100], first_tof61[100], first_tof62[100], first_tof63[100], first_tof64[100], first_tof71[100], first_tof72[100], first_tof73[100], first_tof74[100], first_tof81[100], first_tof82[100], first_tof83[100], first_tof84[100], first_tof91[100], first_tof92[100], first_tof93[100], first_tof94[100];
Float_t first_amp01[100], first_amp02[100], first_amp03[100], first_amp04[100], first_amp11[100], first_amp12[100], first_amp13[100], first_amp14[100], first_amp21[100], first_amp22[100], first_amp23[100], first_amp24[100], first_amp31[100], first_amp32[100], first_amp33[100], first_amp34[100], first_amp41[100], first_amp42[100], first_amp43[100], first_amp44[100], first_amp51[100], first_amp52[100], first_amp53[100], first_amp54[100], first_amp61[100], first_amp62[100], first_amp63[100], first_amp64[100], first_amp71[100], first_amp72[100], first_amp73[100], first_amp74[100], first_amp81[100], first_amp82[100], first_amp83[100], first_amp84[100], first_amp91[100], first_amp92[100], first_amp93[100], first_amp94[100];
Int_t first_mult01 = 0, first_mult02 = 0, first_mult03 = 0, first_mult04 = 0, first_mult11 = 0, first_mult12 = 0, first_mult13 = 0, first_mult14 = 0, first_mult21 = 0, first_mult22 = 0, first_mult23 = 0, first_mult24 = 0, first_mult31 = 0, first_mult32 = 0, first_mult33 = 0, first_mult34 = 0, first_mult41 = 0, first_mult42 = 0, first_mult43 = 0, first_mult44 = 0, first_mult51 = 0, first_mult52 = 0, first_mult53 = 0, first_mult54 = 0, first_mult61 = 0, first_mult62 = 0, first_mult63 = 0, first_mult64 = 0, first_mult71 = 0, first_mult72 = 0, first_mult73 = 0, first_mult74 = 0, first_mult81 = 0, first_mult82 = 0, first_mult83 = 0, first_mult84 = 0, first_mult91 = 0, first_mult92 = 0, first_mult93 = 0, first_mult94 = 0;
Double_t neutron_energy;
Bool_t addCoincidences(Int_t the_detn) {
    if(the_detn<0 || the_detn>39) {        
        std::cout << "ERROR in addCoincidences: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(the_detn==0) { mult01++; detn_cathode[0]++; }
    if(the_detn==1) { mult02++; detn_cathode[1]++; }
    if(the_detn==2) { mult03++; detn_cathode[2]++; }
    if(the_detn==3) { mult04++; detn_cathode[3]++; }
    if(the_detn==4) { mult11++; detn_cathode[4]++; }
    if(the_detn==5) { mult12++; detn_cathode[5]++; }
    if(the_detn==6) { mult13++; detn_cathode[6]++; }
    if(the_detn==7) { mult14++; detn_cathode[7]++; }
    if(the_detn==8) { mult21++; detn_cathode[8]++; }
    if(the_detn==9) { mult22++; detn_cathode[9]++; }
    if(the_detn==10) { mult23++; detn_cathode[10]++; }
    if(the_detn==11) { mult24++; detn_cathode[11]++; }
    if(the_detn==12) { mult31++; detn_cathode[12]++; }
    if(the_detn==13) { mult32++; detn_cathode[13]++; }
    if(the_detn==14) { mult33++; detn_cathode[14]++; }
    if(the_detn==15) { mult34++; detn_cathode[15]++; }
    if(the_detn==16) { mult41++; detn_cathode[16]++; }
    if(the_detn==17) { mult42++; detn_cathode[17]++; }
    if(the_detn==18) { mult43++; detn_cathode[18]++; }
    if(the_detn==19) { mult44++; detn_cathode[19]++; }
    if(the_detn==20) { mult51++; detn_cathode[20]++; }
    if(the_detn==21) { mult52++; detn_cathode[21]++; }
    if(the_detn==22) { mult53++; detn_cathode[22]++; }
    if(the_detn==23) { mult54++; detn_cathode[23]++; }
    if(the_detn==24) { mult61++; detn_cathode[24]++; }
    if(the_detn==25) { mult62++; detn_cathode[25]++; }
    if(the_detn==26) { mult63++; detn_cathode[26]++; }
    if(the_detn==27) { mult64++; detn_cathode[27]++; }
    if(the_detn==28) { mult71++; detn_cathode[28]++; }
    if(the_detn==29) { mult72++; detn_cathode[29]++; }
    if(the_detn==30) { mult73++; detn_cathode[30]++; }
    if(the_detn==31) { mult74++; detn_cathode[31]++; }
    if(the_detn==32) { mult81++; detn_cathode[32]++; }
    if(the_detn==33) { mult82++; detn_cathode[33]++; }
    if(the_detn==34) { mult83++; detn_cathode[34]++; }
    if(the_detn==35) { mult84++; detn_cathode[35]++; }
    if(the_detn==36) { mult91++; detn_cathode[36]++; }
    if(the_detn==37) { mult92++; detn_cathode[37]++; }
    if(the_detn==38) { mult93++; detn_cathode[38]++; }
    if(the_detn==39) { mult94++; detn_cathode[39]++; }
    return true;
}

Bool_t addTof(Int_t detn_index, Double_t the_tof) {
    if(detn_cathode[detn_index] < 1 || detn_cathode[detn_index] > 20) {         
        std::cout << "ERROR in addTof: detn_all[" << detn_index << "] value: " << detn_cathode[detn_index] << " not valid" << std::endl;
        return false;
    }
    switch (detn_index) {
        case 0: if (mult01 > 0 && mult01 < 20) tof01[mult01 - 1] = the_tof; break;
        case 1: if (mult02 > 0 && mult02 < 20) tof02[mult02 - 1] = the_tof; break;
        case 2: if (mult03 > 0 && mult03 < 20) tof03[mult03 - 1] = the_tof; break;
        case 3: if (mult04 > 0 && mult04 < 20) tof04[mult04 - 1] = the_tof; break;
        case 4: if (mult11 > 0 && mult11 < 20) tof11[mult11 - 1] = the_tof; break;
        case 5: if (mult12 > 0 && mult12 < 20) tof12[mult12 - 1] = the_tof; break;
        case 6: if (mult13 > 0 && mult13 < 20) tof13[mult13 - 1] = the_tof; break;
        case 7: if (mult14 > 0 && mult14 < 20) tof14[mult14 - 1] = the_tof; break;
        case 8: if (mult21 > 0 && mult21 < 20) tof21[mult21 - 1] = the_tof; break;
        case 9: if (mult22 > 0 && mult22 < 20) tof22[mult22 - 1] = the_tof; break;
        case 10: if (mult23 > 0 && mult23 < 20) tof23[mult23 - 1] = the_tof; break;
        case 11: if (mult24 > 0 && mult24 < 20) tof24[mult24 - 1] = the_tof; break;
        case 12: if (mult31 > 0 && mult31 < 20) tof31[mult31 - 1] = the_tof; break;
        case 13: if (mult32 > 0 && mult32 < 20) tof32[mult32 - 1] = the_tof; break;
        case 14: if (mult33 > 0 && mult33 < 20) tof33[mult33 - 1] = the_tof; break;
        case 15: if (mult34 > 0 && mult34 < 20) tof34[mult34 - 1] = the_tof; break;
        case 16: if (mult41 > 0 && mult41 < 20) tof41[mult41 - 1] = the_tof; break;
        case 17: if (mult42 > 0 && mult42 < 20) tof42[mult42 - 1] = the_tof; break;
        case 18: if (mult43 > 0 && mult43 < 20) tof43[mult43 - 1] = the_tof; break;
        case 19: if (mult44 > 0 && mult44 < 20) tof44[mult44 - 1] = the_tof; break;
        case 20: if (mult51 > 0 && mult51 < 20) tof51[mult51 - 1] = the_tof; break;
        case 21: if (mult52 > 0 && mult52 < 20) tof52[mult52 - 1] = the_tof; break;
        case 22: if (mult53 > 0 && mult53 < 20) tof53[mult53 - 1] = the_tof; break;
        case 23: if (mult54 > 0 && mult54 < 20) tof54[mult54 - 1] = the_tof; break;
        case 24: if (mult61 > 0 && mult61 < 20) tof61[mult61 - 1] = the_tof; break;
        case 25: if (mult62 > 0 && mult62 < 20) tof62[mult62 - 1] = the_tof; break;
        case 26: if (mult63 > 0 && mult63 < 20) tof63[mult63 - 1] = the_tof; break;
        case 27: if (mult64 > 0 && mult64 < 20) tof64[mult64 - 1] = the_tof; break;
        case 28: if (mult71 > 0 && mult71 < 20) tof71[mult71 - 1] = the_tof; break;
        case 29: if (mult72 > 0 && mult72 < 20) tof72[mult72 - 1] = the_tof; break;
        case 30: if (mult73 > 0 && mult73 < 20) tof73[mult73 - 1] = the_tof; break;
        case 31: if (mult74 > 0 && mult74 < 20) tof74[mult74 - 1] = the_tof; break;
        case 32: if (mult81 > 0 && mult81 < 20) tof81[mult81 - 1] = the_tof; break;
        case 33: if (mult82 > 0 && mult82 < 20) tof82[mult82 - 1] = the_tof; break;
        case 34: if (mult83 > 0 && mult83 < 20) tof83[mult83 - 1] = the_tof; break;
        case 35: if (mult84 > 0 && mult84 < 20) tof84[mult84 - 1] = the_tof; break;
        case 36: if (mult91 > 0 && mult91 < 20) tof91[mult91 - 1] = the_tof; break;
        case 37: if (mult92 > 0 && mult92 < 20) tof92[mult92 - 1] = the_tof; break;
        case 38: if (mult93 > 0 && mult93 < 20) tof93[mult93 - 1] = the_tof; break;
        case 39: if (mult94 > 0 && mult94 < 20) tof94[mult94 - 1] = the_tof; break;
        default:
            std::cout << "ERROR in addTof: Unexpected the_detn value: " << detn_index << std::endl;
            return false;
    }
    return true;
}
Bool_t addAmp(Int_t detn_index, Float_t the_amp) {
    if (detn_cathode[detn_index] < 1 || detn_cathode[detn_index] > 20) {         
        std::cout << "ERROR in addAmp: detn_all[" << detn_index << "] value: " << detn_cathode[detn_index] << " not valid" << std::endl;
        return false;
    }
    switch (detn_index) {
        case 0: if (mult01 > 0 && mult01 < 20) amp01[mult01 - 1] = the_amp; break;
        case 1: if (mult02 > 0 && mult02 < 20) amp02[mult02 - 1] = the_amp; break;
        case 2: if (mult03 > 0 && mult03 < 20) amp03[mult03 - 1] = the_amp; break;
        case 3: if (mult04 > 0 && mult04 < 20) amp04[mult04 - 1] = the_amp; break;
        case 4: if (mult11 > 0 && mult11 < 20) amp11[mult11 - 1] = the_amp; break;
        case 5: if (mult12 > 0 && mult12 < 20) amp12[mult12 - 1] = the_amp; break;
        case 6: if (mult13 > 0 && mult13 < 20) amp13[mult13 - 1] = the_amp; break;
        case 7: if (mult14 > 0 && mult14 < 20) amp14[mult14 - 1] = the_amp; break;
        case 8: if (mult21 > 0 && mult21 < 20) amp21[mult21 - 1] = the_amp; break;
        case 9: if (mult22 > 0 && mult22 < 20) amp22[mult22 - 1] = the_amp; break;
        case 10: if (mult23 > 0 && mult23 < 20) amp23[mult23 - 1] = the_amp; break;
        case 11: if (mult24 > 0 && mult24 < 20) amp24[mult24 - 1] = the_amp; break;
        case 12: if (mult31 > 0 && mult31 < 20) amp31[mult31 - 1] = the_amp; break;
        case 13: if (mult32 > 0 && mult32 < 20) amp32[mult32 - 1] = the_amp; break;
        case 14: if (mult33 > 0 && mult33 < 20) amp33[mult33 - 1] = the_amp; break;
        case 15: if (mult34 > 0 && mult34 < 20) amp34[mult34 - 1] = the_amp; break;
        case 16: if (mult41 > 0 && mult41 < 20) amp41[mult41 - 1] = the_amp; break;
        case 17: if (mult42 > 0 && mult42 < 20) amp42[mult42 - 1] = the_amp; break;
        case 18: if (mult43 > 0 && mult43 < 20) amp43[mult43 - 1] = the_amp; break;
        case 19: if (mult44 > 0 && mult44 < 20) amp44[mult44 - 1] = the_amp; break;
        case 20: if (mult51 > 0 && mult51 < 20) amp51[mult51 - 1] = the_amp; break;
        case 21: if (mult52 > 0 && mult52 < 20) amp52[mult52 - 1] = the_amp; break;
        case 22: if (mult53 > 0 && mult53 < 20) amp53[mult53 - 1] = the_amp; break;
        case 23: if (mult54 > 0 && mult54 < 20) amp54[mult54 - 1] = the_amp; break;
        case 24: if (mult61 > 0 && mult61 < 20) amp61[mult61 - 1] = the_amp; break;
        case 25: if (mult62 > 0 && mult62 < 20) amp62[mult62 - 1] = the_amp; break;
        case 26: if (mult63 > 0 && mult63 < 20) amp63[mult63 - 1] = the_amp; break;
        case 27: if (mult64 > 0 && mult64 < 20) amp64[mult64 - 1] = the_amp; break;
        case 28: if (mult71 > 0 && mult71 < 20) amp71[mult71 - 1] = the_amp; break;
        case 29: if (mult72 > 0 && mult72 < 20) amp72[mult72 - 1] = the_amp; break;
        case 30: if (mult73 > 0 && mult73 < 20) amp73[mult73 - 1] = the_amp; break;
        case 31: if (mult74 > 0 && mult74 < 20) amp74[mult74 - 1] = the_amp; break;
        case 32: if (mult81 > 0 && mult81 < 20) amp81[mult81 - 1] = the_amp; break;
        case 33: if (mult82 > 0 && mult82 < 20) amp82[mult82 - 1] = the_amp; break;
        case 34: if (mult83 > 0 && mult83 < 20) amp83[mult83 - 1] = the_amp; break;
        case 35: if (mult84 > 0 && mult84 < 20) amp84[mult84 - 1] = the_amp; break;
        case 36: if (mult91 > 0 && mult91 < 20) amp91[mult91 - 1] = the_amp; break;
        case 37: if (mult92 > 0 && mult92 < 20) amp92[mult92 - 1] = the_amp; break;
        case 38: if (mult93 > 0 && mult93 < 20) amp93[mult93 - 1] = the_amp; break;
        case 39: if (mult94 > 0 && mult94 < 20) amp94[mult94 - 1] = the_amp; break;
        default:
            std::cout << "ERROR in addAmp: Unexpected detn_index value: " << detn_index << std::endl;
            return false;
    }
    return true;
}

void anode_cathode(Int_t run_number) {
    std::string filename = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/anodes_final/out_run" + std::to_string(run_number) + ".root";
    TFile *file = TFile::Open(filename.c_str()) ;//coincidences file
    TTree *intree = (TTree*)file->Get("coincidences");
    std::string filename2 = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Cathodes/collection_cathode" + std::to_string(run_number) + ".root";
    TFile *file2 = TFile::Open(filename2.c_str());// data from the pickup and the cathode
    TTree *cathodeTree = (TTree*)file2->Get("cathode_preselection");


    Long64_t nentries = intree->GetEntries(); // getting the number of entries to loop over them
    Long64_t ncathodeEntries = cathodeTree->GetEntries();
    std::cout << "INFO: Processing run_number " << run_number << " with " << nentries << " entries" << std::endl;
    std::cout << "INFO: Processing run_number " << run_number << " with " << ncathodeEntries << " entries" << std::endl;
    //the first thing we are going to do (after loading the branches) is to callibrate the time of the cathodes by the time of the pickup.
    //then we are going to loop over the entries of the coincidences tree and collect for each event the detectors for with detn_all ==1.
    // by making a switch we will select the corresponding tof and amp variables of the anodes and open a coincidence window with the cathodes. 
    //we will select cathode events that are between t0 and t0+300 being t0 the time of the anode. we will impose a threshold in amplitudes of 0.3*amp0.
    
    
    intree->SetBranchStatus("*", 1); // turn on all branches
    intree->SetBranchAddress("RunNumber", &RunNumber);
    intree->SetBranchAddress("time", &eventTime);
    intree->SetBranchAddress("psTime", &psTime);
    intree->SetBranchAddress("BunchNumber", &BunchNumber);
    intree->SetBranchAddress("PSpulse", &PSpulse);
    intree->SetBranchAddress("PulseIntensity", &PulseIntensity);
    intree->SetBranchAddress("detn", &detn_all);
    intree->SetBranchAddress("mult0", &mult0);
    intree->SetBranchAddress("mult1", &mult1);
    intree->SetBranchAddress("mult2", &mult2);
    intree->SetBranchAddress("mult3", &mult3);
    intree->SetBranchAddress("mult4", &mult4);
    intree->SetBranchAddress("mult5", &mult5);
    intree->SetBranchAddress("mult6", &mult6);
    intree->SetBranchAddress("mult7", &mult7);
    intree->SetBranchAddress("mult8", &mult8);
    intree->SetBranchAddress("mult9", &mult9);
    intree->SetBranchAddress("tof0", tof0);
    intree->SetBranchAddress("tof1", tof1);
    intree->SetBranchAddress("tof2", tof2);
    intree->SetBranchAddress("tof3", tof3);
    intree->SetBranchAddress("tof4", tof4);
    intree->SetBranchAddress("tof5", tof5);
    intree->SetBranchAddress("tof6", tof6);
    intree->SetBranchAddress("tof7", tof7);
    intree->SetBranchAddress("tof8", tof8);
    intree->SetBranchAddress("tof9", tof9);
    intree->SetBranchAddress("amp0", amp0);
    intree->SetBranchAddress("amp1", amp1);
    intree->SetBranchAddress("amp2", amp2);
    intree->SetBranchAddress("amp3", amp3);
    intree->SetBranchAddress("amp4", amp4);
    intree->SetBranchAddress("amp5", amp5);
    intree->SetBranchAddress("amp6", amp6);
    intree->SetBranchAddress("amp7", amp7);
    intree->SetBranchAddress("amp8", amp8);
    intree->SetBranchAddress("amp9", amp9);
    intree->SetBranchAddress("mult", &mult);
    intree->SetBranchAddress("neutron_energy", &neutron_energy);
    intree->SetBranchAddress("gamma_flash0", &gamma_flash0);
    intree->SetBranchAddress("gamma_flash1", &gamma_flash1);
    intree->SetBranchAddress("gamma_flash2", &gamma_flash2);
    intree->SetBranchAddress("gamma_flash3", &gamma_flash3);
    intree->SetBranchAddress("gamma_flash4", &gamma_flash4);
    intree->SetBranchAddress("gamma_flash5", &gamma_flash5);
    intree->SetBranchAddress("gamma_flash6", &gamma_flash6);
    intree->SetBranchAddress("gamma_flash7", &gamma_flash7);
    intree->SetBranchAddress("gamma_flash8", &gamma_flash8);
    intree->SetBranchAddress("gamma_flash9", &gamma_flash9);


    cathodeTree->SetBranchStatus("*", 1); // turn on all branches
    cathodeTree->SetBranchAddress("RunNumber", &cathode_RunNumber);
    cathodeTree->SetBranchAddress("time", &cathode_time);
    cathodeTree->SetBranchAddress("psTime", &cathode_psTime);
    cathodeTree->SetBranchAddress("BunchNumber", &cathode_BunchNumber);
    cathodeTree->SetBranchAddress("PSpulse", &cathode_PSpulse);
    cathodeTree->SetBranchAddress("PulseIntensity", &cathode_PulseIntensity);
    cathodeTree->SetBranchAddress("first_tof01", &first_tof01);
    cathodeTree->SetBranchAddress("first_tof02", &first_tof02);
    cathodeTree->SetBranchAddress("first_tof03", &first_tof03);
    cathodeTree->SetBranchAddress("first_tof04", &first_tof04);
    cathodeTree->SetBranchAddress("first_tof11", &first_tof11);
    cathodeTree->SetBranchAddress("first_tof12", &first_tof12);
    cathodeTree->SetBranchAddress("first_tof13", &first_tof13);
    cathodeTree->SetBranchAddress("first_tof14", &first_tof14);
    cathodeTree->SetBranchAddress("first_tof21", &first_tof21);
    cathodeTree->SetBranchAddress("first_tof22", &first_tof22);
    cathodeTree->SetBranchAddress("first_tof23", &first_tof23);
    cathodeTree->SetBranchAddress("first_tof24", &first_tof24);
    cathodeTree->SetBranchAddress("first_tof31", &first_tof31);
    cathodeTree->SetBranchAddress("first_tof32", &first_tof32);
    cathodeTree->SetBranchAddress("first_tof33", &first_tof33);
    cathodeTree->SetBranchAddress("first_tof34", &first_tof34);
    cathodeTree->SetBranchAddress("first_tof41", &first_tof41);
    cathodeTree->SetBranchAddress("first_tof42", &first_tof42);
    cathodeTree->SetBranchAddress("first_tof43", &first_tof43);
    cathodeTree->SetBranchAddress("first_tof44", &first_tof44);
    cathodeTree->SetBranchAddress("first_tof51", &first_tof51);
    cathodeTree->SetBranchAddress("first_tof52", &first_tof52);
    cathodeTree->SetBranchAddress("first_tof53", &first_tof53);
    cathodeTree->SetBranchAddress("first_tof54", &first_tof54);
    cathodeTree->SetBranchAddress("first_tof61", &first_tof61);
    cathodeTree->SetBranchAddress("first_tof62", &first_tof62);
    cathodeTree->SetBranchAddress("first_tof63", &first_tof63);
    cathodeTree->SetBranchAddress("first_tof64", &first_tof64);
    cathodeTree->SetBranchAddress("first_tof71", &first_tof71);
    cathodeTree->SetBranchAddress("first_tof72", &first_tof72);
    cathodeTree->SetBranchAddress("first_tof73", &first_tof73);
    cathodeTree->SetBranchAddress("first_tof74", &first_tof74);
    cathodeTree->SetBranchAddress("first_tof81", &first_tof81);
    cathodeTree->SetBranchAddress("first_tof82", &first_tof82);
    cathodeTree->SetBranchAddress("first_tof83", &first_tof83);
    cathodeTree->SetBranchAddress("first_tof84", &first_tof84);
    cathodeTree->SetBranchAddress("first_tof91", &first_tof91);
    cathodeTree->SetBranchAddress("first_tof92", &first_tof92);
    cathodeTree->SetBranchAddress("first_tof93", &first_tof93);
    cathodeTree->SetBranchAddress("first_tof94", &first_tof94);

    cathodeTree->SetBranchAddress("first_amp01", &first_amp01);
    cathodeTree->SetBranchAddress("first_amp02", &first_amp02);
    cathodeTree->SetBranchAddress("first_amp03", &first_amp03);
    cathodeTree->SetBranchAddress("first_amp04", &first_amp04);
    cathodeTree->SetBranchAddress("first_amp11", &first_amp11);
    cathodeTree->SetBranchAddress("first_amp12", &first_amp12);
    cathodeTree->SetBranchAddress("first_amp13", &first_amp13);
    cathodeTree->SetBranchAddress("first_amp14", &first_amp14);
    cathodeTree->SetBranchAddress("first_amp21", &first_amp21);
    cathodeTree->SetBranchAddress("first_amp22", &first_amp22);
    cathodeTree->SetBranchAddress("first_amp23", &first_amp23);
    cathodeTree->SetBranchAddress("first_amp24", &first_amp24);
    cathodeTree->SetBranchAddress("first_amp31", &first_amp31);
    cathodeTree->SetBranchAddress("first_amp32", &first_amp32);
    cathodeTree->SetBranchAddress("first_amp33", &first_amp33);
    cathodeTree->SetBranchAddress("first_amp34", &first_amp34);
    cathodeTree->SetBranchAddress("first_amp41", &first_amp41);
    cathodeTree->SetBranchAddress("first_amp42", &first_amp42);
    cathodeTree->SetBranchAddress("first_amp43", &first_amp43);
    cathodeTree->SetBranchAddress("first_amp44", &first_amp44);
    cathodeTree->SetBranchAddress("first_amp51", &first_amp51);
    cathodeTree->SetBranchAddress("first_amp52", &first_amp52);
    cathodeTree->SetBranchAddress("first_amp53", &first_amp53);
    cathodeTree->SetBranchAddress("first_amp54", &first_amp54);
    cathodeTree->SetBranchAddress("first_amp61", &first_amp61);
    cathodeTree->SetBranchAddress("first_amp62", &first_amp62);
    cathodeTree->SetBranchAddress("first_amp63", &first_amp63);
    cathodeTree->SetBranchAddress("first_amp64", &first_amp64);
    cathodeTree->SetBranchAddress("first_amp71", &first_amp71);
    cathodeTree->SetBranchAddress("first_amp72", &first_amp72);
    cathodeTree->SetBranchAddress("first_amp73", &first_amp73);
    cathodeTree->SetBranchAddress("first_amp74", &first_amp74);
    cathodeTree->SetBranchAddress("first_amp81", &first_amp81);
    cathodeTree->SetBranchAddress("first_amp82", &first_amp82);
    cathodeTree->SetBranchAddress("first_amp83", &first_amp83);
    cathodeTree->SetBranchAddress("first_amp84", &first_amp84);
    cathodeTree->SetBranchAddress("first_amp91", &first_amp91);
    cathodeTree->SetBranchAddress("first_amp92", &first_amp92);
    cathodeTree->SetBranchAddress("first_amp93", &first_amp93);
    cathodeTree->SetBranchAddress("first_amp94", &first_amp94);

    cathodeTree->SetBranchAddress("first_mult01", &first_mult01);
    cathodeTree->SetBranchAddress("first_mult02", &first_mult02);
    cathodeTree->SetBranchAddress("first_mult03", &first_mult03);
    cathodeTree->SetBranchAddress("first_mult04", &first_mult04);
    cathodeTree->SetBranchAddress("first_mult11", &first_mult11);
    cathodeTree->SetBranchAddress("first_mult12", &first_mult12);
    cathodeTree->SetBranchAddress("first_mult13", &first_mult13);
    cathodeTree->SetBranchAddress("first_mult14", &first_mult14);
    cathodeTree->SetBranchAddress("first_mult21", &first_mult21);
    cathodeTree->SetBranchAddress("first_mult22", &first_mult22);
    cathodeTree->SetBranchAddress("first_mult23", &first_mult23);
    cathodeTree->SetBranchAddress("first_mult24", &first_mult24);
    cathodeTree->SetBranchAddress("first_mult31", &first_mult31);
    cathodeTree->SetBranchAddress("first_mult32", &first_mult32);
    cathodeTree->SetBranchAddress("first_mult33", &first_mult33);
    cathodeTree->SetBranchAddress("first_mult34", &first_mult34);
    cathodeTree->SetBranchAddress("first_mult41", &first_mult41);
    cathodeTree->SetBranchAddress("first_mult42", &first_mult42);
    cathodeTree->SetBranchAddress("first_mult43", &first_mult43);
    cathodeTree->SetBranchAddress("first_mult44", &first_mult44);
    cathodeTree->SetBranchAddress("first_mult51", &first_mult51);
    cathodeTree->SetBranchAddress("first_mult52", &first_mult52);
    cathodeTree->SetBranchAddress("first_mult53", &first_mult53);
    cathodeTree->SetBranchAddress("first_mult54", &first_mult54);
    cathodeTree->SetBranchAddress("first_mult61", &first_mult61);
    cathodeTree->SetBranchAddress("first_mult62", &first_mult62);
    cathodeTree->SetBranchAddress("first_mult63", &first_mult63);
    cathodeTree->SetBranchAddress("first_mult64", &first_mult64);
    cathodeTree->SetBranchAddress("first_mult71", &first_mult71);
    cathodeTree->SetBranchAddress("first_mult72", &first_mult72);
    cathodeTree->SetBranchAddress("first_mult73", &first_mult73);
    cathodeTree->SetBranchAddress("first_mult74", &first_mult74);
    cathodeTree->SetBranchAddress("first_mult81", &first_mult81);
    cathodeTree->SetBranchAddress("first_mult82", &first_mult82);
    cathodeTree->SetBranchAddress("first_mult83", &first_mult83);
    cathodeTree->SetBranchAddress("first_mult84", &first_mult84);
    cathodeTree->SetBranchAddress("first_mult91", &first_mult91);
    cathodeTree->SetBranchAddress("first_mult92", &first_mult92);
    cathodeTree->SetBranchAddress("first_mult93", &first_mult93);
    cathodeTree->SetBranchAddress("first_mult94", &first_mult94);
    cathodeTree->SetBranchAddress("detn_cathode", &first_detn_all);

    // Create output file and tree
    std::string out_filename = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Cathodes/cathode_coincidences" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(out_filename.c_str(), "RECREATE");
    TTree *outtree = new TTree("cathode_coincidences", "cathode_coincidences");
    outtree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    outtree->Branch("BunchNumber", &BunchNumber, "BunchNumber/I");
    outtree->Branch("PSpulse", &PSpulse, "PSpulse/I");
    outtree->Branch("time", &eventTime, "time/I");
    outtree->Branch("psTime", &psTime, "psTime/D");
    outtree->Branch("PulseIntensity", &PulseIntensity, "PulseIntensity/F");
    outtree->Branch("cathode_detn", &detn_cathode, "detn_cathode[40]/I");
    outtree->Branch("mult0", &mult0, "mult0/I");
    outtree->Branch("mult1", &mult1, "mult1/I");
    outtree->Branch("mult2", &mult2, "mult2/I");
    outtree->Branch("mult3", &mult3, "mult3/I");
    outtree->Branch("mult4", &mult4, "mult4/I");
    outtree->Branch("mult5", &mult5, "mult5/I");
    outtree->Branch("mult6", &mult6, "mult6/I");
    outtree->Branch("mult7", &mult7, "mult7/I");
    outtree->Branch("mult8", &mult8, "mult8/I");
    outtree->Branch("mult9", &mult9, "mult9/I");
    outtree->Branch("amp0", amp0, "amp0[mult0]/F");
    outtree->Branch("amp1", amp1, "amp1[mult1]/F");
    outtree->Branch("amp2", amp2, "amp2[mult2]/F");
    outtree->Branch("amp3", amp3, "amp3[mult3]/F");
    outtree->Branch("amp4", amp4, "amp4[mult4]/F");
    outtree->Branch("amp5", amp5, "amp5[mult5]/F");
    outtree->Branch("amp6", amp6, "amp6[mult6]/F");
    outtree->Branch("amp7", amp7, "amp7[mult7]/F");
    outtree->Branch("amp8", amp8, "amp8[mult8]/F");
    outtree->Branch("amp9", amp9, "amp9[mult9]/F");
    outtree->Branch("detn", &detn_all, "detn[10]/I");
    outtree->Branch("mult", &mult, "mult/I");
    outtree->Branch("tof0", tof0, "tof0[mult0]/D");
    outtree->Branch("tof1", tof1, "tof1[mult1]/D");
    outtree->Branch("tof2", tof2, "tof2[mult2]/D");
    outtree->Branch("tof3", tof3, "tof3[mult3]/D");
    outtree->Branch("tof4", tof4, "tof4[mult4]/D");
    outtree->Branch("tof5", tof5, "tof5[mult5]/D");
    outtree->Branch("tof6", tof6, "tof6[mult6]/D");
    outtree->Branch("tof7", tof7, "tof7[mult7]/D");
    outtree->Branch("tof8", tof8, "tof8[mult8]/D");
    outtree->Branch("tof9", tof9, "tof9[mult9]/D");
    outtree->Branch("mult01", &mult01, "mult01/I");
    outtree->Branch("mult02", &mult02, "mult02/I");
    outtree->Branch("mult03", &mult03, "mult03/I");
    outtree->Branch("mult04", &mult04, "mult04/I");
    outtree->Branch("mult11", &mult11, "mult11/I");
    outtree->Branch("mult12", &mult12, "mult12/I");
    outtree->Branch("mult13", &mult13, "mult13/I");
    outtree->Branch("mult14", &mult14, "mult14/I");
    outtree->Branch("mult21", &mult21, "mult21/I");
    outtree->Branch("mult22", &mult22, "mult22/I");
    outtree->Branch("mult23", &mult23, "mult23/I");
    outtree->Branch("mult24", &mult24, "mult24/I");
    outtree->Branch("mult31", &mult31, "mult31/I");
    outtree->Branch("mult32", &mult32, "mult32/I");
    outtree->Branch("mult33", &mult33, "mult33/I");
    outtree->Branch("mult34", &mult34, "mult34/I");
    outtree->Branch("mult41", &mult41, "mult41/I");
    outtree->Branch("mult42", &mult42, "mult42/I");
    outtree->Branch("mult43", &mult43, "mult43/I");
    outtree->Branch("mult44", &mult44, "mult44/I");
    outtree->Branch("mult51", &mult51, "mult51/I");
    outtree->Branch("mult52", &mult52, "mult52/I");
    outtree->Branch("mult53", &mult53, "mult53/I");
    outtree->Branch("mult54", &mult54, "mult54/I");
    outtree->Branch("mult61", &mult61, "mult61/I");
    outtree->Branch("mult62", &mult62, "mult62/I");
    outtree->Branch("mult63", &mult63, "mult63/I");
    outtree->Branch("mult64", &mult64, "mult64/I");
    outtree->Branch("mult71", &mult71, "mult71/I");
    outtree->Branch("mult72", &mult72, "mult72/I");
    outtree->Branch("mult73", &mult73, "mult73/I");
    outtree->Branch("mult74", &mult74, "mult74/I");
    outtree->Branch("mult81", &mult81, "mult81/I");
    outtree->Branch("mult82", &mult82, "mult82/I");
    outtree->Branch("mult83", &mult83, "mult83/I");
    outtree->Branch("mult84", &mult84, "mult84/I");
    outtree->Branch("mult91", &mult91, "mult91/I");
    outtree->Branch("mult92", &mult92, "mult92/I");
    outtree->Branch("mult93", &mult93, "mult93/I");
    outtree->Branch("mult94", &mult94, "mult94/I");
    outtree->Branch("tof01", tof01, "tof01[mult01]/D");
    outtree->Branch("tof02", tof02, "tof02[mult02]/D");
    outtree->Branch("tof03", tof03, "tof03[mult03]/D");
    outtree->Branch("tof04", tof04, "tof04[mult04]/D");
    outtree->Branch("tof11", tof11, "tof11[mult11]/D");
    outtree->Branch("tof12", tof12, "tof12[mult12]/D");
    outtree->Branch("tof13", tof13, "tof13[mult13]/D");
    outtree->Branch("tof14", tof14, "tof14[mult14]/D");
    outtree->Branch("tof21", tof21, "tof21[mult21]/D");
    outtree->Branch("tof22", tof22, "tof22[mult22]/D");
    outtree->Branch("tof23", tof23, "tof23[mult23]/D");
    outtree->Branch("tof24", tof24, "tof24[mult24]/D");
    outtree->Branch("tof31", tof31, "tof31[mult31]/D");
    outtree->Branch("tof32", tof32, "tof32[mult32]/D");
    outtree->Branch("tof33", tof33, "tof33[mult33]/D");
    outtree->Branch("tof34", tof34, "tof34[mult34]/D");
    outtree->Branch("tof41", tof41, "tof41[mult41]/D");
    outtree->Branch("tof42", tof42, "tof42[mult42]/D");
    outtree->Branch("tof43", tof43, "tof43[mult43]/D");
    outtree->Branch("tof44", tof44, "tof44[mult44]/D");
    outtree->Branch("tof51", tof51, "tof51[mult51]/D");
    outtree->Branch("tof52", tof52, "tof52[mult52]/D");
    outtree->Branch("tof53", tof53, "tof53[mult53]/D");
    outtree->Branch("tof54", tof54, "tof54[mult54]/D");
    outtree->Branch("tof61", tof61, "tof61[mult61]/D");
    outtree->Branch("tof62", tof62, "tof62[mult62]/D");
    outtree->Branch("tof63", tof63, "tof63[mult63]/D");
    outtree->Branch("tof64", tof64, "tof64[mult64]/D");
    outtree->Branch("tof71", tof71, "tof71[mult71]/D");
    outtree->Branch("tof72", tof72, "tof72[mult72]/D");
    outtree->Branch("tof73", tof73, "tof73[mult73]/D");
    outtree->Branch("tof74", tof74, "tof74[mult74]/D");
    outtree->Branch("tof81", tof81, "tof81[mult81]/D");
    outtree->Branch("tof82", tof82, "tof82[mult82]/D");
    outtree->Branch("tof83", tof83, "tof83[mult83]/D");
    outtree->Branch("tof84", tof84, "tof84[mult84]/D");
    outtree->Branch("tof91", tof91, "tof91[mult91]/D");
    outtree->Branch("tof92", tof92, "tof92[mult92]/D");
    outtree->Branch("tof93", tof93, "tof93[mult93]/D");
    outtree->Branch("tof94", tof94, "tof94[mult94]/D");
    outtree->Branch("amp01", amp01, "amp01[mult01]/F");
    outtree->Branch("amp02", amp02, "amp02[mult02]/F");
    outtree->Branch("amp03", amp03, "amp03[mult03]/F");
    outtree->Branch("amp04", amp04, "amp04[mult04]/F");
    outtree->Branch("amp11", amp11, "amp11[mult11]/F");
    outtree->Branch("amp12", amp12, "amp12[mult12]/F");
    outtree->Branch("amp13", amp13, "amp13[mult13]/F");
    outtree->Branch("amp14", amp14, "amp14[mult14]/F");
    outtree->Branch("amp21", amp21, "amp21[mult21]/F");
    outtree->Branch("amp22", amp22, "amp22[mult22]/F");
    outtree->Branch("amp23", amp23, "amp23[mult23]/F");
    outtree->Branch("amp24", amp24, "amp24[mult24]/F");
    outtree->Branch("amp31", amp31, "amp31[mult31]/F");
    outtree->Branch("amp32", amp32, "amp32[mult32]/F");
    outtree->Branch("amp33", amp33, "amp33[mult33]/F");
    outtree->Branch("amp34", amp34, "amp34[mult34]/F");
    outtree->Branch("amp41", amp41, "amp41[mult41]/F");
    outtree->Branch("amp42", amp42, "amp42[mult42]/F");
    outtree->Branch("amp43", amp43, "amp43[mult43]/F");
    outtree->Branch("amp44", amp44, "amp44[mult44]/F");
    outtree->Branch("amp51", amp51, "amp51[mult51]/F");
    outtree->Branch("amp52", amp52, "amp52[mult52]/F");
    outtree->Branch("amp53", amp53, "amp53[mult53]/F");
    outtree->Branch("amp54", amp54, "amp54[mult54]/F");
    outtree->Branch("amp61", amp61, "amp61[mult61]/F");
    outtree->Branch("amp62", amp62, "amp62[mult62]/F");
    outtree->Branch("amp63", amp63, "amp63[mult63]/F");
    outtree->Branch("amp64", amp64, "amp64[mult64]/F");
    outtree->Branch("amp71", amp71, "amp71[mult71]/F");
    outtree->Branch("amp72", amp72, "amp72[mult72]/F");
    outtree->Branch("amp73", amp73, "amp73[mult73]/F");
    outtree->Branch("amp74", amp74, "amp74[mult74]/F");
    outtree->Branch("amp81", amp81, "amp81[mult81]/F");
    outtree->Branch("amp82", amp82, "amp82[mult82]/F");
    outtree->Branch("amp83", amp83, "amp83[mult83]/F");
    outtree->Branch("amp84", amp84, "amp84[mult84]/F");
    outtree->Branch("amp91", amp91, "amp91[mult91]/F");
    outtree->Branch("amp92", amp92, "amp92[mult92]/F");
    outtree->Branch("amp93", amp93, "amp93[mult93]/F");
    outtree->Branch("amp94", amp94, "amp94[mult94]/F");
    outtree->Branch("neutron_energy", &neutron_energy, "neutron_energy/D");
    for (Long64_t i = 0; i < nentries; i++) {
        intree->GetEntry(i);
        std::cout << "Event: " << i << std::endl;
        for (Long64_t j = 0; j < ncathodeEntries; j++) {
            cathodeTree->GetEntry(j);
            if ( cathode_BunchNumber== BunchNumber && cathode_PSpulse==PSpulse && cathode_RunNumber==RunNumber){
                if (cathode_psTime== psTime && cathode_time==eventTime && cathode_PulseIntensity==PulseIntensity) {
                    
                    for (int k = 0; k < 10; k++) { //look for coincidences in the 10 detectors imposing that the anode has triggered (detn_all[k] != 0)
                        if (detn_all[k] != 0) {
                            Double_t *anode_tof = nullptr; // pointer to hold tof, mult and amp values for the anode
                            Float_t *anode_amp = nullptr;
                            Int_t anode_mult;
                            switch (k) {
                            case 0: anode_tof = new Double_t[mult0]; anode_amp = new Float_t[mult0]; std::copy(tof0, tof0 + mult0, anode_tof); std::copy(amp0, amp0 + mult0, anode_amp); anode_mult = mult0; break;
                            case 1: anode_tof = new Double_t[mult1]; anode_amp = new Float_t[mult1]; std::copy(tof1, tof1 + mult1, anode_tof); std::copy(amp1, amp1 + mult1, anode_amp); anode_mult = mult1; break;
                            case 2: anode_tof = new Double_t[mult2]; anode_amp = new Float_t[mult2]; std::copy(tof2, tof2 + mult2, anode_tof); std::copy(amp2, amp2 + mult2, anode_amp); anode_mult = mult2; break;
                            case 3: anode_tof = new Double_t[mult3]; anode_amp = new Float_t[mult3]; std::copy(tof3, tof3 + mult3, anode_tof); std::copy(amp3, amp3 + mult3, anode_amp); anode_mult = mult3; break;
                            case 4: anode_tof = new Double_t[mult4]; anode_amp = new Float_t[mult4]; std::copy(tof4, tof4 + mult4, anode_tof); std::copy(amp4, amp4 + mult4, anode_amp); anode_mult = mult4; break;
                            case 5: anode_tof = new Double_t[mult5]; anode_amp = new Float_t[mult5]; std::copy(tof5, tof5 + mult5, anode_tof); std::copy(amp5, amp5 + mult5, anode_amp); anode_mult = mult5; break;
                            case 6: anode_tof = new Double_t[mult6]; anode_amp = new Float_t[mult6]; std::copy(tof6, tof6 + mult6, anode_tof); std::copy(amp6, amp6 + mult6, anode_amp); anode_mult = mult6; break;
                            case 7: anode_tof = new Double_t[mult7]; anode_amp = new Float_t[mult7]; std::copy(tof7, tof7 + mult7, anode_tof); std::copy(amp7, amp7 + mult7, anode_amp); anode_mult = mult7; break;
                            case 8: anode_tof = new Double_t[mult8]; anode_amp = new Float_t[mult8]; std::copy(tof8, tof8 + mult8, anode_tof); std::copy(amp8, amp8 + mult8, anode_amp); anode_mult = mult8; break;
                            case 9: anode_tof = new Double_t[mult9]; anode_amp = new Float_t[mult9]; std::copy(tof9, tof9 + mult9, anode_tof); std::copy(amp9, amp9 + mult9, anode_amp); anode_mult = mult9; break;
                            }
                            for (int l = 0; l < anode_mult; l++) { // loop over anode hits
                                for (int subdet = 1; subdet < 5; subdet++) {
                                    Double_t *first_tof = nullptr; // pointer to hold tof, mult and amp values for the cathode detector
                                    Int_t first_mult = 0;
                                    Float_t *first_amp = nullptr;
                                    Int_t cathode_detn_index = -1;
                                    if (k == 0) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult01]; first_amp = new Float_t[first_mult01]; std::copy(first_tof01, first_tof01 + first_mult01, first_tof); std::copy(first_amp01, first_amp01 + first_mult01, first_amp); first_mult = first_mult01; cathode_detn_index = 0; break;
                                            case 2: first_tof = new Double_t[first_mult02]; first_amp = new Float_t[first_mult02]; std::copy(first_tof02, first_tof02 + first_mult02, first_tof); std::copy(first_amp02, first_amp02 + first_mult02, first_amp); first_mult = first_mult02; cathode_detn_index = 1; break;
                                            case 3: first_tof = new Double_t[first_mult03]; first_amp = new Float_t[first_mult03]; std::copy(first_tof03, first_tof03 + first_mult03, first_tof); std::copy(first_amp03, first_amp03 + first_mult03, first_amp); first_mult = first_mult03; cathode_detn_index = 2; break;
                                            case 4: first_tof = new Double_t[first_mult04]; first_amp = new Float_t[first_mult04]; std::copy(first_tof04, first_tof04 + first_mult04, first_tof); std::copy(first_amp04, first_amp04 + first_mult04, first_amp); first_mult = first_mult04; cathode_detn_index = 3; break;
                                        }
                                    } else if (k == 1) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult11]; first_amp = new Float_t[first_mult11]; std::copy(first_tof11, first_tof11 + first_mult11, first_tof); std::copy(first_amp11, first_amp11 + first_mult11, first_amp); first_mult = first_mult11; cathode_detn_index = 4; break;
                                            case 2: first_tof = new Double_t[first_mult12]; first_amp = new Float_t[first_mult12]; std::copy(first_tof12, first_tof12 + first_mult12, first_tof); std::copy(first_amp12, first_amp12 + first_mult12, first_amp); first_mult = first_mult12; cathode_detn_index = 5; break;
                                            case 3: first_tof = new Double_t[first_mult13]; first_amp = new Float_t[first_mult13]; std::copy(first_tof13, first_tof13 + first_mult13, first_tof); std::copy(first_amp13, first_amp13 + first_mult13, first_amp); first_mult = first_mult13; cathode_detn_index = 6; break;
                                            case 4: first_tof = new Double_t[first_mult14]; first_amp = new Float_t[first_mult14]; std::copy(first_tof14, first_tof14 + first_mult14, first_tof); std::copy(first_amp14, first_amp14 + first_mult14, first_amp); first_mult = first_mult14; cathode_detn_index = 7; break;
                                        }
                                    } else if (k == 2) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult21]; first_amp = new Float_t[first_mult21]; std::copy(first_tof21, first_tof21 + first_mult21, first_tof); std::copy(first_amp21, first_amp21 + first_mult21, first_amp); first_mult = first_mult21; cathode_detn_index = 8; break;
                                            case 2: first_tof = new Double_t[first_mult22]; first_amp = new Float_t[first_mult22]; std::copy(first_tof22, first_tof22 + first_mult22, first_tof); std::copy(first_amp22, first_amp22 + first_mult22, first_amp); first_mult = first_mult22; cathode_detn_index = 9; break;
                                            case 3: first_tof = new Double_t[first_mult23]; first_amp = new Float_t[first_mult23]; std::copy(first_tof23, first_tof23 + first_mult23, first_tof); std::copy(first_amp23, first_amp23 + first_mult23, first_amp); first_mult = first_mult23; cathode_detn_index = 10; break;
                                            case 4: first_tof = new Double_t[first_mult24]; first_amp = new Float_t[first_mult24]; std::copy(first_tof24, first_tof24 + first_mult24, first_tof); std::copy(first_amp24, first_amp24 + first_mult24, first_amp); first_mult = first_mult24; cathode_detn_index = 11; break;
                                        }
                                    } else if (k == 3) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult31]; first_amp = new Float_t[first_mult31]; std::copy(first_tof31, first_tof31 + first_mult31, first_tof); std::copy(first_amp31, first_amp31 + first_mult31, first_amp); first_mult = first_mult31; cathode_detn_index = 12; break;
                                            case 2: first_tof = new Double_t[first_mult32]; first_amp = new Float_t[first_mult32]; std::copy(first_tof32, first_tof32 + first_mult32, first_tof); std::copy(first_amp32, first_amp32 + first_mult32, first_amp); first_mult = first_mult32; cathode_detn_index = 13; break;
                                            case 3: first_tof = new Double_t[first_mult33]; first_amp = new Float_t[first_mult33]; std::copy(first_tof33, first_tof33 + first_mult33, first_tof); std::copy(first_amp33, first_amp33 + first_mult33, first_amp); first_mult = first_mult33; cathode_detn_index = 14; break;
                                            case 4: first_tof = new Double_t[first_mult34]; first_amp = new Float_t[first_mult34]; std::copy(first_tof34, first_tof34 + first_mult34, first_tof); std::copy(first_amp34, first_amp34 + first_mult34, first_amp); first_mult = first_mult34; cathode_detn_index = 15; break;
                                        }
                                    } else if (k == 4) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult41]; first_amp = new Float_t[first_mult41]; std::copy(first_tof41, first_tof41 + first_mult41, first_tof); std::copy(first_amp41, first_amp41 + first_mult41, first_amp); first_mult = first_mult41; cathode_detn_index = 16; break;
                                            case 2: first_tof = new Double_t[first_mult42]; first_amp = new Float_t[first_mult42]; std::copy(first_tof42, first_tof42 + first_mult42, first_tof); std::copy(first_amp42, first_amp42 + first_mult42, first_amp); first_mult = first_mult42; cathode_detn_index = 17; break;
                                            case 3: first_tof = new Double_t[first_mult43]; first_amp = new Float_t[first_mult43]; std::copy(first_tof43, first_tof43 + first_mult43, first_tof); std::copy(first_amp43, first_amp43 + first_mult43, first_amp); first_mult = first_mult43; cathode_detn_index = 18; break;
                                            case 4: first_tof = new Double_t[first_mult44]; first_amp = new Float_t[first_mult44]; std::copy(first_tof44, first_tof44 + first_mult44, first_tof); std::copy(first_amp44, first_amp44 + first_mult44, first_amp); first_mult = first_mult44; cathode_detn_index = 19; break;
                                        }
                                    } else if (k == 5) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult51]; first_amp = new Float_t[first_mult51]; std::copy(first_tof51, first_tof51 + first_mult51, first_tof); std::copy(first_amp51, first_amp51 + first_mult51, first_amp); first_mult = first_mult51; cathode_detn_index = 20; break;
                                            case 2: first_tof = new Double_t[first_mult52]; first_amp = new Float_t[first_mult52]; std::copy(first_tof52, first_tof52 + first_mult52, first_tof); std::copy(first_amp52, first_amp52 + first_mult52, first_amp); first_mult = first_mult52; cathode_detn_index = 21; break;
                                            case 3: first_tof = new Double_t[first_mult53]; first_amp = new Float_t[first_mult53]; std::copy(first_tof53, first_tof53 + first_mult53, first_tof); std::copy(first_amp53, first_amp53 + first_mult53, first_amp); first_mult = first_mult53; cathode_detn_index = 22; break;
                                            case 4: first_tof = new Double_t[first_mult54]; first_amp = new Float_t[first_mult54]; std::copy(first_tof54, first_tof54 + first_mult54, first_tof); std::copy(first_amp54, first_amp54 + first_mult54, first_amp); first_mult = first_mult54; cathode_detn_index = 23; break;
                                        }
                                    } else if (k == 6) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult61]; first_amp = new Float_t[first_mult61]; std::copy(first_tof61, first_tof61 + first_mult61, first_tof); std::copy(first_amp61, first_amp61 + first_mult61, first_amp); first_mult = first_mult61; cathode_detn_index = 24; break;
                                            case 2: first_tof = new Double_t[first_mult62]; first_amp = new Float_t[first_mult62]; std::copy(first_tof62, first_tof62 + first_mult62, first_tof); std::copy(first_amp62, first_amp62 + first_mult62, first_amp); first_mult = first_mult62; cathode_detn_index = 25; break;
                                            case 3: first_tof = new Double_t[first_mult63]; first_amp = new Float_t[first_mult63]; std::copy(first_tof63, first_tof63 + first_mult63, first_tof); std::copy(first_amp63, first_amp63 + first_mult63, first_amp); first_mult = first_mult63; cathode_detn_index = 26; break;
                                            case 4: first_tof = new Double_t[first_mult64]; first_amp = new Float_t[first_mult64]; std::copy(first_tof64, first_tof64 + first_mult64, first_tof); std::copy(first_amp64, first_amp64 + first_mult64, first_amp); first_mult = first_mult64; cathode_detn_index = 27; break;
                                        }
                                    } else if (k == 7) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult71]; first_amp = new Float_t[first_mult71]; std::copy(first_tof71, first_tof71 + first_mult71, first_tof); std::copy(first_amp71, first_amp71 + first_mult71, first_amp); first_mult = first_mult71; cathode_detn_index = 28; break;
                                            case 2: first_tof = new Double_t[first_mult72]; first_amp = new Float_t[first_mult72]; std::copy(first_tof72, first_tof72 + first_mult72, first_tof); std::copy(first_amp72, first_amp72 + first_mult72, first_amp); first_mult = first_mult72; cathode_detn_index = 29; break;
                                            case 3: first_tof = new Double_t[first_mult73]; first_amp = new Float_t[first_mult73]; std::copy(first_tof73, first_tof73 + first_mult73, first_tof); std::copy(first_amp73, first_amp73 + first_mult73, first_amp); first_mult = first_mult73; cathode_detn_index = 30; break;
                                            case 4: first_tof = new Double_t[first_mult74]; first_amp = new Float_t[first_mult74]; std::copy(first_tof74, first_tof74 + first_mult74, first_tof); std::copy(first_amp74, first_amp74 + first_mult74, first_amp); first_mult = first_mult74; cathode_detn_index = 31; break;
                                        }
                                    } else if (k == 8) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult81]; first_amp = new Float_t[first_mult81]; std::copy(first_tof81, first_tof81 + first_mult81, first_tof); std::copy(first_amp81, first_amp81 + first_mult81, first_amp); first_mult = first_mult81; cathode_detn_index = 32; break;
                                            case 2: first_tof = new Double_t[first_mult82]; first_amp = new Float_t[first_mult82]; std::copy(first_tof82, first_tof82 + first_mult82, first_tof); std::copy(first_amp82, first_amp82 + first_mult82, first_amp); first_mult = first_mult82; cathode_detn_index = 33; break;
                                            case 3: first_tof = new Double_t[first_mult83]; first_amp = new Float_t[first_mult83]; std::copy(first_tof83, first_tof83 + first_mult83, first_tof); std::copy(first_amp83, first_amp83 + first_mult83, first_amp); first_mult = first_mult83; cathode_detn_index = 34; break;
                                            case 4: first_tof = new Double_t[first_mult84]; first_amp = new Float_t[first_mult84]; std::copy(first_tof84, first_tof84 + first_mult84, first_tof); std::copy(first_amp84, first_amp84 + first_mult84, first_amp); first_mult = first_mult84; cathode_detn_index = 35; break;
                                        }
                                    } else if (k == 9) {
                                        switch (subdet) {
                                            case 1: first_tof = new Double_t[first_mult91]; first_amp = new Float_t[first_mult91]; std::copy(first_tof91, first_tof91 + first_mult91, first_tof); std::copy(first_amp91, first_amp91 + first_mult91, first_amp); first_mult = first_mult91; cathode_detn_index = 36; break;
                                            case 2: first_tof = new Double_t[first_mult92]; first_amp = new Float_t[first_mult92]; std::copy(first_tof92, first_tof92 + first_mult92, first_tof); std::copy(first_amp92, first_amp92 + first_mult92, first_amp); first_mult = first_mult92; cathode_detn_index = 37; break;
                                            case 3: first_tof = new Double_t[first_mult93]; first_amp = new Float_t[first_mult93]; std::copy(first_tof93, first_tof93 + first_mult93, first_tof); std::copy(first_amp93, first_amp93 + first_mult93, first_amp); first_mult = first_mult93; cathode_detn_index = 38; break;
                                            case 4: first_tof = new Double_t[first_mult94]; first_amp = new Float_t[first_mult94]; std::copy(first_tof94, first_tof94 + first_mult94, first_tof); std::copy(first_amp94, first_amp94 + first_mult94, first_amp); first_mult = first_mult94; cathode_detn_index = 39; break;
                                        }
                                    }
                                if (first_mult > 0) {
                                for (int m = 0; m < first_mult; m++) { // loop over cathode hits
                                    if (first_tof[m] >  anode_tof[l] && first_tof[m] < anode_tof[l]+time_for_coincidence && first_amp[m]>150 && cathode_detn_index>=4*k && cathode_detn_index<4*(k+1)) { 
                                        addCoincidences(cathode_detn_index);
                                        addTof(cathode_detn_index, first_tof[m]);
                                        addAmp(cathode_detn_index, first_amp[m]);
                                        }
                                    }
                                    delete [] first_tof; // free the memory allocated for the first tof
                                    delete [] first_amp; // free the memory allocated for the first amp
                                }
                            }
                            }
                            delete [] anode_tof; // free the memory allocated for the anode tof
                            delete [] anode_amp; // free the memory allocated for the anode amp
                        }
                    }
                }
            }
        }  
        outtree->Fill();
        // Reset the variables for the next entry
        mult01 = 0, mult02 = 0, mult03 = 0, mult04 = 0, mult11 = 0, mult12 = 0, mult13 = 0, mult14 = 0, mult21 = 0, mult22 = 0, mult23 = 0, mult24 = 0, mult31 = 0, mult32 = 0, mult33 = 0, mult34 = 0, mult41 = 0, mult42 = 0, mult43 = 0, mult44 = 0, mult51 = 0, mult52 = 0, mult53 = 0, mult54 = 0, mult61 = 0, mult62 = 0, mult63 = 0, mult64 = 0, mult71 = 0, mult72 = 0, mult73 = 0, mult74 = 0, mult81 = 0, mult82 = 0, mult83 = 0, mult84 = 0, mult91 = 0, mult92 = 0, mult93 = 0, mult94 = 0;
        mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
       std::fill(std::begin(detn_cathode), std::end(detn_cathode), 0);
    }
        // Write the output tree to the file
outtree->Write();
outfile->Close();

}
