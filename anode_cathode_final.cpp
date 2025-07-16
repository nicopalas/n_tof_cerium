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
Int_t time_for_coincidence = 10; //time difference between signals CHECK IF TIME IS IN ns!!!!
Int_t max_mul = 5; //max number of signals in coincidence in a single detector
Int_t mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
Double_t tof0[5],tof1[5],tof2[5],tof3[5],tof4[5],tof5[5],tof6[5],tof7[5],tof8[5],tof9[5]; //must be hardcoded if arrays are defined at file scope
Float_t amp0[5],amp1[5],amp2[5],amp3[5],amp4[5],amp5[5],amp6[5],amp7[5],amp8[5],amp9[5];  //must be hardcoded if arrays are defined at file scope
Double_t gamma_flash0[5], gamma_flash1[5], gamma_flash2[5], gamma_flash3[5], gamma_flash4[5], gamma_flash5[5], gamma_flash6[5], gamma_flash7[5], gamma_flash8[5], gamma_flash9[5];
Int_t detn_all[10]; 
Int_t cathode_detn[40];
Double_t tof01[20], tof02[20], tof03[20], tof04[20], tof11[20], tof12[20], tof13[20], tof14[20], tof21[20], tof22[20], tof23[20], tof24[20], tof31[20], tof32[20], tof33[20], tof34[20], tof41[20], tof42[20], tof43[20], tof44[20], tof51[20], tof52[20], tof53[20], tof54[20], tof61[20], tof62[20], tof63[20], tof64[20], tof71[20], tof72[20], tof73[20], tof74[20], tof81[20], tof82[20], tof83[20], tof84[20], tof91[20], tof92[20], tof93[20], tof94[20];
Float_t amp01[20], amp02[20], amp03[20], amp04[20], amp11[20], amp12[20], amp13[20], amp14[20], amp21[20], amp22[20], amp23[20], amp24[20], amp31[20], amp32[20], amp33[20], amp34[20], amp41[20], amp42[20], amp43[20], amp44[20], amp51[20], amp52[20], amp53[20], amp54[20], amp61[20], amp62[20], amp63[20], amp64[20], amp71[20], amp72[20], amp73[20], amp74[20], amp81[20], amp82[20], amp83[20], amp84[20], amp91[20], amp92[20], amp93[20], amp94[20];
Int_t mult01 = 0, mult02 = 0, mult03 = 0, mult04 = 0, mult11 = 0, mult12 = 0, mult13 = 0, mult14 = 0, mult21 = 0, mult22 = 0, mult23 = 0, mult24 = 0, mult31 = 0, mult32 = 0, mult33 = 0, mult34 = 0, mult41 = 0, mult42 = 0, mult43 = 0, mult44 = 0, mult51 = 0, mult52 = 0, mult53 = 0, mult54 = 0, mult61 = 0, mult62 = 0, mult63 = 0, mult64 = 0, mult71 = 0, mult72 = 0, mult73 = 0, mult74 = 0, mult81 = 0, mult82 = 0, mult83 = 0, mult84 = 0, mult91 = 0, mult92 = 0, mult93 = 0, mult94 = 0;
Int_t cathode_RunNumber, cathode_time, cathode_BunchNumber, cathode_PSpulse;
Double_t cathode_psTime;
Float_t cathode_PulseIntensity;
Int_t first_detn_all[40]; 
Int_t detn_cathode[40];
Double_t neutron_energy;
Bool_t addCoincidences(Int_t the_detn) {
    if(the_detn<0 || the_detn>94 || the_detn%10>4) {        
        std::cout << "ERROR in addCoincidences: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(the_detn==01) { mult01++; detn_cathode[0]++; }
    if(the_detn==02) { mult02++; detn_cathode[1]++; }
    if(the_detn==03) { mult03++; detn_cathode[2]++; }
    if(the_detn==04) { mult04++; detn_cathode[3]++; }
    if(the_detn==11)  { mult11++; detn_cathode[4]++; }
    if(the_detn==12) { mult12++; detn_cathode[5]++; }
    if(the_detn==13) { mult13++; detn_cathode[6]++; }
    if(the_detn==14) { mult14++; detn_cathode[7]++; }
    if(the_detn==21) { mult21++; detn_cathode[8]++; }
    if(the_detn==22) { mult22++; detn_cathode[9]++; }
    if(the_detn==23) { mult23++; detn_cathode[10]++; }
    if(the_detn==24) { mult24++; detn_cathode[11]++; }
    if(the_detn==31) { mult31++; detn_cathode[12]++; }
    if(the_detn==32) { mult32++; detn_cathode[13]++; }
    if(the_detn==33) { mult33++; detn_cathode[14]++; }
    if(the_detn==34) { mult34++; detn_cathode[15]++; }
    if(the_detn==41) { mult41++; detn_cathode[16]++; }
    if(the_detn==42) { mult42++; detn_cathode[17]++; }
    if(the_detn==43) { mult43++; detn_cathode[18]++; }
    if(the_detn==44) { mult44++; detn_cathode[19]++; }
    if(the_detn==51) { mult51++; detn_cathode[20]++; }
    if(the_detn==52) { mult52++; detn_cathode[21]++; }
    if(the_detn==53) { mult53++; detn_cathode[22]++; }
    if(the_detn==54) { mult54++; detn_cathode[23]++; }
    if(the_detn==61) { mult61++; detn_cathode[24]++; }
    if(the_detn==62) { mult62++; detn_cathode[25]++; }
    if(the_detn==63) { mult63++; detn_cathode[26]++; }
    if(the_detn==64) { mult64++; detn_cathode[27]++; }
    if(the_detn==71) { mult71++; detn_cathode[28]++; }
    if(the_detn==72) { mult72++; detn_cathode[29]++; }
    if(the_detn==73) { mult73++; detn_cathode[30]++; }
    if(the_detn==74) { mult74++; detn_cathode[31]++; }
    if(the_detn==81) { mult81++; detn_cathode[32]++; }
    if(the_detn==82) { mult82++; detn_cathode[33]++; }
    if(the_detn==83) { mult83++; detn_cathode[34]++; }
    if(the_detn==84) { mult84++; detn_cathode[35]++; }
    if(the_detn==91) { mult91++; detn_cathode[36]++; }
    if(the_detn==92) { mult92++; detn_cathode[37]++; }
    if(the_detn==93) { mult93++; detn_cathode[38]++; }
    if(the_detn==94) { mult94++; detn_cathode[39]++; }
    return true;
}

Bool_t addTof(Int_t detn_index, Double_t the_tof) {
    if(detn_index<1 ||detn_index > 94 || detn_index % 10 > 4) {         
        std::cout << "ERROR in addTof: detn" << detn_index << "value: " <<detn_index << " not valid" << std::endl;
        return false;
    }
    switch (detn_index) {
        case 01: if (mult01 > 0 && mult01 < 10) tof01[mult01 - 1] = the_tof; else if (mult01 >= 10) { std::cout << "ERROR in addTof: mult01 exceeds 10" << std::endl; return false; } break;
        case 02: if (mult02 > 0 && mult02 < 10) tof02[mult02 - 1] = the_tof; else if (mult02 >= 10) { std::cout << "ERROR in addTof: mult02 exceeds 10" << std::endl; return false; } break;
        case 03: if (mult03 > 0 && mult03 < 10) tof03[mult03 - 1] = the_tof; else if (mult03 >= 10) { std::cout << "ERROR in addTof: mult03 exceeds 10" << std::endl; return false; } break;
        case 04: if (mult04 > 0 && mult04 < 10) tof04[mult04 - 1] = the_tof; else if (mult04 >= 10) { std::cout << "ERROR in addTof: mult04 exceeds 10" << std::endl; return false; } break;
        case 11: if (mult11 > 0 && mult11 < 10) tof11[mult11 - 1] = the_tof; else if (mult11 >= 10) { std::cout << "ERROR in addTof: mult11 exceeds 10" << std::endl; return false; } break;
        case 12: if (mult12 > 0 && mult12 < 10) tof12[mult12 - 1] = the_tof; else if (mult12 >= 10) { std::cout << "ERROR in addTof: mult12 exceeds 10" << std::endl; return false; } break;
        case 13: if (mult13 > 0 && mult13 < 10) tof13[mult13 - 1] = the_tof; else if (mult13 >= 10) { std::cout << "ERROR in addTof: mult13 exceeds 10" << std::endl; return false; } break;
        case 14: if (mult14 > 0 && mult14 < 10) tof14[mult14 - 1] = the_tof; else if (mult14 >= 10) { std::cout << "ERROR in addTof: mult14 exceeds 10" << std::endl; return false; } break;
        case 21: if (mult21 > 0 && mult21 < 10) tof21[mult21 - 1] = the_tof; else if (mult21 >= 10) { std::cout << "ERROR in addTof: mult21 exceeds 10" << std::endl; return false; } break;
        case 22: if (mult22 > 0 && mult22 < 10) tof22[mult22 - 1] = the_tof; else if (mult22 >= 10) { std::cout << "ERROR in addTof: mult22 exceeds 10" << std::endl; return false; } break;
        case 23: if (mult23 > 0 && mult23 < 10) tof23[mult23 - 1] = the_tof; else if (mult23 >= 10) { std::cout << "ERROR in addTof: mult23 exceeds 10" << std::endl; return false; } break;
        case 24: if (mult24 > 0 && mult24 < 10) tof24[mult24 - 1] = the_tof; else if (mult24 >= 10) { std::cout << "ERROR in addTof: mult24 exceeds 10" << std::endl; return false; } break;
        case 31: if (mult31 > 0 && mult31 < 10) tof31[mult31 - 1] = the_tof; else if (mult31 >= 10) { std::cout << "ERROR in addTof: mult31 exceeds 10" << std::endl; return false; } break;
        case 32: if (mult32 > 0 && mult32 < 10) tof32[mult32 - 1] = the_tof; else if (mult32 >= 10) { std::cout << "ERROR in addTof: mult32 exceeds 10" << std::endl; return false; } break;
        case 33: if (mult33 > 0 && mult33 < 10) tof33[mult33 - 1] = the_tof; else if (mult33 >= 10) { std::cout << "ERROR in addTof: mult33 exceeds 10" << std::endl; return false; } break;
        case 34: if (mult34 > 0 && mult34 < 10) tof34[mult34 - 1] = the_tof; else if (mult34 >= 10) { std::cout << "ERROR in addTof: mult34 exceeds 10" << std::endl; return false; } break;
        case 41: if (mult41 > 0 && mult41 < 10) tof41[mult41 - 1] = the_tof; else if (mult41 >= 10) { std::cout << "ERROR in addTof: mult41 exceeds 10" << std::endl; return false; } break;
        case 42: if (mult42 > 0 && mult42 < 10) tof42[mult42 - 1] = the_tof; else if (mult42 >= 10) { std::cout << "ERROR in addTof: mult42 exceeds 10" << std::endl; return false; } break;
        case 43: if (mult43 > 0 && mult43 < 10) tof43[mult43 - 1] = the_tof; else if (mult43 >= 10) { std::cout << "ERROR in addTof: mult43 exceeds 10" << std::endl; return false; } break;
        case 44: if (mult44 > 0 && mult44 < 10) tof44[mult44 - 1] = the_tof; else if (mult44 >= 10) { std::cout << "ERROR in addTof: mult44 exceeds 10" << std::endl; return false; } break;
        case 51: if (mult51 > 0 && mult51 < 10) tof51[mult51 - 1] = the_tof; else if (mult51 >= 10) { std::cout << "ERROR in addTof: mult51 exceeds 10" << std::endl; return false; } break;
        case 52: if (mult52 > 0 && mult52 < 10) tof52[mult52 - 1] = the_tof; else if (mult52 >= 10) { std::cout << "ERROR in addTof: mult52 exceeds 10" << std::endl; return false; } break;
        case 53: if (mult53 > 0 && mult53 < 10) tof53[mult53 - 1] = the_tof; else if (mult53 >= 10) { std::cout << "ERROR in addTof: mult53 exceeds 10" << std::endl; return false; } break;
        case 54: if (mult54 > 0 && mult54 < 10) tof54[mult54 - 1] = the_tof; else if (mult54 >= 10) { std::cout << "ERROR in addTof: mult54 exceeds 10" << std::endl; return false; } break;
        case 61: if (mult61 > 0 && mult61 < 10) tof61[mult61 - 1] = the_tof; else if (mult61 >= 10) { std::cout << "ERROR in addTof: mult61 exceeds 10" << std::endl; return false; } break;
        case 62: if (mult62 > 0 && mult62 < 10) tof62[mult62 - 1] = the_tof; else if (mult62 >= 10) { std::cout << "ERROR in addTof: mult62 exceeds 10" << std::endl; return false; } break;
        case 63: if (mult63 > 0 && mult63 < 10) tof63[mult63 - 1] = the_tof; else if (mult63 >= 10) { std::cout << "ERROR in addTof: mult63 exceeds 10" << std::endl; return false; } break;
        case 64: if (mult64 > 0 && mult64 < 10) tof64[mult64 - 1] = the_tof; else if (mult64 >= 10) { std::cout << "ERROR in addTof: mult64 exceeds 10" << std::endl; return false; } break;
        case 71: if (mult71 > 0 && mult71 < 10) tof71[mult71 - 1] = the_tof; else if (mult71 >= 10) { std::cout << "ERROR in addTof: mult71 exceeds 10" << std::endl; return false; } break;
        case 72: if (mult72 > 0 && mult72 < 10) tof72[mult72 - 1] = the_tof; else if (mult72 >= 10) { std::cout << "ERROR in addTof: mult72 exceeds 10" << std::endl; return false; } break;
        case 73: if (mult73 > 0 && mult73 < 10) tof73[mult73 - 1] = the_tof; else if (mult73 >= 10) { std::cout << "ERROR in addTof: mult73 exceeds 10" << std::endl; return false; } break;
        case 74: if (mult74 > 0 && mult74 < 10) tof74[mult74 - 1] = the_tof; else if (mult74 >= 10) { std::cout << "ERROR in addTof: mult74 exceeds 10" << std::endl; return false; } break;
        case 81: if (mult81 > 0 && mult81 < 10) tof81[mult81 - 1] = the_tof; else if (mult81 >= 10) { std::cout << "ERROR in addTof: mult81 exceeds 10" << std::endl; return false; } break;
        case 82: if (mult82 > 0 && mult82 < 10) tof82[mult82 - 1] = the_tof; else if (mult82 >= 10) { std::cout << "ERROR in addTof: mult82 exceeds 10" << std::endl; return false; } break;
        case 83: if (mult83 > 0 && mult83 < 10) tof83[mult83 - 1] = the_tof; else if (mult83 >= 10) { std::cout << "ERROR in addTof: mult83 exceeds 10" << std::endl; return false; } break;
        case 84: if (mult84 > 0 && mult84 < 10) tof84[mult84 - 1] = the_tof; else if (mult84 >= 10) { std::cout << "ERROR in addTof: mult84 exceeds 10" << std::endl; return false; } break;
        case 91: if (mult91 > 0 && mult91 < 10) tof91[mult91 - 1] = the_tof; else if (mult91 >= 10) { std::cout << "ERROR in addTof: mult91 exceeds 10" << std::endl; return false; } break;
        case 92: if (mult92 > 0 && mult92 < 10) tof92[mult92 - 1] = the_tof; else if (mult92 >= 10) { std::cout << "ERROR in addTof: mult92 exceeds 10" << std::endl; return false; } break;
        case 93: if (mult93 > 0 && mult93 < 10) tof93[mult93 - 1] = the_tof; else if (mult93 >= 10) { std::cout << "ERROR in addTof: mult93 exceeds 10" << std::endl; return false; } break;
        case 94: if (mult94 > 0 && mult94 < 10) tof94[mult94 - 1] = the_tof; else if (mult94 >= 10) { std::cout << "ERROR in addTof: mult94 exceeds 10" << std::endl; return false; } break;
        default:
            std::cout << "ERROR in addTof: Unexpected the_detn value: " << detn_index << std::endl;
            return false;
    }
    return true;
}
Bool_t addAmp(Int_t detn_index, Float_t the_amp) {
    if(detn_index<1 ||detn_index > 94 || detn_index % 10 > 4) {         
        std::cout << "ERROR in addAmp: detn" << detn_index << "value: " <<detn_index << " not valid" << std::endl;
        return false;
    }
    switch (detn_index) {
        case 01: if (mult01 > 0 && mult01 < 10) amp01[mult01 - 1] = the_amp; else if (mult01 >= 10) { std::cout << "ERROR in addAmp: mult01 exceeds 10" << std::endl; return false; } break;
        case 02: if (mult02 > 0 && mult02 < 10) amp02[mult02 - 1] = the_amp; else if (mult02 >= 10) { std::cout << "ERROR in addAmp: mult02 exceeds 10" << std::endl; return false; } break;
        case 03: if (mult03 > 0 && mult03 < 10) amp03[mult03 - 1] = the_amp; else if (mult03 >= 10) { std::cout << "ERROR in addAmp: mult03 exceeds 10" << std::endl; return false; } break;
        case 04: if (mult04 > 0 && mult04 < 10) amp04[mult04 - 1] = the_amp; else if (mult04 >= 10) { std::cout << "ERROR in addAmp: mult04 exceeds 10" << std::endl; return false; } break;
        case 11: if (mult11 > 0 && mult11 < 10) amp11[mult11 - 1] = the_amp; else if (mult11 >= 10) { std::cout << "ERROR in addAmp: mult11 exceeds 10" << std::endl; return false; } break;
        case 12: if (mult12 > 0 && mult12 < 10) amp12[mult12 - 1] = the_amp; else if (mult12 >= 10) { std::cout << "ERROR in addAmp: mult12 exceeds 10" << std::endl; return false; } break;
        case 13: if (mult13 > 0 && mult13 < 10) amp13[mult13 - 1] = the_amp; else if (mult13 >= 10) { std::cout << "ERROR in addAmp: mult13 exceeds 10" << std::endl; return false; } break;
        case 14: if (mult14 > 0 && mult14 < 10) amp14[mult14 - 1] = the_amp; else if (mult14 >= 10) { std::cout << "ERROR in addAmp: mult14 exceeds 10" << std::endl; return false; } break;
        case 21: if (mult21 > 0 && mult21 < 10) amp21[mult21 - 1] = the_amp; else if (mult21 >= 10) { std::cout << "ERROR in addAmp: mult21 exceeds 10" << std::endl; return false; } break;
        case 22: if (mult22 > 0 && mult22 < 10) amp22[mult22 - 1] = the_amp; else if (mult22 >= 10) { std::cout << "ERROR in addAmp: mult22 exceeds 10" << std::endl; return false; } break;
        case 23: if (mult23 > 0 && mult23 < 10) amp23[mult23 - 1] = the_amp; else if (mult23 >= 10) { std::cout << "ERROR in addAmp: mult23 exceeds 10" << std::endl; return false; } break;
        case 24: if (mult24 > 0 && mult24 < 10) amp24[mult24 - 1] = the_amp; else if (mult24 >= 10) { std::cout << "ERROR in addAmp: mult24 exceeds 10" << std::endl; return false; } break;
        case 31: if (mult31 > 0 && mult31 < 10) amp31[mult31 - 1] = the_amp; else if (mult31 >= 10) { std::cout << "ERROR in addAmp: mult31 exceeds 10" << std::endl; return false; } break;
        case 32: if (mult32 > 0 && mult32 < 10) amp32[mult32 - 1] = the_amp; else if (mult32 >= 10) { std::cout << "ERROR in addAmp: mult32 exceeds 10" << std::endl; return false; } break;
        case 33: if (mult33 > 0 && mult33 < 10) amp33[mult33 - 1] = the_amp; else if (mult33 >= 10) { std::cout << "ERROR in addAmp: mult33 exceeds 10" << std::endl; return false; } break;
        case 34: if (mult34 > 0 && mult34 < 10) amp34[mult34 - 1] = the_amp; else if (mult34 >= 10) { std::cout << "ERROR in addAmp: mult34 exceeds 10" << std::endl; return false; } break;
        case 41: if (mult41 > 0 && mult41 < 10) amp41[mult41 - 1] = the_amp; else if (mult41 >= 10) { std::cout << "ERROR in addAmp: mult41 exceeds 10" << std::endl; return false; } break;
        case 42: if (mult42 > 0 && mult42 < 10) amp42[mult42 - 1] = the_amp; else if (mult42 >= 10) { std::cout << "ERROR in addAmp: mult42 exceeds 10" << std::endl; return false; } break;
        case 43: if (mult43 > 0 && mult43 < 10) amp43[mult43 - 1] = the_amp; else if (mult43 >= 10) { std::cout << "ERROR in addAmp: mult43 exceeds 10" << std::endl; return false; } break;
        case 44: if (mult44 > 0 && mult44 < 10) amp44[mult44 - 1] = the_amp; else if (mult44 >= 10) { std::cout << "ERROR in addAmp: mult44 exceeds 10" << std::endl; return false; } break;
        case 51: if (mult51 > 0 && mult51 < 10) amp51[mult51 - 1] = the_amp; else if (mult51 >= 10) { std::cout << "ERROR in addAmp: mult51 exceeds 10" << std::endl; return false; } break;
        case 52: if (mult52 > 0 && mult52 < 10) amp52[mult52 - 1] = the_amp; else if (mult52 >= 10) { std::cout << "ERROR in addAmp: mult52 exceeds 10" << std::endl; return false; } break;
        case 53: if (mult53 > 0 && mult53 < 10) amp53[mult53 - 1] = the_amp; else if (mult53 >= 10) { std::cout << "ERROR in addAmp: mult53 exceeds 10" << std::endl; return false; } break;
        case 54: if (mult54 > 0 && mult54 < 10) amp54[mult54 - 1] = the_amp; else if (mult54 >= 10) { std::cout << "ERROR in addAmp: mult54 exceeds 10" << std::endl; return false; } break;
        case 61: if (mult61 > 0 && mult61 < 10) amp61[mult61 - 1] = the_amp; else if (mult61 >= 10) { std::cout << "ERROR in addAmp: mult61 exceeds 10" << std::endl; return false; } break;
        case 62: if (mult62 > 0 && mult62 < 10) amp62[mult62 - 1] = the_amp; else if (mult62 >= 10) { std::cout << "ERROR in addAmp: mult62 exceeds 10" << std::endl; return false; } break;
        case 63: if (mult63 > 0 && mult63 < 10) amp63[mult63 - 1] = the_amp; else if (mult63 >= 10) { std::cout << "ERROR in addAmp: mult63 exceeds 10" << std::endl; return false; } break;
        case 64: if (mult64 > 0 && mult64 < 10) amp64[mult64 - 1] = the_amp; else if (mult64 >= 10) { std::cout << "ERROR in addAmp: mult64 exceeds 10" << std::endl; return false; } break;
        case 71: if (mult71 > 0 && mult71 < 10) amp71[mult71 - 1] = the_amp; else if (mult71 >= 10) { std::cout << "ERROR in addAmp: mult71 exceeds 10" << std::endl; return false; } break;
        case 72: if (mult72 > 0 && mult72 < 10) amp72[mult72 - 1] = the_amp; else if (mult72 >= 10) { std::cout << "ERROR in addAmp: mult72 exceeds 10" << std::endl; return false; } break;
        case 73: if (mult73 > 0 && mult73 < 10) amp73[mult73 - 1] = the_amp; else if (mult73 >= 10) { std::cout << "ERROR in addAmp: mult73 exceeds 10" << std::endl; return false; } break;
        case 74: if (mult74 > 0 && mult74 < 10) amp74[mult74 - 1] = the_amp; else if (mult74 >= 10) { std::cout << "ERROR in addAmp: mult74 exceeds 10" << std::endl; return false; } break;
        case 81: if (mult81 > 0 && mult81 < 10) amp81[mult81 - 1] = the_amp; else if (mult81 >= 10) { std::cout << "ERROR in addAmp: mult81 exceeds 10" << std::endl; return false; } break;
        case 82: if (mult82 > 0 && mult82 < 10) amp82[mult82 - 1] = the_amp; else if (mult82 >= 10) { std::cout << "ERROR in addAmp: mult82 exceeds 10" << std::endl; return false; } break;
        case 83: if (mult83 > 0 && mult83 < 10) amp83[mult83 - 1] = the_amp; else if (mult83 >= 10) { std::cout << "ERROR in addAmp: mult83 exceeds 10" << std::endl; return false; } break;
        case 84: if (mult84 > 0 && mult84 < 10) amp84[mult84 - 1] = the_amp; else if (mult84 >= 10) { std::cout << "ERROR in addAmp: mult84 exceeds 10" << std::endl; return false; } break;
        case 91: if (mult91 > 0 && mult91 < 10) amp91[mult91 - 1] = the_amp; else if (mult91 >= 10) { std::cout << "ERROR in addAmp: mult91 exceeds 10" << std::endl; return false; } break;
        case 92: if (mult92 > 0 && mult92 < 10) amp92[mult92 - 1] = the_amp; else if (mult92 >= 10) { std::cout << "ERROR in addAmp: mult92 exceeds 10" << std::endl; return false; } break;
        case 93: if (mult93 > 0 && mult93 < 10) amp93[mult93 - 1] = the_amp; else if (mult93 >= 10) { std::cout << "ERROR in addAmp: mult93 exceeds 10" << std::endl; return false; } break;
        case 94: if (mult94 > 0 && mult94 < 10) amp94[mult94 - 1] = the_amp; else if (mult94 >= 10) { std::cout << "ERROR in addAmp: mult94 exceeds 10" << std::endl; return false; } break;
        default:
            std::cout << "ERROR in addAmp: Unexpected detn_index value: " << detn_index << std::endl;
            return false;
    }
    return true;
}

void anode_cathode_final(Int_t run_number) {
    std::string filename = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Anodes/anodes_final/out_run" + std::to_string(run_number) + ".root";
    TFile *file = TFile::Open(filename.c_str()) ;//coincidences file
    TTree *intree = (TTree*)file->Get("coincidences");
    std::string filename2 = "/nucl_lustre/n_tof_INTC_P_665/DATA/run" + std::to_string(run_number) + ".root";
    TFile *file2 = TFile::Open(filename2.c_str());// data from the pickup and the cathode
    TTree *cathodeTree = (TTree*)file2->Get("PPAC");
    TTree *pickupTree = (TTree*)file2->Get("PKUP");


    Int_t cathode_RunNumber, cathode_time, cathode_BunchNumber, cathode_PSpulse;
    Double_t cathode_psTime;
    Float_t cathode_PulseIntensity;
    Float_t cathode_amp;
    Double_t cathode_tof;
    Int_t cathode_detn;

    Int_t pickup_RunNumber, pickup_time, pickup_BunchNumber, pickup_PSpulse;
    Double_t pickup_psTime;
    Float_t pickup_PulseIntensity;
    Float_t pickup_amp;
    Double_t pickup_tof;


    Long64_t nentries = intree->GetEntries(); // getting the number of entries to loop over them
    Long64_t ncathodeEntries = cathodeTree->GetEntries();
    Long64_t npickupEntries = pickupTree->GetEntries();
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
    intree->SetBranchAddress("gamma_flash0", gamma_flash0);
    intree->SetBranchAddress("gamma_flash1", gamma_flash1);
    intree->SetBranchAddress("gamma_flash2", gamma_flash2);
    intree->SetBranchAddress("gamma_flash3", gamma_flash3);
    intree->SetBranchAddress("gamma_flash4", gamma_flash4);
    intree->SetBranchAddress("gamma_flash5", gamma_flash5);
    intree->SetBranchAddress("gamma_flash6", gamma_flash6);
    intree->SetBranchAddress("gamma_flash7", gamma_flash7);
    intree->SetBranchAddress("gamma_flash8", gamma_flash8);
    intree->SetBranchAddress("gamma_flash9", gamma_flash9);


    cathodeTree->SetBranchStatus("*", 0); // turn off all branches
    cathodeTree->SetBranchStatus("RunNumber", 1); // Int_t
    cathodeTree->SetBranchStatus("time", 1); // Int_t
    cathodeTree->SetBranchStatus("psTime", 1); // Double_t
    cathodeTree->SetBranchStatus("BunchNumber", 1); // Int_t
    cathodeTree->SetBranchStatus("PSpulse", 1); // Int_t
    cathodeTree->SetBranchStatus("PulseIntensity", 1); // Float_t
    cathodeTree->SetBranchStatus("detn", 1); // Int_t
    cathodeTree->SetBranchStatus("tof", 1); // Double_t
    cathodeTree->SetBranchStatus("amp", 1); // Float_t

    cathodeTree->SetBranchAddress("RunNumber", &cathode_RunNumber);
    cathodeTree->SetBranchAddress("time", &cathode_time);
    cathodeTree->SetBranchAddress("psTime", &cathode_psTime);
    cathodeTree->SetBranchAddress("BunchNumber", &cathode_BunchNumber);
    cathodeTree->SetBranchAddress("PSpulse", &cathode_PSpulse);
    cathodeTree->SetBranchAddress("PulseIntensity", &cathode_PulseIntensity);
    cathodeTree->SetBranchAddress("detn", &cathode_detn);
    cathodeTree->SetBranchAddress("tof", &cathode_tof);
    cathodeTree->SetBranchAddress("amp", &cathode_amp);

    pickupTree->SetBranchStatus("*", 0); // turn off all branches
    pickupTree->SetBranchStatus("RunNumber", 1); // Int_t
    pickupTree->SetBranchStatus("time", 1); // Int_t
    pickupTree->SetBranchStatus("psTime", 1); // Double_t
    pickupTree->SetBranchStatus("BunchNumber", 1); // Int_t
    pickupTree->SetBranchStatus("PSpulse", 1); // Int_t
    pickupTree->SetBranchStatus("PulseIntensity", 1); // Float_t
    pickupTree->SetBranchStatus("tof", 1); // Double_t
    pickupTree->SetBranchAddress("RunNumber", &pickup_RunNumber);
    pickupTree->SetBranchAddress("time", &pickup_time);
    pickupTree->SetBranchAddress("psTime", &pickup_psTime);
    pickupTree->SetBranchAddress("BunchNumber", &pickup_BunchNumber);
    pickupTree->SetBranchAddress("PSpulse", &pickup_PSpulse);
    pickupTree->SetBranchAddress("PulseIntensity", &pickup_PulseIntensity);
    pickupTree->SetBranchAddress("tof", &pickup_tof);

    // Create output file and tree
    std::string out_filename = "nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Cathodes/cathode_coincidences_final" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(out_filename.c_str(), "RECREATE");
    TTree *outtree = new TTree("cathode_coincidences", "cathode_coincidences");
    outtree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    outtree->Branch("BunchNumber", &BunchNumber, "BunchNumber/I");
    outtree->Branch("PSpulse", &PSpulse, "PSpulse/I");
    outtree->Branch("time", &eventTime, "time/I");
    outtree->Branch("psTime", &psTime, "psTime/D");
    outtree->Branch("PulseIntensity", &PulseIntensity, "PulseIntensity/F");
    outtree->Branch("cathode_detn", cathode_detn, "cathode_detn[40]/I");
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
    outtree->Branch("detn", detn_all, "detn[20]/I");
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

    std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Float_t, Double_t>> cathode_signals;
    std::vector<std::vector<std::tuple<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>, Double_t, Int_t>>> configurations;


    for (Long64_t i = 0; i < ncathodeEntries; ++i) {
        cathodeTree->GetEntry(i);
        cathode_signals.emplace_back(cathode_RunNumber, cathode_time, cathode_psTime, cathode_BunchNumber, cathode_PSpulse, cathode_PulseIntensity, cathode_detn, cathode_tof, cathode_amp); //insert tuple at the end

    }
    std::cout <<"INFO: number of signals: "<< cathode_signals.size() << std::endl; // print number of signals
    std::vector<bool> used(cathode_signals.size(), false); // track used signals to avoid double counting
    for (size_t j = 0; j < npickupEntries; ++j) {
        pickupTree->GetEntry(j);
        for(size_t i = 0; i < cathode_signals.size(); ++i){
            if (pickup_RunNumber == std::get<0>(cathode_signals[i]) && 
                pickup_time == std::get<1>(cathode_signals[i]) && 
                pickup_psTime == std::get<2>(cathode_signals[i]) && 
                pickup_BunchNumber == std::get<3>(cathode_signals[i]) && 
                pickup_PSpulse == std::get<4>(cathode_signals[i]) &&
                pickup_PulseIntensity== std::get<5>(cathode_signals[i])) {
                    
                    std::get<7>(cathode_signals[i])=std::get<7>(cathode_signals[i])-pickup_tof;
            }
        }
    }   
    std::cout << "INFO: signals calibrated with the pickup!" << std::endl; // print after calibration

    for (Long64_t i = 0; i < 1000; i++) {
        intree->GetEntry(i);
        std::cout << "Event: " << i << std::endl;
        std::vector<std::tuple<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>, Double_t, Int_t>> coincidences_cathode;
        for (Long64_t j = 0; j < cathode_signals.size(); j++) {
            if (!used[j]){ // skip already used signals
            if (std::get<0>(cathode_signals[j]) == RunNumber && std::get<1>(cathode_signals[j]) == eventTime && std::get<2>(cathode_signals[j]) == psTime && std::get<3>(cathode_signals[j]) == BunchNumber && std::get<4>(cathode_signals[j]) == PSpulse && std::get<5>(cathode_signals[j]) == PulseIntensity) {
                    
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
                                for (int l = 0; l < anode_mult; l++) { // loop over cathode hits
                                    if (std::get<7>(cathode_signals[j]) >  anode_tof[l]-50 && std::get<7>(cathode_signals[j]) < anode_tof[l]+400 && std::get<8>(cathode_signals[j])<=2*anode_amp[l] && std::get<8>(cathode_signals[j])>300 && std::get<6>(cathode_signals[j])>=10*k && std::get<6>(cathode_signals[j])<10*k+4) { 
                                        used[j] = true; // mark the signal as used)
                                        coincidences_cathode.push_back(std::make_tuple(cathode_signals[j], anode_tof[l], k)); // add the signal to the coincidences
                                        }
                                    }
                                
                            delete [] anode_tof; // free the memory allocated for the anode tof
                            delete [] anode_amp; // free the memory allocated for the anode amp
                        }
                    }
                }
            }
        }
        std::cout << "INFO: number of coincidences in the configuration: " << coincidences_cathode.size() << std::endl; // print number of coincidences
        if(coincidences_cathode.size() == 0) continue; // skip empty coincidences
        configurations.push_back(coincidences_cathode); // add the coincidences to the configurations
    }
    std::cout << "INFO: number of coincidences in the first step: " << configurations.size() << std::endl; // print number of coincidences
    Int_t lonelyCoincidence = 0;
    for (size_t i = 0; i < cathode_signals.size(); ++i) {
        if (used[i] == false && std::get<8>(cathode_signals[i]) > 300) {
            for (size_t j=0; j<configurations.size();j++){
                for (size_t k = 0; k < configurations[j].size(); ++k) {
                    if (std::get<0>(cathode_signals[i]) == std::get<0>(std::get<0>(configurations[j][k])) && std::get<8>(std::get<0>(configurations[j][k])) > 300 && std::get<5>(cathode_signals[i]) == std::get<5>(std::get<0>(configurations[j][k])) && std::get<4>(cathode_signals[i]) == std::get<4>(std::get<0>(configurations[j][k]))) {
                        if (std::get<1>(cathode_signals[i]) == std::get<1>(std::get<0>(configurations[j][k])) && std::get<3>(cathode_signals[i]) == std::get<3>(std::get<0>(configurations[j][k])) && std::get<4>(cathode_signals[i]) == std::get<4>(std::get<0>(configurations[j][k])) && std::get<2>(cathode_signals[i]) == std::get<2>(std::get<0>(configurations[j][k]))) {
                            if ( std::get<6>(cathode_signals[i])<=std::get<2>(configurations[j][k])*10+4 && std::get<6>(cathode_signals[i])>=std::get<2>(configurations[j][k])*10 && abs(std::get<7>(cathode_signals[i]) - std::get<7>(std::get<0>(configurations[j][k]))) < 20) {
                                lonelyCoincidence++;
                                configurations[j].push_back(std::make_tuple(cathode_signals[i], std::get<1>(configurations[j][k]), std::get<2>(configurations[j][k]))); // Add the signal with a placeholder value
                                used[i] = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (used[i]==true && std::get<8>(cathode_signals[i])>300) {
            Int_t coincidence_index = -1;
            for (Int_t j = 0; j < configurations.size(); ++j) {
                for (Int_t k = 0; k < configurations[j].size(); ++k) {
                    if (std::get<0>(cathode_signals[i]) == std::get<0>(std::get<0>(configurations[j][k])) && std::get<1>(std::get<0>(configurations[j][k])) == std::get<1>(cathode_signals[i]) && std::get<2>(std::get<0>(configurations[j][k])) == std::get<2>(cathode_signals[i])) {
                        if (std::get<3>(cathode_signals[i]) == std::get<3>(std::get<0>(configurations[j][k])) && std::get<4>(cathode_signals[i]) == std::get<4>(std::get<0>(configurations[j][k])) && std::get<5>(cathode_signals[i]) == std::get<5>(std::get<0>(configurations[j][k]))) {
                            if (std::get<6>(cathode_signals[i]) == std::get<6>(std::get<0>(configurations[j][k])) && std::get<7>(cathode_signals[i]) == std::get<7>(std::get<0>(configurations[j][k])) && std::get<8>(cathode_signals[i]) == std::get<8>(std::get<0>(configurations[j][k]))) {
                                coincidence_index = j;
                                break;
                            }
                        }
                    }
                }
                if (coincidence_index != -1) { //if you find the index stop the loop
                    break;
                }
            }
            for (Int_t k=0; k<configurations.size();k++){ //loop to find coincidences
                if (k!= coincidence_index){ //ensure that the coincidence we are merging is not the same
                    for (Int_t l=0; l<configurations[k].size(); l++){ // loop over the elements of the coincidence to search for possible candidates
                        if (std::get<0>(std::get<0>(configurations[k][l])) == std::get<0>(cathode_signals[i]) && std::get<1>(std::get<0>(configurations[k][l])) == std::get<1>(cathode_signals[i]) && std::get<2>(std::get<0>(configurations[k][l])) == std::get<2>(cathode_signals[i])) {                                    
                            if (std::get<8>(std::get<0>(configurations[k][l]))>300 && std::get<0>(cathode_signals[i]) == std::get<0>(std::get<0>(configurations[k][l])) && std::get<1>(std::get<0>(configurations[k][l])) == std::get<1>(cathode_signals[i]) && std::get<2>(std::get<0>(configurations[k][l])) == std::get<2>(cathode_signals[i])) {
                                if (std::get<3>(cathode_signals[i]) == std::get<3>(std::get<0>(configurations[k][l])) && std::get<4>(cathode_signals[i]) == std::get<4>(std::get<0>(configurations[k][l])) && std::get<5>(cathode_signals[i]) == std::get<5>(std::get<0>(configurations[k][l])) && std::get<6>(cathode_signals[i])== std::get<6>(std::get<0>(configurations[k][l]))) {
                                   for (size_t m=0; m<configurations[coincidence_index].size();m++){
                                        if (std::get<2>(configurations[coincidence_index][m])==std::get<2>(configurations[k][l]) && std::get<1>(configurations[coincidence_index][m])== std::get<1>(configurations[k][l])){ //check if the signal is already in the coincidence
                                            break;                                                
                                        }                                        
                                    }
                                    if (abs(std::get<7>(cathode_signals[i])- std::get<7>(std::get<0>(configurations[k][l]))) < 20) {
                                        for (size_t m=0; m<configurations[k].size(); m++){
                                            configurations[coincidence_index].push_back(configurations[k][m]); //add the signals from the other coincidence to the first one
                                        }
                                        configurations.erase(configurations.begin()+k);
                                        break; //remove the element merged to avoid double counting
                                    }
                                }
                            }                                
                        }
                    }                                
                }
            }
        }
    }
    // Print the entire configurations vector
std::cout << "INFO: number of coincidences in the second step: " << configurations.size() << std::endl; // print number of coincidences
std::cout << "INFO: Printing all configurations:" << std::endl;
for (size_t i = 0; i < configurations.size(); ++i) {
    std::cout << "Configuration " << i + 1 << ":" << std::endl;
    std::cout << "  Number of coincidences: " << configurations[i].size() << std::endl;

    // Sort the configuration by cathode detn
    std::sort(configurations[i].begin(), configurations[i].end(), [](const auto& a, const auto& b) {
        return std::get<6>(std::get<0>(a)) < std::get<6>(std::get<0>(b));
    });

    for (size_t j = 0; j < configurations[i].size(); ++j) {
        const auto& cathode_signal = std::get<0>(configurations[i][j]);
        Double_t anode_tof = std::get<1>(configurations[i][j]);
        Int_t detector_index = std::get<2>(configurations[i][j]);

        std::cout << "  Cathode Signal: "
                  << "RunNumber=" << std::get<0>(cathode_signal) << ", "
                  << "time=" << std::get<1>(cathode_signal) << ", "
                  << "psTime=" << std::get<2>(cathode_signal) << ", "
                  << "BunchNumber=" << std::get<3>(cathode_signal) << ", "
                  << "PSpulse=" << std::get<4>(cathode_signal) << ", "
                  << "PulseIntensity=" << std::get<5>(cathode_signal) << ", "
                  << "detn=" << std::get<6>(cathode_signal) << ", "
                  << "tof=" << std::get<7>(cathode_signal)-anode_tof << ", "
                  << "amp=" << std::get<8>(cathode_signal) << std::endl;
    }
}
}
           
