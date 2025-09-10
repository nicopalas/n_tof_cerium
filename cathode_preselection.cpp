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
Int_t max_mul = 5; //max number of signals in coincidence in a single detector
Int_t mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
Double_t tof0[5],tof1[5],tof2[5],tof3[5],tof4[5],tof5[5],tof6[5],tof7[5],tof8[5],tof9[5]; //must be hardcoded if arrays are defined at file scope
Float_t amp0[5],amp1[5],amp2[5],amp3[5],amp4[5],amp5[5],amp6[5],amp7[5],amp8[5],amp9[5];  //must be hardcoded if arrays are defined at file scope
Double_t gamma_flash0[5], gamma_flash1[5], gamma_flash2[5], gamma_flash3[5], gamma_flash4[5], gamma_flash5[5], gamma_flash6[5], gamma_flash7[5], gamma_flash8[5], gamma_flash9[5];
Int_t detn_cathode[40]; 
Double_t first_tof01[100], first_tof02[100], first_tof03[100], first_tof04[100], first_tof11[100], first_tof12[100], first_tof13[100], first_tof14[100], first_tof21[100], first_tof22[100], first_tof23[100], first_tof24[100], first_tof31[100], first_tof32[100], first_tof33[100], first_tof34[100], first_tof41[100], first_tof42[100], first_tof43[100], first_tof44[100], first_tof51[100], first_tof52[100], first_tof53[100], first_tof54[100], first_tof61[100], first_tof62[100], first_tof63[100], first_tof64[100], first_tof71[100], first_tof72[100], first_tof73[100], first_tof74[100], first_tof81[100], first_tof82[100], first_tof83[100], first_tof84[100], first_tof91[100], first_tof92[100], first_tof93[100], first_tof94[100];
Float_t first_amp01[100], first_amp02[100], first_amp03[100], first_amp04[100], first_amp11[100], first_amp12[100], first_amp13[100], first_amp14[100], first_amp21[100], first_amp22[100], first_amp23[100], first_amp24[100], first_amp31[100], first_amp32[100], first_amp33[100], first_amp34[100], first_amp41[100], first_amp42[100], first_amp43[100], first_amp44[100], first_amp51[100], first_amp52[100], first_amp53[100], first_amp54[100], first_amp61[100], first_amp62[100], first_amp63[100], first_amp64[100], first_amp71[100], first_amp72[100], first_amp73[100], first_amp74[100], first_amp81[100], first_amp82[100], first_amp83[100], first_amp84[100], first_amp91[100], first_amp92[100], first_amp93[100], first_amp94[100];
Int_t first_mult01 = 0, first_mult02 = 0, first_mult03 = 0, first_mult04 = 0, first_mult11 = 0, first_mult12 = 0, first_mult13 = 0, first_mult14 = 0, first_mult21 = 0, first_mult22 = 0, first_mult23 = 0, first_mult24 = 0, first_mult31 = 0, first_mult32 = 0, first_mult33 = 0, first_mult34 = 0, first_mult41 = 0, first_mult42 = 0, first_mult43 = 0, first_mult44 = 0, first_mult51 = 0, first_mult52 = 0, first_mult53 = 0, first_mult54 = 0, first_mult61 = 0, first_mult62 = 0, first_mult63 = 0, first_mult64 = 0, first_mult71 = 0, first_mult72 = 0, first_mult73 = 0, first_mult74 = 0, first_mult81 = 0, first_mult82 = 0, first_mult83 = 0, first_mult84 = 0, first_mult91 = 0, first_mult92 = 0, first_mult93 = 0, first_mult94 = 0;
Double_t neutron_energy;
Bool_t addCoincidences(Int_t the_detn) {
    if(the_detn<0 || the_detn>94 || the_detn%10>4) {        
        std::cout << "ERROR in addCoincidences: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(the_detn==01) { first_mult01++; detn_cathode[0]++; }
    if(the_detn==02) { first_mult02++; detn_cathode[1]++; }
    if(the_detn==03) { first_mult03++; detn_cathode[2]++; }
    if(the_detn==04) { first_mult04++; detn_cathode[3]++; }
    if(the_detn==11) { first_mult11++; detn_cathode[4]++; }
    if(the_detn==12) { first_mult12++; detn_cathode[5]++; }
    if(the_detn==13) { first_mult13++; detn_cathode[6]++; }
    if(the_detn==14) { first_mult14++; detn_cathode[7]++; }
    if(the_detn==21) { first_mult21++; detn_cathode[8]++; }
    if(the_detn==22) { first_mult22++; detn_cathode[9]++; }
    if(the_detn==23) { first_mult23++; detn_cathode[10]++; }
    if(the_detn==24) { first_mult24++; detn_cathode[11]++; }
    if(the_detn==31) { first_mult31++; detn_cathode[12]++; }
    if(the_detn==32) { first_mult32++; detn_cathode[13]++; }
    if(the_detn==33) { first_mult33++; detn_cathode[14]++; }
    if(the_detn==34) { first_mult34++; detn_cathode[15]++; }
    if(the_detn==41) { first_mult41++; detn_cathode[16]++; }
    if(the_detn==42) { first_mult42++; detn_cathode[17]++; }
    if(the_detn==43) { first_mult43++; detn_cathode[18]++; }
    if(the_detn==44) { first_mult44++; detn_cathode[19]++; }
    if(the_detn==51) { first_mult51++; detn_cathode[20]++; }
    if(the_detn==52) { first_mult52++; detn_cathode[21]++; }
    if(the_detn==53) { first_mult53++; detn_cathode[22]++; }
    if(the_detn==54) { first_mult54++; detn_cathode[23]++; }
    if(the_detn==61) { first_mult61++; detn_cathode[24]++; }
    if(the_detn==62) { first_mult62++; detn_cathode[25]++; }
    if(the_detn==63) { first_mult63++; detn_cathode[26]++; }
    if(the_detn==64) { first_mult64++; detn_cathode[27]++; }
    if(the_detn==71) { first_mult71++; detn_cathode[28]++; }
    if(the_detn==72) { first_mult72++; detn_cathode[29]++; }
    if(the_detn==73) { first_mult73++; detn_cathode[30]++; }
    if(the_detn==74) { first_mult74++; detn_cathode[31]++; }
    if(the_detn==81) { first_mult81++; detn_cathode[32]++; }
    if(the_detn==82) { first_mult82++; detn_cathode[33]++; }
    if(the_detn==83) { first_mult83++; detn_cathode[34]++; }
    if(the_detn==84) { first_mult84++; detn_cathode[35]++; }
    if(the_detn==91) { first_mult91++; detn_cathode[36]++; }
    if(the_detn==92) { first_mult92++; detn_cathode[37]++; }
    if(the_detn==93) { first_mult93++; detn_cathode[38]++; }
    if(the_detn==94) { first_mult94++; detn_cathode[39]++; }
    return true;
}

Bool_t addTof(Int_t the_detn, Double_t the_tof) {
    Int_t detn_index = -1;
    switch (the_detn) {
        case 01: detn_index = 0; break;
        case 02: detn_index = 1; break;
        case 03: detn_index = 2; break;
        case 04: detn_index = 3; break;
        case 11: detn_index = 4; break;
        case 12: detn_index = 5; break;
        case 13: detn_index = 6; break;
        case 14: detn_index = 7; break;
        case 21: detn_index = 8; break;
        case 22: detn_index = 9; break;
        case 23: detn_index = 10; break;
        case 24: detn_index = 11; break;
        case 31: detn_index = 12; break;
        case 32: detn_index = 13; break;
        case 33: detn_index = 14; break;
        case 34: detn_index = 15; break;
        case 41: detn_index = 16; break;
        case 42: detn_index = 17; break;
        case 43: detn_index = 18; break;
        case 44: detn_index = 19; break;
        case 51: detn_index = 20; break;
        case 52: detn_index = 21; break;
        case 53: detn_index = 22; break;
        case 54: detn_index = 23; break;
        case 61: detn_index = 24; break;
        case 62: detn_index = 25; break;
        case 63: detn_index = 26; break;
        case 64: detn_index = 27; break;
        case 71: detn_index = 28; break;
        case 72: detn_index = 29; break;
        case 73: detn_index = 30; break;
        case 74: detn_index = 31; break;
        case 81: detn_index = 32; break;
        case 82: detn_index = 33; break;
        case 83: detn_index = 34; break;
        case 84: detn_index = 35; break;
        case 91: detn_index = 36; break;
        case 92: detn_index = 37; break;
        case 93: detn_index = 38; break;
        case 94: detn_index = 39; break;
        default:
            std::cout << "ERROR in addTof: Unexpected the_detn value: " << the_detn << std::endl;
            return false;
    }
    if(detn_cathode[detn_index] < 1 || detn_cathode[detn_index] > 100) {         
        std::cout << "ERROR in addTof: detn_all[" << detn_index << "] value: " << detn_cathode[detn_index] << " not valid" << std::endl;
        return false;
    }
    switch (the_detn) {
        case 01: if (first_mult01 > 0 && first_mult01 < 100) first_tof01[first_mult01 - 1] = the_tof; break;
        case 02: if (first_mult02 > 0 && first_mult02 < 100) first_tof02[first_mult02 - 1] = the_tof; break;
        case 03: if (first_mult03 > 0 && first_mult03 < 100) first_tof03[first_mult03 - 1] = the_tof; break;
        case 04: if (first_mult04 > 0 && first_mult04 < 100) first_tof04[first_mult04 - 1] = the_tof; break;
        case 11: if (first_mult11 > 0 && first_mult11 < 100) first_tof11[first_mult11 - 1] = the_tof; break;
        case 12: if (first_mult12 > 0 && first_mult12 < 100) first_tof12[first_mult12 - 1] = the_tof; break;
        case 13: if (first_mult13 > 0 && first_mult13 < 100) first_tof13[first_mult13 - 1] = the_tof; break;
        case 14: if (first_mult14 > 0 && first_mult14 < 100) first_tof14[first_mult14 - 1] = the_tof; break;
        case 21: if (first_mult21 > 0 && first_mult21 < 100) first_tof21[first_mult21 - 1] = the_tof; break;
        case 22: if (first_mult22 > 0 && first_mult22 < 100) first_tof22[first_mult22 - 1] = the_tof; break;
        case 23: if (first_mult23 > 0 && first_mult23 < 100) first_tof23[first_mult23 - 1] = the_tof; break;
        case 24: if (first_mult24 > 0 && first_mult24 < 100) first_tof24[first_mult24 - 1] = the_tof; break;
        case 31: if (first_mult31 > 0 && first_mult31 < 100) first_tof31[first_mult31 - 1] = the_tof; break;
        case 32: if (first_mult32 > 0 && first_mult32 < 100) first_tof32[first_mult32 - 1] = the_tof; break;
        case 33: if (first_mult33 > 0 && first_mult33 < 100) first_tof33[first_mult33 - 1] = the_tof; break;
        case 34: if (first_mult34 > 0 && first_mult34 < 100) first_tof34[first_mult34 - 1] = the_tof; break;
        case 41: if (first_mult41 > 0 && first_mult41 < 100) first_tof41[first_mult41 - 1] = the_tof; break;
        case 42: if (first_mult42 > 0 && first_mult42 < 100) first_tof42[first_mult42 - 1] = the_tof; break;
        case 43: if (first_mult43 > 0 && first_mult43 < 100) first_tof43[first_mult43 - 1] = the_tof; break;
        case 44: if (first_mult44 > 0 && first_mult44 < 100) first_tof44[first_mult44 - 1] = the_tof; break;
        case 51: if (first_mult51 > 0 && first_mult51 < 100) first_tof51[first_mult51 - 1] = the_tof; break;
        case 52: if (first_mult52 > 0 && first_mult52 < 100) first_tof52[first_mult52 - 1] = the_tof; break;
        case 53: if (first_mult53 > 0 && first_mult53 < 100) first_tof53[first_mult53 - 1] = the_tof; break;
        case 54: if (first_mult54 > 0 && first_mult54 < 100) first_tof54[first_mult54 - 1] = the_tof; break;
        case 61: if (first_mult61 > 0 && first_mult61 < 100) first_tof61[first_mult61 - 1] = the_tof; break;
        case 62: if (first_mult62 > 0 && first_mult62 < 100) first_tof62[first_mult62 - 1] = the_tof; break;
        case 63: if (first_mult63 > 0 && first_mult63 < 100) first_tof63[first_mult63 - 1] = the_tof; break;
        case 64: if (first_mult64 > 0 && first_mult64 < 100) first_tof64[first_mult64 - 1] = the_tof; break;
        case 71: if (first_mult71 > 0 && first_mult71 < 100) first_tof71[first_mult71 - 1] = the_tof; break;
        case 72: if (first_mult72 > 0 && first_mult72 < 100) first_tof72[first_mult72 - 1] = the_tof; break;
        case 73: if (first_mult73 > 0 && first_mult73 < 100) first_tof73[first_mult73 - 1] = the_tof; break;
        case 74: if (first_mult74 > 0 && first_mult74 < 100) first_tof74[first_mult74 - 1] = the_tof; break;
        case 81: if (first_mult81 > 0 && first_mult81 < 100) first_tof81[first_mult81 - 1] = the_tof; break;
        case 82: if (first_mult82 > 0 && first_mult82 < 100) first_tof82[first_mult82 - 1] = the_tof; break;
        case 83: if (first_mult83 > 0 && first_mult83 < 100) first_tof83[first_mult83 - 1] = the_tof; break;
        case 84: if (first_mult84 > 0 && first_mult84 < 100) first_tof84[first_mult84 - 1] = the_tof; break;
        case 91: if (first_mult91 > 0 && first_mult91 < 100) first_tof91[first_mult91 - 1] = the_tof; break;
        case 92: if (first_mult92 > 0 && first_mult92 < 100) first_tof92[first_mult92 - 1] = the_tof; break;
        case 93: if (first_mult93 > 0 && first_mult93 < 100) first_tof93[first_mult93 - 1] = the_tof; break;
        case 94: if (first_mult94 > 0 && first_mult94 < 100) first_tof94[first_mult94 - 1] = the_tof; break;
        default:
            std::cout << "ERROR in addTof: Unexpected the_detn value: " << the_detn << std::endl;
            return false;
    }
    return true;
}

Bool_t addAmp(Int_t the_detn, Float_t the_amp) {
    Int_t detn_index = -1;
    switch (the_detn) {
        case 01: detn_index = 0; break;
        case 02: detn_index = 1; break;
        case 03: detn_index = 2; break;
        case 04: detn_index = 3; break;
        case 11: detn_index = 4; break;
        case 12: detn_index = 5; break;
        case 13: detn_index = 6; break;
        case 14: detn_index = 7; break;
        case 21: detn_index = 8; break;
        case 22: detn_index = 9; break;
        case 23: detn_index = 10; break;
        case 24: detn_index = 11; break;
        case 31: detn_index = 12; break;
        case 32: detn_index = 13; break;
        case 33: detn_index = 14; break;
        case 34: detn_index = 15; break;
        case 41: detn_index = 16; break;
        case 42: detn_index = 17; break;
        case 43: detn_index = 18; break;
        case 44: detn_index = 19; break;
        case 51: detn_index = 20; break;
        case 52: detn_index = 21; break;
        case 53: detn_index = 22; break;
        case 54: detn_index = 23; break;
        case 61: detn_index = 24; break;
        case 62: detn_index = 25; break;
        case 63: detn_index = 26; break;
        case 64: detn_index = 27; break;
        case 71: detn_index = 28; break;
        case 72: detn_index = 29; break;
        case 73: detn_index = 30; break;
        case 74: detn_index = 31; break;
        case 81: detn_index = 32; break;
        case 82: detn_index = 33; break;
        case 83: detn_index = 34; break;
        case 84: detn_index = 35; break;
        case 91: detn_index = 36; break;
        case 92: detn_index = 37; break;
        case 93: detn_index = 38; break;
        case 94: detn_index = 39; break;
        default:
            std::cout << "ERROR in addAmp: Unexpected the_detn value: " << the_detn << std::endl;
            return false;
    }
    if(detn_cathode[detn_index] < 1 || detn_cathode[detn_index] > 100) {        
        std::cout << "ERROR in addAmp: detn_cathode[" << detn_index << "] value: " << detn_cathode[detn_index] << " not valid" << std::endl;
        return false;
    }
    switch (the_detn) {
        case 01: if (first_mult01 > 0 && first_mult01 < 100) first_amp01[first_mult01 - 1] = the_amp; break;
        case 02: if (first_mult02 > 0 && first_mult02 < 100) first_amp02[first_mult02 - 1] = the_amp; break;
        case 03: if (first_mult03 > 0 && first_mult03 < 100) first_amp03[first_mult03 - 1] = the_amp; break;
        case 04: if (first_mult04 > 0 && first_mult04 < 100) first_amp04[first_mult04 - 1] = the_amp; break;
        case 11: if (first_mult11 > 0 && first_mult11 < 100) first_amp11[first_mult11 - 1] = the_amp; break;
        case 12: if (first_mult12 > 0 && first_mult12 < 100) first_amp12[first_mult12 - 1] = the_amp; break;
        case 13: if (first_mult13 > 0 && first_mult13 < 100) first_amp13[first_mult13 - 1] = the_amp; break;
        case 14: if (first_mult14 > 0 && first_mult14 < 100) first_amp14[first_mult14 - 1] = the_amp; break;
        case 21: if (first_mult21 > 0 && first_mult21 < 100) first_amp21[first_mult21 - 1] = the_amp; break;
        case 22: if (first_mult22 > 0 && first_mult22 < 100) first_amp22[first_mult22 - 1] = the_amp; break;
        case 23: if (first_mult23 > 0 && first_mult23 < 100) first_amp23[first_mult23 - 1] = the_amp; break;
        case 24: if (first_mult24 > 0 && first_mult24 < 100) first_amp24[first_mult24 - 1] = the_amp; break;
        case 31: if (first_mult31 > 0 && first_mult31 < 100) first_amp31[first_mult31 - 1] = the_amp; break;
        case 32: if (first_mult32 > 0 && first_mult32 < 100) first_amp32[first_mult32 - 1] = the_amp; break;
        case 33: if (first_mult33 > 0 && first_mult33 < 100) first_amp33[first_mult33 - 1] = the_amp; break;
        case 34: if (first_mult34 > 0 && first_mult34 < 100) first_amp34[first_mult34 - 1] = the_amp; break;
        case 41: if (first_mult41 > 0 && first_mult41 < 100) first_amp41[first_mult41 - 1] = the_amp; break;
        case 42: if (first_mult42 > 0 && first_mult42 < 100) first_amp42[first_mult42 - 1] = the_amp; break;
        case 43: if (first_mult43 > 0 && first_mult43 < 100) first_amp43[first_mult43 - 1] = the_amp; break;
        case 44: if (first_mult44 > 0 && first_mult44 < 100) first_amp44[first_mult44 - 1] = the_amp; break;
        case 51: if (first_mult51 > 0 && first_mult51 < 100) first_amp51[first_mult51 - 1] = the_amp; break;
        case 52: if (first_mult52 > 0 && first_mult52 < 100) first_amp52[first_mult52 - 1] = the_amp; break;
        case 53: if (first_mult53 > 0 && first_mult53 < 100) first_amp53[first_mult53 - 1] = the_amp; break;
        case 54: if (first_mult54 > 0 && first_mult54 < 100) first_amp54[first_mult54 - 1] = the_amp; break;
        case 61: if (first_mult61 > 0 && first_mult61 < 100) first_amp61[first_mult61 - 1] = the_amp; break;
        case 62: if (first_mult62 > 0 && first_mult62 < 100) first_amp62[first_mult62 - 1] = the_amp; break;
        case 63: if (first_mult63 > 0 && first_mult63 < 100) first_amp63[first_mult63 - 1] = the_amp; break;
        case 64: if (first_mult64 > 0 && first_mult64 < 100) first_amp64[first_mult64 - 1] = the_amp; break;
        case 71: if (first_mult71 > 0 && first_mult71 < 100) first_amp71[first_mult71 - 1] = the_amp; break;
        case 72: if (first_mult72 > 0 && first_mult72 < 100) first_amp72[first_mult72 - 1] = the_amp; break;
        case 73: if (first_mult73 > 0 && first_mult73 < 100) first_amp73[first_mult73 - 1] = the_amp; break;
        case 74: if (first_mult74 > 0 && first_mult74 < 100) first_amp74[first_mult74 - 1] = the_amp; break;
        case 81: if (first_mult81 > 0 && first_mult81 < 100) first_amp81[first_mult81 - 1] = the_amp; break;
        case 82: if (first_mult82 > 0 && first_mult82 < 100) first_amp82[first_mult82 - 1] = the_amp; break;
        case 83: if (first_mult83 > 0 && first_mult83 < 100) first_amp83[first_mult83 - 1] = the_amp; break;
        case 84: if (first_mult84 > 0 && first_mult84 < 100) first_amp84[first_mult84 - 1] = the_amp; break;
        case 91: if (first_mult91 > 0 && first_mult91 < 100) first_amp91[first_mult91 - 1] = the_amp; break;
        case 92: if (first_mult92 > 0 && first_mult92 < 100) first_amp92[first_mult92 - 1] = the_amp; break;
        case 93: if (first_mult93 > 0 && first_mult93 < 100) first_amp93[first_mult93 - 1] = the_amp; break;
        case 94: if (first_mult94 > 0 && first_mult94 < 100) first_amp94[first_mult94 - 1] = the_amp; break;
        default:
            std::cout << "ERROR in addAmp: Unexpected the_detn value: " << the_detn << std::endl;
            return false;
    }
    return true;
}
void cathode_preselection(Int_t run_number, Int_t time_for_coincidence, Float_t amp_threshold) {
    std::string filename = "/nucl_lustre/n_tof_INTC_P_665/DATA/run" + std::to_string(run_number) + ".root";
    TFile *file = TFile::Open(filename.c_str());// data from the pickup and the cathode
    TTree *cathodeTree = (TTree*)file->Get("PPAC");
    TTree *pickupTree = (TTree*)file->Get("PKUP");


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

    Long64_t ncathodeEntries = cathodeTree->GetEntries();
    Long64_t npickupEntries = pickupTree->GetEntries();
    std::cout << "INFO: Processing run_number " << run_number << " with " << ncathodeEntries << " entries" << std::endl;

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
    pickupTree->SetBranchStatus("amp", 1); // Float_t
    pickupTree->SetBranchAddress("RunNumber", &pickup_RunNumber);
    pickupTree->SetBranchAddress("time", &pickup_time);
    pickupTree->SetBranchAddress("psTime", &pickup_psTime);
    pickupTree->SetBranchAddress("BunchNumber", &pickup_BunchNumber);
    pickupTree->SetBranchAddress("PSpulse", &pickup_PSpulse);
    pickupTree->SetBranchAddress("PulseIntensity", &pickup_PulseIntensity);
    pickupTree->SetBranchAddress("tof", &pickup_tof);
    pickupTree->SetBranchAddress("amp", &pickup_amp);

    // Create output file and tree
    std::string out_filename = "/nucl_lustre/n_tof_INTC_P_665/Analysis/Output/Cathodes/preselection/collection_cathode" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(out_filename.c_str(), "RECREATE");
    TTree *outtree = new TTree("cathode_preselection", "cathode_preselection");
    outtree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    outtree->Branch("BunchNumber", &BunchNumber, "BunchNumber/I");
    outtree->Branch("PSpulse", &PSpulse, "PSpulse/I");
    outtree->Branch("time", &eventTime, "time/I");
    outtree->Branch("psTime", &psTime, "psTime/D");
    outtree->Branch("PulseIntensity", &PulseIntensity, "PulseIntensity/F");
    outtree->Branch("detn_cathode", detn_cathode, "detn_cathode[40]/I");
    outtree->Branch("first_mult01", &first_mult01, "first_mult01/I");
    outtree->Branch("first_mult02", &first_mult02, "first_mult02/I");
    outtree->Branch("first_mult03", &first_mult03, "first_mult03/I");
    outtree->Branch("first_mult04", &first_mult04, "first_mult04/I");
    outtree->Branch("first_mult11", &first_mult11, "first_mult11/I");
    outtree->Branch("first_mult12", &first_mult12, "first_mult12/I");
    outtree->Branch("first_mult13", &first_mult13, "first_mult13/I");
    outtree->Branch("first_mult14", &first_mult14, "first_mult14/I");
    outtree->Branch("first_mult21", &first_mult21, "first_mult21/I");
    outtree->Branch("first_mult22", &first_mult22, "first_mult22/I");
    outtree->Branch("first_mult23", &first_mult23, "first_mult23/I");
    outtree->Branch("first_mult24", &first_mult24, "first_mult24/I");
    outtree->Branch("first_mult31", &first_mult31, "first_mult31/I");
    outtree->Branch("first_mult32", &first_mult32, "first_mult32/I");
    outtree->Branch("first_mult33", &first_mult33, "first_mult33/I");
    outtree->Branch("first_mult34", &first_mult34, "first_mult34/I");
    outtree->Branch("first_mult41", &first_mult41, "first_mult41/I");
    outtree->Branch("first_mult42", &first_mult42, "first_mult42/I");
    outtree->Branch("first_mult43", &first_mult43, "first_mult43/I");
    outtree->Branch("first_mult44", &first_mult44, "first_mult44/I");
    outtree->Branch("first_mult51", &first_mult51, "first_mult51/I");
    outtree->Branch("first_mult52", &first_mult52, "first_mult52/I");
    outtree->Branch("first_mult53", &first_mult53, "first_mult53/I");
    outtree->Branch("first_mult54", &first_mult54, "first_mult54/I");
    outtree->Branch("first_mult61", &first_mult61, "first_mult61/I");
    outtree->Branch("first_mult62", &first_mult62, "first_mult62/I");
    outtree->Branch("first_mult63", &first_mult63, "first_mult63/I");
    outtree->Branch("first_mult64", &first_mult64, "first_mult64/I");
    outtree->Branch("first_mult71", &first_mult71, "first_mult71/I");
    outtree->Branch("first_mult72", &first_mult72, "first_mult72/I");
    outtree->Branch("first_mult73", &first_mult73, "first_mult73/I");
    outtree->Branch("first_mult74", &first_mult74, "first_mult74/I");
    outtree->Branch("first_mult81", &first_mult81, "first_mult81/I");
    outtree->Branch("first_mult82", &first_mult82, "first_mult82/I");
    outtree->Branch("first_mult83", &first_mult83, "first_mult83/I");
    outtree->Branch("first_mult84", &first_mult84, "first_mult84/I");
    outtree->Branch("first_mult91", &first_mult91, "first_mult91/I");
    outtree->Branch("first_mult92", &first_mult92, "first_mult92/I");
    outtree->Branch("first_mult93", &first_mult93, "first_mult93/I");
    outtree->Branch("first_mult94", &first_mult94, "first_mult94/I");
    outtree->Branch("first_tof01", first_tof01, "first_tof01[first_mult01]/D");
    outtree->Branch("first_tof02", first_tof02, "first_tof02[first_mult02]/D");
    outtree->Branch("first_tof03", first_tof03, "first_tof03[first_mult03]/D");
    outtree->Branch("first_tof04", first_tof04, "first_tof04[first_mult04]/D");
    outtree->Branch("first_tof11", first_tof11, "first_tof11[first_mult11]/D");
    outtree->Branch("first_tof12", first_tof12, "first_tof12[first_mult12]/D");
    outtree->Branch("first_tof13", first_tof13, "first_tof13[first_mult13]/D");
    outtree->Branch("first_tof14", first_tof14, "first_tof14[first_mult14]/D");
    outtree->Branch("first_tof21", first_tof21, "first_tof21[first_mult21]/D");
    outtree->Branch("first_tof22", first_tof22, "first_tof22[first_mult22]/D");
    outtree->Branch("first_tof23", first_tof23, "first_tof23[first_mult23]/D");
    outtree->Branch("first_tof24", first_tof24, "first_tof24[first_mult24]/D");
    outtree->Branch("first_tof31", first_tof31, "first_tof31[first_mult31]/D");
    outtree->Branch("first_tof32", first_tof32, "first_tof32[first_mult32]/D");
    outtree->Branch("first_tof33", first_tof33, "first_tof33[first_mult33]/D");
    outtree->Branch("first_tof34", first_tof34, "first_tof34[first_mult34]/D");
    outtree->Branch("first_tof41", first_tof41, "first_tof41[first_mult41]/D");
    outtree->Branch("first_tof42", first_tof42, "first_tof42[first_mult42]/D");
    outtree->Branch("first_tof43", first_tof43, "first_tof43[first_mult43]/D");
    outtree->Branch("first_tof44", first_tof44, "first_tof44[first_mult44]/D");
    outtree->Branch("first_tof51", first_tof51, "first_tof51[first_mult51]/D");
    outtree->Branch("first_tof52", first_tof52, "first_tof52[first_mult52]/D");
    outtree->Branch("first_tof53", first_tof53, "first_tof53[first_mult53]/D");
    outtree->Branch("first_tof54", first_tof54, "first_tof54[first_mult54]/D");
    outtree->Branch("first_tof61", first_tof61, "first_tof61[first_mult61]/D");
    outtree->Branch("first_tof62", first_tof62, "first_tof62[first_mult62]/D");
    outtree->Branch("first_tof63", first_tof63, "first_tof63[first_mult63]/D");
    outtree->Branch("first_tof64", first_tof64, "first_tof64[first_mult64]/D");
    outtree->Branch("first_tof71", first_tof71, "first_tof71[first_mult71]/D");
    outtree->Branch("first_tof72", first_tof72, "first_tof72[first_mult72]/D");
    outtree->Branch("first_tof73", first_tof73, "first_tof73[first_mult73]/D");
    outtree->Branch("first_tof74", first_tof74, "first_tof74[first_mult74]/D");
    outtree->Branch("first_tof81", first_tof81, "first_tof81[first_mult81]/D");
    outtree->Branch("first_tof82", first_tof82, "first_tof82[first_mult82]/D");
    outtree->Branch("first_tof83", first_tof83, "first_tof83[first_mult83]/D");
    outtree->Branch("first_tof84", first_tof84, "first_tof84[first_mult84]/D");
    outtree->Branch("first_tof91", first_tof91, "first_tof91[first_mult91]/D");
    outtree->Branch("first_tof92", first_tof92, "first_tof92[first_mult92]/D");
    outtree->Branch("first_tof93", first_tof93, "first_tof93[first_mult93]/D");
    outtree->Branch("first_tof94", first_tof94, "first_tof94[first_mult94]/D");
    outtree->Branch("first_amp01", first_amp01, "first_amp01[first_mult01]/F");
    outtree->Branch("first_amp02", first_amp02, "first_amp02[first_mult02]/F");
    outtree->Branch("first_amp03", first_amp03, "first_amp03[first_mult03]/F");
    outtree->Branch("first_amp04", first_amp04, "first_amp04[first_mult04]/F");
    outtree->Branch("first_amp11", first_amp11, "first_amp11[first_mult11]/F");
    outtree->Branch("first_amp12", first_amp12, "first_amp12[first_mult12]/F");
    outtree->Branch("first_amp13", first_amp13, "first_amp13[first_mult13]/F");
    outtree->Branch("first_amp14", first_amp14, "first_amp14[first_mult14]/F");
    outtree->Branch("first_amp21", first_amp21, "first_amp21[first_mult21]/F");
    outtree->Branch("first_amp22", first_amp22, "first_amp22[first_mult22]/F");
    outtree->Branch("first_amp23", first_amp23, "first_amp23[first_mult23]/F");
    outtree->Branch("first_amp24", first_amp24, "first_amp24[first_mult24]/F");
    outtree->Branch("first_amp31", first_amp31, "first_amp31[first_mult31]/F");
    outtree->Branch("first_amp32", first_amp32, "first_amp32[first_mult32]/F");
    outtree->Branch("first_amp33", first_amp33, "first_amp33[first_mult33]/F");
    outtree->Branch("first_amp34", first_amp34, "first_amp34[first_mult34]/F");
    outtree->Branch("first_amp41", first_amp41, "first_amp41[first_mult41]/F");
    outtree->Branch("first_amp42", first_amp42, "first_amp42[first_mult42]/F");
    outtree->Branch("first_amp43", first_amp43, "first_amp43[first_mult43]/F");
    outtree->Branch("first_amp44", first_amp44, "first_amp44[first_mult44]/F");
    outtree->Branch("first_amp51", first_amp51, "first_amp51[first_mult51]/F");
    outtree->Branch("first_amp52", first_amp52, "first_amp52[first_mult52]/F");
    outtree->Branch("first_amp53", first_amp53, "first_amp53[first_mult53]/F");
    outtree->Branch("first_amp54", first_amp54, "first_amp54[first_mult54]/F");
    outtree->Branch("first_amp61", first_amp61, "first_amp61[first_mult61]/F");
    outtree->Branch("first_amp62", first_amp62, "first_amp62[first_mult62]/F");
    outtree->Branch("first_amp63", first_amp63, "first_amp63[first_mult63]/F");
    outtree->Branch("first_amp64", first_amp64, "first_amp64[first_mult64]/F");
    outtree->Branch("first_amp71", first_amp71, "first_amp71[first_mult71]/F");
    outtree->Branch("first_amp72", first_amp72, "first_amp72[first_mult72]/F");
    outtree->Branch("first_amp73", first_amp73, "first_amp73[first_mult73]/F");
    outtree->Branch("first_amp74", first_amp74, "first_amp74[first_mult74]/F");
    outtree->Branch("first_amp81", first_amp81, "first_amp81[first_mult81]/F");
    outtree->Branch("first_amp82", first_amp82, "first_amp82[first_mult82]/F");
    outtree->Branch("first_amp83", first_amp83, "first_amp83[first_mult83]/F");
    outtree->Branch("first_amp84", first_amp84, "first_amp84[first_mult84]/F");
    outtree->Branch("first_amp91", first_amp91, "first_amp91[first_mult91]/F");
    outtree->Branch("first_amp92", first_amp92, "first_amp92[first_mult92]/F");
    outtree->Branch("first_amp93", first_amp93, "first_amp93[first_mult93]/F");
    outtree->Branch("first_amp94", first_amp94, "first_amp94[first_mult94]/F");


    std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Float_t, Double_t>> signals;
    std::vector<std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>>> configurations;


    for (Long64_t i = 0; i < ncathodeEntries; ++i) {
        cathodeTree->GetEntry(i);
        signals.emplace_back(cathode_RunNumber, cathode_time, cathode_psTime, cathode_BunchNumber, cathode_PSpulse, cathode_PulseIntensity, cathode_detn, cathode_tof, cathode_amp); //insert tuple at the end

    }
    std::cout <<"INFO: number of signals: "<< signals.size() << std::endl; // print number of signals
    std::vector<bool> used(signals.size(), false); // track used signals to avoid double counting
    for (size_t i = 0; i < signals.size(); ++i) {
        for(size_t j = 0; j < npickupEntries; ++j){
            pickupTree->GetEntry(j);
            if (pickup_RunNumber == std::get<0>(signals[i]) && 
                pickup_time == std::get<1>(signals[i]) && 
                pickup_psTime == std::get<2>(signals[i]) && 
                pickup_BunchNumber == std::get<3>(signals[i]) && 
                pickup_PSpulse == std::get<4>(signals[i]) &&
                pickup_PulseIntensity== std::get<5>(signals[i])) {
                    std::get<7>(signals[i])=std::get<7>(signals[i])-pickup_tof;
                    std::cout << std::get<7>(signals[i]) << std::endl;
                    std::cout << "INFO: pickup signal found for cathode signal " << i << std::endl;
                }
            }
    }

    for (size_t i = 0; i < signals.size(); ++i) {
        if (used[i] == false && std::get<8>(signals[i]) > amp_threshold) { 
            std::vector<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t, Int_t, Double_t, Float_t>> coincidences;
            for (size_t j = i + 1; j <signals.size(); ++j) { 
                if (std::get<8>(signals[j]) > amp_threshold && used[j] == false) {
                    if (std::get<0>(signals[j]) == std::get<0>(signals[i]) && 
                        std::get<5>(signals[j]) == std::get<5>(signals[i]) && 
                        std::get<4>(signals[j]) == std::get<4>(signals[i])) { 
                        
                        if (std::get<3>(signals[j]) == std::get<3>(signals[i]) && 
                            std::get<2>(signals[j]) == std::get<2>(signals[i]) && 
                            std::get<1>(signals[j]) == std::get<1>(signals[i])) {
                            
                            if (abs(std::get<7>(signals[j]) - std::get<7>(signals[i])) < time_for_coincidence) {
                               

                                used[j] = true;
                                used[i] = true;

                                if (coincidences.size() == 0) {

                                    coincidences.push_back(signals[i]); 
                                }
                                coincidences.push_back(signals[j]);
                            }
                        }
                    }
                }
            }
            if (!coincidences.empty()) {
                configurations.push_back(coincidences); 
            }
        }
    }
    std::cout << "INFO: number of coincidences: " << configurations.size() << std::endl; // print number of coincidences
    

    for (size_t i = 0; i < configurations.size(); ++i) {
            RunNumber = std::get<0>(configurations[i][0]);
            eventTime = std::get<1>(configurations[i][0]);
            psTime = std::get<2>(configurations[i][0]);
            BunchNumber = std::get<3>(configurations[i][0]);
            PSpulse = std::get<4>(configurations[i][0]);
            PulseIntensity = std::get<5>(configurations[i][0]); //store every coincidence in the tree
        for (size_t j = 0; j < configurations[i].size(); ++j) {//loop over configurations and over each signal in the configurations
            
            cathode_detn = std::get<6>(configurations[i][j]);
            cathode_tof = std::get<7>(configurations[i][j]);
            cathode_amp = std::get<8>(configurations[i][j]);
            addCoincidences(cathode_detn);
            addTof(cathode_detn, cathode_tof);
            addAmp(cathode_detn, cathode_amp);
        }
        mult = configurations[i].size(); //mult will be the number of signals in coincidence

        outtree->Fill(); //fill the tree
    
        mult=1;
        first_mult01 = 0;
        first_mult02 = 0;
        first_mult03 = 0;
        first_mult04 = 0;
        first_mult11 = 0;
        first_mult12 = 0;
        first_mult13 = 0;
        first_mult14 = 0;
        first_mult21 = 0;
        first_mult22 = 0;
        first_mult23 = 0;
        first_mult24 = 0;
        first_mult31 = 0;
        first_mult32 = 0;
        first_mult33 = 0;
        first_mult34 = 0;
        first_mult41 = 0;
        first_mult42 = 0;
        first_mult43 = 0;
        first_mult44 = 0;
        first_mult51 = 0;
        first_mult52 = 0;
        first_mult53 = 0;
        first_mult54 = 0;
        first_mult61 = 0;
        first_mult62 = 0;
        first_mult63 = 0;
        first_mult64 = 0;
        first_mult71 = 0;
        first_mult72 = 0;
        first_mult73 = 0;
        first_mult74 = 0;
        first_mult81 = 0;
        first_mult82 = 0;
        first_mult83 = 0;
        first_mult84 = 0;
        first_mult91 = 0;
        first_mult92 = 0;
        first_mult93 = 0;
        first_mult94 = 0;
        std::fill(std::begin(detn_cathode), std::end(detn_cathode), 0);
    }
        outtree->Write();
        outfile->Close();
    }
    
