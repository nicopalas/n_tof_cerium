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

Int_t RunNumber, BunchNumber, PSpulse, detn, eventTime;
Float_t PulseIntensity;
Double_t psTime;
Double_t tof0[20], tof1[20], tof2[20], tof3[20], tof4[20], tof5[20];
Double_t tof6[20], tof7[20], tof8[20], tof9[20];
Float_t amp0[20], amp1[20], amp2[20], amp3[20], amp4[20];
Float_t amp5[20], amp6[20], amp7[20], amp8[20], amp9[20];
Int_t new_mult0, new_mult1, new_mult2, new_mult3, new_mult4;
Int_t new_mult5, new_mult6, new_mult7, new_mult8, new_mult9;
Int_t new_mult;
Int_t new_detn_all[10];

// Amplitudes (each has 20 elements)
Float_t new_amp0[20], new_amp1[20], new_amp2[20], new_amp3[20], new_amp4[20];
Float_t new_amp5[20], new_amp6[20], new_amp7[20], new_amp8[20], new_amp9[20];

// Times of flight (each has 20 elements)
Double_t new_tof0[20], new_tof1[20], new_tof2[20], new_tof3[20], new_tof4[20];
Double_t new_tof5[20], new_tof6[20], new_tof7[20], new_tof8[20], new_tof9[20];
Double_t neutron_energy_0, neutron_energy_1, neutron_energy_2, neutron_energy_3, neutron_energy_4;
Double_t neutron_energy_5, neutron_energy_6, neutron_energy_7, neutron_energy_8;
Int_t mult0, mult1, mult2, mult3, mult4, mult5, mult6, mult7, mult8, mult9;
Int_t mult;
Int_t mult01 = 0, mult02 = 0, mult03 = 0, mult04 = 0, mult11 = 0, mult12 = 0, mult13 = 0, mult14 = 0, mult21 = 0, mult22 = 0, mult23 = 0, mult24 = 0, mult31 = 0, mult32 = 0, mult33 = 0, mult34 = 0, mult41 = 0, mult42 = 0, mult43 = 0, mult44 = 0, mult51 = 0, mult52 = 0, mult53 = 0, mult54 = 0, mult61 = 0, mult62 = 0, mult63 = 0, mult64 = 0, mult71 = 0, mult72 = 0, mult73 = 0, mult74 = 0, mult81 = 0, mult82 = 0, mult83 = 0, mult84 = 0, mult91 = 0, mult92 = 0, mult93 = 0, mult94 = 0;
Int_t cathode_RunNumber, cathode_time, cathode_BunchNumber, cathode_PSpulse;
Double_t cathode_psTime;
Float_t cathode_PulseIntensity;
Int_t detn_cathode[40];
Int_t detn_all[10];
Double_t tof01[80], tof02[80], tof03[80], tof04[80], tof11[80], tof12[80], tof13[80], tof14[80], tof21[80], tof22[80], tof23[80], tof24[80], tof31[80], tof32[80], tof33[80], tof34[80], tof41[80], tof42[80], tof43[80], tof44[80], tof51[80], tof52[80], tof53[80], tof54[80], tof61[80], tof62[80], tof63[80], tof64[80], tof71[80], tof72[80], tof73[80], tof74[80], tof81[80], tof82[80], tof83[80], tof84[80], tof91[80], tof92[80], tof93[80], tof94[80];
Float_t amp01[80], amp02[80], amp03[80], amp04[80], amp11[80], amp12[80], amp13[80], amp14[80], amp21[80], amp22[80], amp23[80], amp24[80], amp31[80], amp32[80], amp33[80], amp34[80], amp41[80], amp42[80], amp43[80], amp44[80], amp51[80], amp52[80], amp53[80], amp54[80], amp61[80], amp62[80], amp63[80], amp64[80], amp71[80], amp72[80], amp73[80], amp74[80], amp81[80], amp82[80], amp83[80], amp84[80], amp91[80], amp92[80], amp93[80], amp94[80];

Bool_t addCoincidences_anodes(Int_t the_detn) {
    if(the_detn<0 || the_detn>9) {        
        std::cout << "ERROR in addTof: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(the_detn==0) { new_mult0++; new_detn_all[the_detn]++; }
    if(the_detn==1) { new_mult1++; new_detn_all[the_detn]++; }
    if(the_detn==2) { new_mult2++; new_detn_all[the_detn]++; }
    if(the_detn==3) { new_mult3++; new_detn_all[the_detn]++; }
    if(the_detn==4) { new_mult4++; new_detn_all[the_detn]++; }
    if(the_detn==5) { new_mult5++; new_detn_all[the_detn]++; }
    if(the_detn==6) { new_mult6++; new_detn_all[the_detn]++; }
    if(the_detn==7) { new_mult7++; new_detn_all[the_detn]++; }
    if(the_detn==8) { new_mult8++; new_detn_all[the_detn]++; }
    if(the_detn==9) { new_mult9++; new_detn_all[the_detn]++; }
    return true;
}

Bool_t addTof_anodes(Int_t the_detn, Double_t the_tof) {
    if(the_detn<0 || the_detn>9) {        
        std::cout << "ERROR in addTof: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(new_detn_all[the_detn] < 1 || new_detn_all[the_detn] > 20) {         //check the mutliplicity is valid
        std::cout << "ERROR in addTof: detn_all[" << the_detn <<"] value: " << new_detn_all[the_detn] << " not valid" << std::endl;
        for (int i = 0; i < new_detn_all[the_detn]; ++i) {
            if (the_detn == 0) std::cout << "tof0[" << i << "] = " << new_tof0[i] << std::endl;
            if (the_detn == 1) std::cout << "tof1[" << i << "] = " << new_tof1[i] << std::endl;
            if (the_detn == 2) std::cout << "tof2[" << i << "] = " << new_tof2[i] << std::endl;
            if (the_detn == 3) std::cout << "tof3[" << i << "] = " << new_tof3[i] << std::endl;
            if (the_detn == 4) std::cout << "tof4[" << i << "] = " << new_tof4[i] << std::endl;
            if (the_detn == 5) std::cout << "tof5[" << i << "] = " << new_tof5[i] << std::endl;
            if (the_detn == 6) std::cout << "tof6[" << i << "] = " << new_tof6[i] << std::endl;
            if (the_detn == 7) std::cout << "tof7[" << i << "] = " << new_tof7[i] << std::endl;
            if (the_detn == 8) std::cout << "tof8[" << i << "] = " << new_tof8[i] << std::endl;
            if (the_detn == 9) std::cout << "tof9[" << i << "] = " << new_tof9[i] << std::endl;
        }
        return false;
    }
    if(the_detn==0) new_tof0[new_mult0 - 1] = the_tof;
    if(the_detn==1) new_tof1[new_mult1 - 1] = the_tof;
    if(the_detn==2) new_tof2[new_mult2 - 1] = the_tof;
    if(the_detn==3) new_tof3[new_mult3 - 1] = the_tof;
    if(the_detn==4) new_tof4[new_mult4 - 1] = the_tof;
    if(the_detn==5) new_tof5[new_mult5 - 1] = the_tof;
    if(the_detn==6) new_tof6[new_mult6 - 1] = the_tof;
    if(the_detn==7) new_tof7[new_mult7 - 1] = the_tof;
    if(the_detn==8) new_tof8[new_mult8 - 1] = the_tof;
    if(the_detn==9) new_tof9[new_mult9 - 1] = the_tof;
    return true;
}

Bool_t addAmp_anodes(Int_t the_detn, Float_t the_amp) {
    if(the_detn<0 || the_detn>9) {        
        std::cout << "ERROR in addAmp: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(new_detn_all[the_detn] < 1 || new_detn_all[the_detn] > 20) {        
        std::cout << "ERROR in addAmp: detn_all[" << the_detn <<"] value: " << new_detn_all[the_detn] << " not valid" << std::endl;
        return false;
    }
    if(the_detn==0) new_amp0[new_mult0 - 1] = the_amp;
    if(the_detn==1) new_amp1[new_mult1 - 1] = the_amp;
    if(the_detn==2) new_amp2[new_mult2 - 1] = the_amp;
    if(the_detn==3) new_amp3[new_mult3 - 1] = the_amp;
    if(the_detn==4) new_amp4[new_mult4 - 1] = the_amp;
    if(the_detn==5) new_amp5[new_mult5 - 1] = the_amp;
    if(the_detn==6) new_amp6[new_mult6 - 1] = the_amp;
    if(the_detn==7) new_amp7[new_mult7 - 1] = the_amp;
    if(the_detn==8) new_amp8[new_mult8 - 1] = the_amp;
    if(the_detn==9) new_amp9[new_mult9 - 1] = the_amp;
    return true;
}

Bool_t addCoincidences(Int_t the_detn) {
    if(the_detn<0 || the_detn>94) {        
        std::cout << "ERROR in addCoincidences: the_detn value: " << the_detn << " not valid" << std::endl;
        return false;
    }
    if(the_detn==01) { mult01++; detn_cathode[0]++; }
    if(the_detn==02) { mult02++; detn_cathode[1]++; }
    if(the_detn==03) { mult03++; detn_cathode[2]++; }
    if(the_detn==04) { mult04++; detn_cathode[3]++; }
    if(the_detn==11) { mult11++; detn_cathode[4]++; }
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

Bool_t addTof(Int_t the_detn, Double_t the_tof) {
    Int_t detn_index;
    switch (the_detn) {
        case 01: if (mult01 > 0 && mult01 < 80) tof01[mult01 - 1] = the_tof; detn_index=0; break;
        case 02: if (mult02 > 0 && mult02 < 80) tof02[mult02 - 1] = the_tof; detn_index=1; break;
        case 03: if (mult03 > 0 && mult03 < 80) tof03[mult03 - 1] = the_tof; detn_index=2; break;
        case 04: if (mult04 > 0 && mult04 < 80) tof04[mult04 - 1] = the_tof; detn_index=3; break;
        case 11: if (mult11 > 0 && mult11 < 80) tof11[mult11 - 1] = the_tof; detn_index=4; break;
        case 12: if (mult12 > 0 && mult12 < 80) tof12[mult12 - 1] = the_tof; detn_index=5; break;
        case 13: if (mult13 > 0 && mult13 < 80) tof13[mult13 - 1] = the_tof; detn_index=6; break;
        case 14: if (mult14 > 0 && mult14 < 80) tof14[mult14 - 1] = the_tof; detn_index=7; break;
        case 21: if (mult21 > 0 && mult21 < 80) tof21[mult21 - 1] = the_tof; detn_index=8; break;
        case 22: if (mult22 > 0 && mult22 < 80) tof22[mult22 - 1] = the_tof; detn_index=9; break;
        case 23: if (mult23 > 0 && mult23 < 80) tof23[mult23 - 1] = the_tof; detn_index=10; break;
        case 24: if (mult24 > 0 && mult24 < 80) tof24[mult24 - 1] = the_tof; detn_index=11; break;
        case 31: if (mult31 > 0 && mult31 < 80) tof31[mult31 - 1] = the_tof; detn_index=12; break;
        case 32: if (mult32 > 0 && mult32 < 80) tof32[mult32 - 1] = the_tof; detn_index=13; break;
        case 33: if (mult33 > 0 && mult33 < 80) tof33[mult33 - 1] = the_tof; detn_index=14; break;
        case 34: if (mult34 > 0 && mult34 < 80) tof34[mult34 - 1] = the_tof; detn_index=15; break;
        case 41: if (mult41 > 0 && mult41 < 80) tof41[mult41 - 1] = the_tof; detn_index=16; break;
        case 42: if (mult42 > 0 && mult42 < 80) tof42[mult42 - 1] = the_tof; detn_index=17; break;
        case 43: if (mult43 > 0 && mult43 < 80) tof43[mult43 - 1] = the_tof; detn_index=18; break;
        case 44: if (mult44 > 0 && mult44 < 80) tof44[mult44 - 1] = the_tof; detn_index=19; break;
        case 51: if (mult51 > 0 && mult51 < 80) tof51[mult51 - 1] = the_tof; detn_index=20; break;
        case 52: if (mult52 > 0 && mult52 < 80) tof52[mult52 - 1] = the_tof; detn_index=21; break;
        case 53: if (mult53 > 0 && mult53 < 80) tof53[mult53 - 1] = the_tof; detn_index=22; break;
        case 54: if (mult54 > 0 && mult54 < 80) tof54[mult54 - 1] = the_tof; detn_index=23; break;
        case 61: if (mult61 > 0 && mult61 < 80) tof61[mult61 - 1] = the_tof; detn_index=24; break;
        case 62: if (mult62 > 0 && mult62 < 80) tof62[mult62 - 1] = the_tof; detn_index=25; break;
        case 63: if (mult63 > 0 && mult63 < 80) tof63[mult63 - 1] = the_tof; detn_index=26; break;
        case 64: if (mult64 > 0 && mult64 < 80) tof64[mult64 - 1] = the_tof; detn_index=27; break;
        case 71: if (mult71 > 0 && mult71 < 80) tof71[mult71 - 1] = the_tof; detn_index=28; break;
        case 72: if (mult72 > 0 && mult72 < 80) tof72[mult72 - 1] = the_tof; detn_index=29; break;
        case 73: if (mult73 > 0 && mult73 < 80) tof73[mult73 - 1] = the_tof; detn_index=30; break;
        case 74: if (mult74 > 0 && mult74 < 80) tof74[mult74 - 1] = the_tof; detn_index=31; break;
        case 81: if (mult81 > 0 && mult81 < 80) tof81[mult81 - 1] = the_tof; detn_index=32; break;
        case 82: if (mult82 > 0 && mult82 < 80) tof82[mult82 - 1] = the_tof; detn_index=33; break;
        case 83: if (mult83 > 0 && mult83 < 80) tof83[mult83 - 1] = the_tof; detn_index=34; break;
        case 84: if (mult84 > 0 && mult84 < 80) tof84[mult84 - 1] = the_tof; detn_index=35; break;
        case 91: if (mult91 > 0 && mult91 < 80) tof91[mult91 - 1] = the_tof; detn_index=36; break;
        case 92: if (mult92 > 0 && mult92 < 80) tof92[mult92 - 1] = the_tof; detn_index=37; break;
        case 93: if (mult93 > 0 && mult93 < 80) tof93[mult93 - 1] = the_tof; detn_index=38; break;
        case 94: if (mult94 > 0 && mult94 < 80) tof94[mult94 - 1] = the_tof; detn_index=39; break;
        default:
            std::cout << "ERROR in addTof: Unexpected the_detn value: " << detn_index << std::endl;
            return false;
    }
    if(detn_cathode[detn_index] < 0 || detn_cathode[detn_index] > 80) {
        std::cout << "ERROR in addTof: detn_all[" << detn_index << "] value: " << detn_cathode[detn_index] << " not valid" << std::endl;
        return false;
    }
    return true;
}
Bool_t addAmp(Int_t the_detn, Float_t the_amp) {
    Int_t detn_index;
    switch (the_detn) {
        case 01: if (mult01 > 0 && mult01 < 80) amp01[mult01 - 1] = the_amp; detn_index=0; break;
        case 02: if (mult02 > 0 && mult02 < 80) amp02[mult02 - 1] = the_amp; detn_index=1; break;
        case 03: if (mult03 > 0 && mult03 < 80) amp03[mult03 - 1] = the_amp; detn_index=2; break;
        case 04: if (mult04 > 0 && mult04 < 80) amp04[mult04 - 1] = the_amp; detn_index=3; break;
        case 11: if (mult11 > 0 && mult11 < 80) amp11[mult11 - 1] = the_amp; detn_index=4; break;
        case 12: if (mult12 > 0 && mult12 < 80) amp12[mult12 - 1] = the_amp; detn_index=5; break;
        case 13: if (mult13 > 0 && mult13 < 80) amp13[mult13 - 1] = the_amp; detn_index=6; break;
        case 14: if (mult14 > 0 && mult14 < 80) amp14[mult14 - 1] = the_amp; detn_index=7; break;
        case 21: if (mult21 > 0 && mult21 < 80) amp21[mult21 - 1] = the_amp; detn_index=8; break;
        case 22: if (mult22 > 0 && mult22 < 80) amp22[mult22 - 1] = the_amp; detn_index=9; break;
        case 23: if (mult23 > 0 && mult23 < 80) amp23[mult23 - 1] = the_amp; detn_index=10; break;
        case 24: if (mult24 > 0 && mult24 < 80) amp24[mult24 - 1] = the_amp; detn_index=11; break;
        case 31: if (mult31 > 0 && mult31 < 80) amp31[mult31 - 1] = the_amp; detn_index=12; break;
        case 32: if (mult32 > 0 && mult32 < 80) amp32[mult32 - 1] = the_amp; detn_index=13; break;
        case 33: if (mult33 > 0 && mult33 < 80) amp33[mult33 - 1] = the_amp; detn_index=14; break;
        case 34: if (mult34 > 0 && mult34 < 80) amp34[mult34 - 1] = the_amp; detn_index=15; break;
        case 41: if (mult41 > 0 && mult41 < 80) amp41[mult41 - 1] = the_amp; detn_index=16; break;
        case 42: if (mult42 > 0 && mult42 < 80) amp42[mult42 - 1] = the_amp; detn_index=17; break;
        case 43: if (mult43 > 0 && mult43 < 80) amp43[mult43 - 1] = the_amp; detn_index=18; break;
        case 44: if (mult44 > 0 && mult44 < 80) amp44[mult44 - 1] = the_amp; detn_index=19; break;
        case 51: if (mult51 > 0 && mult51 < 80) amp51[mult51 - 1] = the_amp; detn_index=20; break;
        case 52: if (mult52 > 0 && mult52 < 80) amp52[mult52 - 1] = the_amp; detn_index=21; break;
        case 53: if (mult53 > 0 && mult53 < 80) amp53[mult53 - 1] = the_amp; detn_index=22; break;
        case 54: if (mult54 > 0 && mult54 < 80) amp54[mult54 - 1] = the_amp; detn_index=23; break;
        case 61: if (mult61 > 0 && mult61 < 80) amp61[mult61 - 1] = the_amp; detn_index=24; break;
        case 62: if (mult62 > 0 && mult62 < 80) amp62[mult62 - 1] = the_amp; detn_index=25; break;
        case 63: if (mult63 > 0 && mult63 < 80) amp63[mult63 - 1] = the_amp; detn_index=26; break;
        case 64: if (mult64 > 0 && mult64 < 80) amp64[mult64 - 1] = the_amp; detn_index=27; break;
        case 71: if (mult71 > 0 && mult71 < 80) amp71[mult71 - 1] = the_amp; detn_index=28; break;
        case 72: if (mult72 > 0 && mult72 < 80) amp72[mult72 - 1] = the_amp; detn_index=29; break;
        case 73: if (mult73 > 0 && mult73 < 80) amp73[mult73 - 1] = the_amp; detn_index=30; break;
        case 74: if (mult74 > 0 && mult74 < 80) amp74[mult74 - 1] = the_amp; detn_index=31; break;
        case 81: if (mult81 > 0 && mult81 < 80) amp81[mult81 - 1] = the_amp; detn_index=32; break;
        case 82: if (mult82 > 0 && mult82 < 80) amp82[mult82 - 1] = the_amp; detn_index=33; break;
        case 83: if (mult83 > 0 && mult83 < 80) amp83[mult83 - 1] = the_amp; detn_index=34; break;
        case 84: if (mult84 > 0 && mult84 < 80) amp84[mult84 - 1] = the_amp; detn_index=35; break;
        case 91: if (mult91 > 0 && mult91 < 80) amp91[mult91 - 1] = the_amp; detn_index=36; break;
        case 92: if (mult92 > 0 && mult92 < 80) amp92[mult92 - 1] = the_amp; detn_index=37; break;
        case 93: if (mult93 > 0 && mult93 < 80) amp93[mult93 - 1] = the_amp; detn_index=38; break;
        case 94: if (mult94 > 0 && mult94 < 80) amp94[mult94 - 1] = the_amp; detn_index=39; break;
        default:
            std::cout << "ERROR in addAmp: Unexpected detn_index value: " << detn_index << std::endl;
            return false;
    }
    if (detn_cathode[detn_index] < 0 || detn_cathode[detn_index] > 80) {         
        std::cout << "ERROR in addAmp: detn_all[" << detn_index << "] value: " << detn_cathode[detn_index] << " not valid" << std::endl;
        return false;
    }
    return true;
}

// =====================================================================
// Structure for signals (common to anode/cathode)
// =====================================================================
struct Signal {
    Int_t RunNumber;
    Int_t time;
    Double_t psTime;
    Int_t BunchNumber;
    Int_t PSpulse;
    Float_t PulseIntensity;
    Int_t detn;
    Double_t tof;
    Float_t amp;
};

// =====================================================================
// Main function
// =====================================================================
void anode_cathode_final(Int_t run_number, Float_t amp_threshold, Double_t time_for_coincidence){
    // -----------------------------------------------------------------
    // Input and output files
    // -----------------------------------------------------------------
    std::string data_file = "/Users/nico/Desktop/Tese/n_TOF_data/run" + std::to_string(run_number) + ".root";
    std::string anode_file = "out_run" + std::to_string(run_number) + ".root";
    std::string out_file = "out_cathodes_" + std::to_string(run_number) + ".root";
    TFile *outfile = new TFile(out_file.c_str(), "RECREATE");
    TTree *outtree = new TTree("cathode_coincidences", "cathode_coincidences");

    TFile *fdata = TFile::Open(data_file.c_str());
    if (!fdata || fdata->IsZombie()) {
        std::cerr << "[ERROR] Cannot open data file " << data_file << std::endl;
        return;
    }

    // =====================================================================
    // Load PKUP tree (for TOF correction)
    // =====================================================================
    TTree *pkupTree = nullptr;
    fdata->GetObject("PKUP", pkupTree);
    if (!pkupTree) {
        std::cerr << "[ERROR] PKUP tree not found\n";
        fdata->Close();
        return;
    }

    Int_t pk_RunNumber, pk_time, pk_BunchNumber, pk_PSpulse;
    Double_t pk_psTime, pk_tof;
    Float_t pk_PulseIntensity;

    pkupTree->SetBranchAddress("RunNumber", &pk_RunNumber);
    pkupTree->SetBranchAddress("time", &pk_time);
    pkupTree->SetBranchAddress("psTime", &pk_psTime);
    pkupTree->SetBranchAddress("BunchNumber", &pk_BunchNumber);
    pkupTree->SetBranchAddress("PSpulse", &pk_PSpulse);
    pkupTree->SetBranchAddress("PulseIntensity", &pk_PulseIntensity);
    pkupTree->SetBranchAddress("tof", &pk_tof);

    std::map<std::tuple<Int_t, Int_t, Double_t, Int_t, Int_t, Float_t>, Double_t> pkup_map;
    Long64_t n_pk = pkupTree->GetEntriesFast();
    for (Long64_t i = 0; i < n_pk; ++i) {
        pkupTree->GetEntry(i);
        auto key = std::make_tuple(pk_RunNumber, pk_time, pk_psTime, pk_BunchNumber, pk_PSpulse, pk_PulseIntensity);
        pkup_map[key] = pk_tof;
    }
    std::cout << "[INFO] PKUP lookup map built with " << pkup_map.size() << " entries." << std::endl;

    // =====================================================================
    // Load CATHODE (PPAC) signals and apply PKUP TOF correction
    // =====================================================================
    TTree *ppacTree = nullptr;
    fdata->GetObject("PPAC", ppacTree);
    if (!ppacTree) {
        std::cerr << "[ERROR] PPAC tree not found\n";
        fdata->Close();
        return;
    }

    Signal ctmp;
    ppacTree->SetBranchAddress("RunNumber", &ctmp.RunNumber);
    ppacTree->SetBranchAddress("time", &ctmp.time);
    ppacTree->SetBranchAddress("psTime", &ctmp.psTime);
    ppacTree->SetBranchAddress("BunchNumber", &ctmp.BunchNumber);
    ppacTree->SetBranchAddress("PSpulse", &ctmp.PSpulse);
    ppacTree->SetBranchAddress("PulseIntensity", &ctmp.PulseIntensity);
    ppacTree->SetBranchAddress("detn", &ctmp.detn);
    ppacTree->SetBranchAddress("tof", &ctmp.tof);
    ppacTree->SetBranchAddress("amp", &ctmp.amp);

    std::vector<Signal> cathode_signals;
    Long64_t n_ppac = ppacTree->GetEntriesFast();
    cathode_signals.reserve(n_ppac);

    for (Long64_t i = 0; i < n_ppac; ++i) {
        ppacTree->GetEntry(i);
        auto key = std::make_tuple(ctmp.RunNumber, ctmp.time, ctmp.psTime, ctmp.BunchNumber, ctmp.PSpulse, ctmp.PulseIntensity);
        auto it = pkup_map.find(key);
        if (it == pkup_map.end()) continue;

        // TOF correction using PKUP
        ctmp.tof -= it->second;
        if (ctmp.detn < 63) ctmp.tof += 15;

        // apply amplitude and TOF cuts
        if (std::fabs(ctmp.amp) > amp_threshold && ctmp.tof > -900)
            cathode_signals.push_back(ctmp);
    }

    std::cout << "[INFO] Cathode signals loaded and corrected: " << cathode_signals.size() << std::endl;

    // Sort cathode signals by beam parameters (exclude amp, detn)
    std::sort(cathode_signals.begin(), cathode_signals.end(),
              [](const Signal &a, const Signal &b) {
                  return std::tie(a.RunNumber, a.time, a.psTime,
                                  a.BunchNumber, a.PSpulse, a.PulseIntensity, a.tof)
                         < std::tie(b.RunNumber, b.time, b.psTime,
                                    b.BunchNumber, b.PSpulse, b.PulseIntensity, b.tof);
              });

    // =====================================================================
    // Load ANODE signals (already corrected)
    // =====================================================================
    TFile *fanode = TFile::Open(anode_file.c_str(), "READ");
    if (!fanode || fanode->IsZombie()) {
        std::cerr << "[ERROR] Cannot open anode file " << anode_file << std::endl;
        fdata->Close();
        return;
    }

    TTree *anodeTree = nullptr;
    fanode->GetObject("coincidences", anodeTree);
    if (!anodeTree) {
        std::cerr << "[ERROR] Anode tree not found in " << anode_file << std::endl;
        fdata->Close();
        fanode->Close();
        return;
    }
    




    anodeTree->SetBranchAddress("RunNumber", &RunNumber);
    anodeTree->SetBranchAddress("BunchNumber", &BunchNumber);
    anodeTree->SetBranchAddress("PSpulse", &PSpulse);
    anodeTree->SetBranchAddress("PulseIntensity", &PulseIntensity);
    anodeTree->SetBranchAddress("psTime", &psTime);
    anodeTree->SetBranchAddress("time", &eventTime);
    anodeTree->SetBranchAddress("mult0", &mult0);
    anodeTree->SetBranchAddress("mult1", &mult1);
    anodeTree->SetBranchAddress("mult2", &mult2);
    anodeTree->SetBranchAddress("mult3", &mult3);
    anodeTree->SetBranchAddress("mult4", &mult4);
    anodeTree->SetBranchAddress("mult5", &mult5);
    anodeTree->SetBranchAddress("mult6", &mult6);
    anodeTree->SetBranchAddress("mult7", &mult7);
    anodeTree->SetBranchAddress("mult8", &mult8);
    anodeTree->SetBranchAddress("mult9", &mult9);
    anodeTree->SetBranchAddress("mult", &mult);
    anodeTree->SetBranchAddress("detn_all", detn_all);
    anodeTree->SetBranchAddress("tof0",tof0);
    anodeTree->SetBranchAddress("tof1", tof1);
    anodeTree->SetBranchAddress("tof2", tof2);
    anodeTree->SetBranchAddress("tof3", tof3);
    anodeTree->SetBranchAddress("tof4", tof4);
    anodeTree->SetBranchAddress("tof5", tof5);
    anodeTree->SetBranchAddress("tof6", tof6);
    anodeTree->SetBranchAddress("tof7", tof7);
    anodeTree->SetBranchAddress("tof8", tof8);
    anodeTree->SetBranchAddress("tof9", tof9);
    anodeTree->SetBranchAddress("amp0", amp0);
    anodeTree->SetBranchAddress("amp1", amp1);
    anodeTree->SetBranchAddress("amp2", amp2);
    anodeTree->SetBranchAddress("amp3", amp3);
    anodeTree->SetBranchAddress("amp4", amp4);
    anodeTree->SetBranchAddress("amp5", amp5);
    anodeTree->SetBranchAddress("amp6", amp6);
    anodeTree->SetBranchAddress("amp7", amp7);
    anodeTree->SetBranchAddress("amp8", amp8);
    anodeTree->SetBranchAddress("amp9", amp9);
    Long64_t n_anode = anodeTree->GetEntriesFast();
    std::cout << "[INFO] Anode entries: " << n_anode << std::endl;

    outtree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    outtree->Branch("BunchNumber", &BunchNumber, "BunchNumber/I");
    outtree->Branch("PSpulse", &PSpulse, "PSpulse/I");
    outtree->Branch("time", &eventTime, "time/I");
    outtree->Branch("psTime", &psTime, "psTime/D");
    outtree->Branch("PulseIntensity", &PulseIntensity, "PulseIntensity/F");
    outtree->Branch("cathode_detn", detn_cathode, "detn_cathode[40]/I");
        // multiplicities
    outtree->Branch("mult0", &new_mult0, "mult0/I");
    outtree->Branch("mult1", &new_mult1, "mult1/I");
    outtree->Branch("mult2", &new_mult2, "mult2/I");
    outtree->Branch("mult3", &new_mult3, "mult3/I");
    outtree->Branch("mult4", &new_mult4, "mult4/I");
    outtree->Branch("mult5", &new_mult5, "mult5/I");
    outtree->Branch("mult6", &new_mult6, "mult6/I");
    outtree->Branch("mult7", &new_mult7, "mult7/I");
    outtree->Branch("mult8", &new_mult8, "mult8/I");
    outtree->Branch("mult9", &new_mult9, "mult9/I");
    outtree->Branch("mult",  &new_mult,  "mult/I");

// amplitudes
    outtree->Branch("amp0", new_amp0, "amp0[mult0]/F");
    outtree->Branch("amp1", new_amp1, "amp1[mult1]/F");
    outtree->Branch("amp2", new_amp2, "amp2[mult2]/F");
    outtree->Branch("amp3", new_amp3, "amp3[mult3]/F");
    outtree->Branch("amp4", new_amp4, "amp4[mult4]/F");
    outtree->Branch("amp5", new_amp5, "amp5[mult5]/F");
    outtree->Branch("amp6", new_amp6, "amp6[mult6]/F");
    outtree->Branch("amp7", new_amp7, "amp7[mult7]/F");
    outtree->Branch("amp8", new_amp8, "amp8[mult8]/F");
    outtree->Branch("amp9", new_amp9, "amp9[mult9]/F");

// detector numbers
    outtree->Branch("detn", new_detn_all, "detn[10]/I");

// times of flight
    outtree->Branch("tof0", new_tof0, "tof0[mult0]/D");
    outtree->Branch("tof1", new_tof1, "tof1[mult1]/D");
    outtree->Branch("tof2", new_tof2, "tof2[mult2]/D");
    outtree->Branch("tof3", new_tof3, "tof3[mult3]/D");
    outtree->Branch("tof4", new_tof4, "tof4[mult4]/D");
    outtree->Branch("tof5", new_tof5, "tof5[mult5]/D");
    outtree->Branch("tof6", new_tof6, "tof6[mult6]/D");
    outtree->Branch("tof7", new_tof7, "tof7[mult7]/D");
    outtree->Branch("tof8", new_tof8, "tof8[mult8]/D");
    outtree->Branch("tof9", new_tof9, "tof9[mult9]/D");

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

    // =====================================================================
    // Coincidence search: anode vs cathode
    // =====================================================================
    for (Long64_t i = 0; i < 100; ++i) {
        anodeTree->GetEntry(i);
        for (Int_t j=0;j <10;j++){
            if (detn_all[j]>0){
                Double_t *anode_tof = nullptr; // pointer to hold tof, mult and amp values for the anode
                Float_t *anode_amp = nullptr;
                Int_t anode_mult;
                switch(j){
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
                for (Int_t k=0;k<anode_mult;k++){
                    Bool_t coincidence_found = false;
                    Double_t tmin = anode_tof[k];
                    Double_t tmax = anode_tof[k] + time_for_coincidence;
                    // Scan through cathode signals to find coincidences

                for (const auto& csignal : cathode_signals) {
                    int alpha = csignal.detn / 10; 
                    int beta  = csignal.detn % 10;
                    if (csignal.RunNumber == RunNumber && csignal.time == eventTime &&
                        std::abs(csignal.psTime - psTime) < 1e-6 && csignal.BunchNumber == BunchNumber &&
                      csignal.PSpulse == PSpulse && std::abs(csignal.PulseIntensity - PulseIntensity) < 1e-6) {
                    if (alpha == j && beta >0 && beta <5){
                        if (csignal.tof>=tmin && csignal.tof<tmax){
                            // Store cathode signal info
                            addCoincidences(csignal.detn);
                            addTof(csignal.detn, csignal.tof);
                            addAmp(csignal.detn, csignal.amp); 

                            if (!coincidence_found){
                            addCoincidences_anodes(j);
                            addTof_anodes(j, anode_tof[k]);
                            addAmp_anodes(j, anode_amp[k]);
                            new_mult++;
                            coincidence_found = true;
                            }
                        }
                        else if (csignal.tof > tmax) {
                            // Since cathode signals are sorted by TOF, we can break early
                            break;
                        }
                        else {
                            // Cathode TOF is less than tmin, continue to next cathode signal
                            continue;
                        }
                    }
                    }
                }
            }
                delete[] anode_tof;
                delete[] anode_amp;
                }
            }
        outtree->Fill();
         // Reset the variables for the next entry
        mult01 = 0, mult02 = 0, mult03 = 0, mult04 = 0, mult11 = 0, mult12 = 0, mult13 = 0, mult14 = 0, mult21 = 0, mult22 = 0, mult23 = 0, mult24 = 0, mult31 = 0, mult32 = 0, mult33 = 0, mult34 = 0, mult41 = 0, mult42 = 0, mult43 = 0, mult44 = 0, mult51 = 0, mult52 = 0, mult53 = 0, mult54 = 0, mult61 = 0, mult62 = 0, mult63 = 0, mult64 = 0, mult71 = 0, mult72 = 0, mult73 = 0, mult74 = 0, mult81 = 0, mult82 = 0, mult83 = 0, mult84 = 0, mult91 = 0, mult92 = 0, mult93 = 0, mult94 = 0;
        mult0 = 0, mult1 = 0, mult2 = 0, mult3 = 0, mult4 = 0, mult5 = 0, mult6 = 0, mult7 = 0, mult8 = 0, mult9 = 0;
        mult= 0; new_mult0 = 0; new_mult1 = 0; new_mult2 = 0; new_mult3 = 0; new_mult4 = 0; new_mult5 = 0; new_mult6 = 0; new_mult7 = 0; new_mult8 = 0; new_mult9 = 0; new_mult = 0;
        std::fill(std::begin(new_detn_all), std::end(new_detn_all), 0);
        std::fill(std::begin(detn_all), std::end(detn_all), 0);
        std::fill(std::begin(detn_cathode), std::end(detn_cathode), 0);
    }

    std::cout << "[INFO] Coincidence search complete." << std::endl;
    fdata->Close();
    fanode->Close();
    outfile->cd(); 
    outtree->Write();
    outfile->Close();
}
