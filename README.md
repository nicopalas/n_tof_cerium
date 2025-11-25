This repository contains all the macros used in the analysis.
1. anode_coincidences.cpp file searchs between the different signals in order to find anode coincidences between PPACs.
2. gamma_flash.cpp selects the gamma-flash signals as those with high multiplicities. A calibration of times between different anodes is performed.
3. neutron_energy.cpp computes the neutron energy of each reaction.
4. cathode_preselection.cpp gathers cathode signals in groups to analyze cathode coincidences in a more efficient way.
5. anode_cathode_coincidence.cpp validates all fission anode signals by imposing restricting criteria for cathode signals.

In addition, we also have the bash script to run the program in the shell and the condor_submit_file to run the code in paralell. All the code is written in C++ (while using ROOT). 
