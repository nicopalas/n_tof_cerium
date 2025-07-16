This repository contains all the macros used in the analysis.
1. anode_analysis.cpp file searchs between the different signals in order to find anode coincidences between PPACs.
2. gamma_flash.cpp selects the gamma-flash signals as those with high multiplicities. A calibration of times between different anodes is performed.
3. neutron_energy.cpp computes the neutron energy of each reaction.
4. cathode_preselection.cpp gathers cathode signals in groups to analyze cathode coincidences in a more efficient way.
5. anode_cathode.cpp finds coincidences in cathode signals by opening a time window with respect to the anode signals. 
