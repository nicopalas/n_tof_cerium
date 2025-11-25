This repository contains all the macros and files used in the analysis.
In the Macros folder, you will find two folders (one for analysis of data and one for simulation).   

  a) Analysis
    
1. anode_coincidences.cpp file searchs between the different signals in order to find anode coincidences between PPACs.
2. gamma_flash.cpp selects the gamma-flash signals as those with high multiplicities. A calibration of times between different anodes is performed.
3. neutron_energy.cpp computes the neutron energy of each reaction.
4. cathode_preselection.cpp gathers cathode signals in groups to analyze cathode coincidences in a more efficient way.
5. anode_cathode_coincidence.cpp validates all fission anode signals by imposing restricting criteria for cathode signals.   

b) Simulation. 

1. cerium_simulation.cpp contains the fission fragment generation for cerium, where results are given both for a symmetric (Z/2) split and an asymmetric splitting (driven by Z=34 shell). Kinetic energies, masses and charges for both fragments are given for each events. 
2. gef_root.cpp is a macro that loads the .lmd file from gef and converts it into a root tree where we store the neutron energy, the kinetic energies, masses and charges of the fission fragments.
3. acceptance_energy_loss.cpp is the macro in which the energy losses and the time of flight of fission fragments are computed through each stage of the setup.
4. plotting.cpp for showing the most important results. 

In addition, we also have the bash script to run the program in the shell and the condor_submit_file to run the code in paralell. All the code is written in C++ (while using ROOT). 

In the LATEX folder, you can read the results extracted from the macros above. By now:
1. acceptance_simulation contains the set of results obtained from the simulations of energy loss and time of flight of the fission fragments from uranium and cerium. 
