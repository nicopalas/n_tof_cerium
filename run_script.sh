#!/bin/bash
start=118558
end=118795

# Define excluded runs
is_excluded=(118578 118645 118646 118649 118650 118651 118652 118653 118654)

if [ -z "$1" ]; then
    echo "ERROR: No run_number provided!"
    exit 1
fi

run_number=$(( $1 + $start + 1 ))  # Get run_number from script arguments
# Check if run_number is in the valid range. Marcos says that this is not an optimal way since excluded runs enter in a core and exit inmediately.
# maybe try and make a macro such that for each excluded number just pass to the next...
if [[ "$run_number" -lt "$start" ]] || [[ "$run_number" -gt "$end" ]]; then
    echo "ERROR: run_number $run_number is out of range ($start - $end)!"
    exit 1
fi

# Check if run_number is excluded
for excluded in "${is_excluded[@]}"; do
    if [[ "$run_number" -eq "$excluded" ]]; then
        echo "ERROR: run_number $run_number is in the exclusion list!"
        exit 1
    fi
done

echo "Processing run_number: $run_number"

# Source the required environment
export HOME=`pwd`
source /cvmfs/lhcb.cern.ch/group_login.sh

# Execute the ROOT function. it must have the same name as the root file
lb-run davinci/latest root "/nucl_lustre/n_tof_INTC_P_665/Analysis/Macros/cathode_preselection.cpp($run_number)"


