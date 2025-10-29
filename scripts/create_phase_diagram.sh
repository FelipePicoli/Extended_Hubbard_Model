#!/bin/bash
# list of parameters 
code='Phase_Diagram'
model='EHM_Itensor'

exec_name=$model'_'${code}

# Gets the directory of this script
script_dir=$(dirname "$(realpath "$0")")
# Define the project root
project_root=$(realpath "$script_dir/..") 

# Use paths relative to the project root
dir_folder="$project_root/"

path_results=${dir_folder}'results/'

path_results_tmp=$path_results'tmp_'${exec_name}
path_results_final=$path_results${exec_name}

exec_file=${dir_folder}'src/code_'${exec_name}'.jl'

exec_script=${dir_folder}'src/script_'${exec_name}'.jl'

echo "Project Root: $project_root"
echo "Results Path: $path_results_tmp"
echo "Results Path: $path_results_final"

# Simulation parameters 
U0=-6.0
Uf=6.0

V0=-4.0
Vf=4.0

# DMRG parameters
Npoints=50
nsweeps=20
m=80

mkdir -p "$path_results"
mkdir -p "$path_results_final"

for L in 4 # 5 6
do
    echo "Creating phase diagram for L="${L} 
    echo "Npoints = "${Npoints}

    # Create/clear tmp folder
    # rm -rf "$path_results_tmp"
    mkdir -p "$path_results_tmp"

    # Running simulations
    #
    echo "Starting simulations"
    julia $exec_file -L $L --U0 $U0 --Uf $Uf --V0 $V0 --Vf $Vf --Npoints $Npoints --results $path_results_tmp --model $model --nsweeps $nsweeps --m $m
    echo "Simulation ended"
    echo "Setting up the data"

    # Organize the data 
    julia $exec_script  -L $L --U0 $U0 --Uf $Uf --V0 $V0 --Vf $Vf --Npoints $Npoints --results $path_results_tmp --model $model --code $code

    echo "Cleaning..."
    # Move .csv files if they exist
    if compgen -G "$path_results_tmp/*.csv" > /dev/null; then
        mv "$path_results_tmp"/*.csv "$path_results_final"
        echo "CSV files moved to: $path_results_final"
    else
        echo "WARNING: No CSV files found in $path_results_tmp. Skipping move and cleanup."
    fi
done

echo "All done."
