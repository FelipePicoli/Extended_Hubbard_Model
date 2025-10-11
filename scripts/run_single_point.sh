#!/bin/bash
# list of parameters 
code='Phase_Diagram_Single_Point'
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

echo "Project Root: $project_root"
echo "Results Path: $path_results_tmp"
echo "Results Path: $path_results_final"


# Simulation parameters 
L=$1
U=$2
V=$3

# DMRG parameters
nsweeps=6
m=10

mkdir -p "$path_results"
mkdir -p "$path_results_final"


echo "Creating phase diagram for L="${L} 

mkdir -p "$path_results_tmp"

echo "Starting simulations"
julia $exec_file -L $L --U0 $U --V0 $V --results $path_results_tmp --model $model --nsweeps $nsweeps --m $m

# Organize the data 
# julia $exec_script  -L $L --U0 $U0 --Uf $Uf --V0 $V0 --Vf $Vf --Npoints $Npoints --results $path_results_tmp --model $model --code $code
echo "All done."
