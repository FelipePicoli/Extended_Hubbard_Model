#!/bin/bash
# list of parameters 
#
code='Phase_Diagram'
model='EHM_Itensor'

exec_name=$model'_'${code}

# Gets the directory of this script
script_dir=$(dirname "$(realpath "$0")")
# Define the project root
project_root=$(realpath "$script_dir/..") 

# Use paths relative to the project root
dir_folder="$project_root/"
path_results=${dir_folder}'results/'${exec_name}
exec_file=${dir_folder}'src/code_'${exec_name}'.jl'

echo "Project Root: $project_root"
echo "Results Path: $path_results"

U0=-6.0
Uf=6.0

V0=-4.0
Vf=4.0

Npoints=80
nsweeps=10
m=10

for L in 15
do
    echo "Creating phase diagram for L="${L} 
    echo "Npoints = "${Npoints}
    julia $exec_file -L $L --U0 $U0 --Uf $Uf --V0 $V0 --Vf $Vf --Npoints $Npoints --results $path_results --model $model --nsweeps $nsweeps --m $m
done
