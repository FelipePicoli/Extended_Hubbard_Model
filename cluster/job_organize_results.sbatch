#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jefter@usp.br
#SBATCH --partition=SP2

#SBATCH --job-name=organize_results_L=4_batch
#SBATCH --output=o_organize_results_L=4_batch.out
#SBATCH --error=e_organize_results_L=4_batch.err

# System parameters
L=$1
U0=$2
Uf=$3
V0=$4
Vf=$5
Npoints=$6

model='EHM_Itensor'
code_name='Phase_Diagram_Single_Point'
script_name='Phase_Diagram'

exec_name=$model'_'${code}
# code for creating a dataframe with the scalar results
script='script_'${model}'_'${script_name}

# Gets the directory of this script
project_root=$7 #  '/temporario2/12559016/ehm_dmrg'
# Use paths relative to the project root
dir_folder="$project_root/"
path_results=${dir_folder}'results/'
path_results_tmp=$path_results'tmp_'${exec_name}
path_results_final=$path_results${exec_name}

exec_script=${dir_folder}'src/'${script}'.jl'

echo "Project Root: $project_root"
echo "Results Path: $path_results_tmp"
echo "Results Path: $path_results_final"
echo "Script Path: $exec_script"
echo "Creating phase diagram for L="${L} 
echo "Npoints = "${Npoints}

echo "Organizing results from simulation into .csv files."
module load julia
# Organize the data 
julia $exec_script  -L $L --U0 $U0 --Uf $Uf --V0 $V0 --Vf $Vf --Npoints $Npoints --tmp_results $path_results_tmp --results $path_results_final --model $model --code $code_name 
