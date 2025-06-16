#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=jefter@usp.br

# System parameters
U0=$1
Uf=$2
V0=$3
Vf=$4
Npoints=$3

# DMRG parameters
nsweeps=20
m=50

sbatch_job_name="EHM_ITensor_Phase_Diagram_L=32_U=${U}_V=${V}"
sbatch_output="o_EHM_ITensor_Phase_Diagram_L=32_U=${U}_V=${V}.sbatch"
sbatch_error="e_EHM_ITensor_Phase_Diagram_L=32_U=${U}_V=${V}.sbatch"

SBATCH --job-name=${sbatch_job_name}
SBATCH --output=${sbatch_output}
SBATCH --error=${sbatch_error}

code='Phase_Diagram'
model='EHM_Itensor'

exec_name=$model'_'${code}

# code for creating a dataframe with the scalar results
script='script_'${model}'_'${code}

# Gets the directory of this script
project_root='/temporario2/12559016/ehm_dmrg'

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

echo "Starting simulation"
module load julia # /1.11.5
# julia -e 'using InteractiveUtils; versioninfo()'

# Organize the data 
julia $exec_script  -L $L --U0 $U0 --Uf $Uf --V0 $V0 --Vf $Vf --Npoints $Npoints --tmp_results $path_results_tmp --results $path_results_final --model $model --code $code --pairs 0

