#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00

#SBATCH --mail-type=END
#SBATCH --mail-user=jefter@usp.br

# Get U and V from command line arguments
U=$1
V=$2
Npoints=$3
L=$4

# Simulation parameters 
# DMRG parameters
nsweeps=35
m=64

# sbatch_job_name="EHM_ITensor_Phase_Diagram_L=32_U=${U}_V=${V}"
# sbatch_output="o_EHM_ITensor_Phase_Diagram_L=32_U=${U}_V=${V}.sbatch"
# sbatch_error="e_EHM_ITensor_Phase_Diagram_L=32_U=${U}_V=${V}.sbatch"

# SBATCH --job-name=${SBATCH_JOB_NAME}
# SBATCH --output=${SBATCH_OUTPUT}
# SBATCH --error=${SBATCH_ERROR}
#
code='Phase_Diagram'

code_name='GS_Pair_Points'
model='EHM_Itensor'

exec_name=$model'_'${code}

exec_name_2=$model'_'${code_name}

# code for creating a dataframe with the scalar results
script='script_'${model}'_'${code}

# Gets the directory of this script
project_root='/temporario2/12559016/ehm_dmrg'

# Use paths relative to the project root
dir_folder="$project_root/"

path_results=${dir_folder}'results/'

path_results_tmp=$path_results'tmp_'${exec_name}
path_results_final=$path_results${exec_name}

exec_file=${dir_folder}'src/code_'${exec_name_2}'.jl'

echo "Project Root: $project_root"
echo "Results Path: $path_results_tmp"
echo "Results Path: $path_results_final"
echo "Exec Path: $exec_file"

echo "Creating phase diagram for L="${L}
echo "Npoints = "${Npoints}

echo "Starting simulation"
module load julia # /1.11.5
# julia -e 'using InteractiveUtils; versioninfo()'

# Run simulations
julia $exec_file -L $L --U0 $U --V0 $V --Npoints $Npoints --results $path_results_tmp --model $model --nsweeps $nsweeps --m $m
