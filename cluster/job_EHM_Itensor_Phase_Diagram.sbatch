#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --time=10:30:00

#SBATCH --job-name=job_EHM_Itensor_Phase_Diagram.sbatch
#SBATCH --output=o_job_EHM_Itensor_Phase_Diagram.sbatch
#SBATCH --error=e_job_EHM_Itensor_Phase_Diagram.sbatch

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jefter@usp.br

# Parameters
L=32

U0=-8.0
Uf=8.0
V0=-6.0
Vf=6.0

N_points=50

model='EHM_Itensor'

# Get the directory of this script
script_dir=$(dirname "$(realpath "$0")")
# Define the project root
project_root=$(realpath "$script_dir/../")

# Files and directories#
path_results=${project_root}'/results/'${model}
job_directory=${project_root}'/cluster/jobs/'

pathexec=${project_root}'/src/'
julia_env=${project_root}'/julia_env/'

exec_file=$pathexec'code_'$model'_Phase_Diagram.jl'

echo $script_dir 
echo $project_root
echo $path_results
echo $job_directory 
echo $pathexec
echo $julia_env 
echo $exec_file


module load julia

# path to the julia_env#
export JULIA_PROJECT=$julia_env
# (Optional sanity check)
julia -e 'using Pkg; println("Active env: ", Pkg.project().path)'

# Run the Julia script
julia $exec_file -L $L --U0 $U0 --Uf $Uf --V0 $V0 --Vf $Vf --Npoints $Npoints --results $path_results --model $model 
