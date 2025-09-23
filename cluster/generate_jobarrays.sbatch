#!/bin/sh
#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:00:30
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
# DMRG parameters
nsweeps=$6
m_linkdims=$7
# Simulation / cluster parameters
N_points=$8
N_max_jobs=$9

model='EHM_Itensor'
code_name='Phase_Diagram_Single_Point'

# Get the directory of this script
project_root='/home/jefter/Workspace/ehm_dmrg/' # '/temporario2/12559016/ehm_dmrg'
# The path where all .txt results are stored. 
pathresults=${project_root}'results/tmp_'${model}'_Phase_Diagram'
execfile='generate_jobarray_EHM_Itensor_Phase_Diagram.py'
job_directory=${project_root}'cluster/jobs/'
code_exec=${project_root}'src/code_'${model}_${code_name}'.jl'
code_path=${project_root}'src/'

python3 ${execfile} -L ${L} -U0 ${U0} -Uf ${Uf} -V0 ${V0} -Vf ${Vf} -nsweeps ${nsweeps} -m ${m_linkdims} -N_points ${N_points} -N_max_jobs ${N_max_jobs} -model ${model} -results_path ${pathresults} -project_path ${project_root} -code_path ${code_path} -code_simulation ${code_exec}

mkdir -p ${job_directory}
mv jobarray_*.sbatch "${job_directory}"
