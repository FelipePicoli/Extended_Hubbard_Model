#!/bin/sh
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
script_dir=$(dirname "$(realpath "$0")")
project_root=$(realpath "$script_dir/../")

# The path where all .txt results are stored. 
pathresults=${project_root}'/results/tmp_'${model}
execfile='code_'${model}'.py'
job_directory=${project_root}'/cluster/jobs/'

echo "Running all jobs for L=${L} in directory: $job_directory"

# Retorna o # de jobs rodando
count_current_jobs() {
    squeue -u $USER | grep -c " R\| PD"
}

mkdir -p ${job_directory}
mv *'.sbatch' ${job_directory}

# Still would be better to look for the results before running.
for sbatch_file in "$job_directory"jobarray_${model}_L=${L}_U=*.sbatch; do
    # Waits current job running to finalize. 
    while [ "$(count_current_jobs)" -ge "$MAX_JOBS" ]; do
        echo "Maximum job limit reached ($MAX_JOBS). Waiting..."
        sleep 30
        # sacct --format="user%10,jobid%10,jobname%30,state,ncpu,start,cputime,elapsed" | tail -n 50
    done
    echo "Running job: $sbatch_file"
    # sbatch "$sbatch_file"
    sleep 30
done

# Simulations ended. Running single-job to organize results 





