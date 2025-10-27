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
project_root='/temporario2/12559016/ehm_dmrg'

# The path where all .txt results are stored. 
pathresults=${project_root}'/results/tmp_'${model}
execfile='code_'${model}'.py'
job_directory=${project_root}'/cluster/jobs/'

# Retorna o # de jobs rodando
count_current_jobs() {
    squeue -u $USER | grep -c " R\| PD"
}

mkdir -p ${job_directory}
mv jobarray_*.sbatch "${job_directory}"

echo "Running all jobs for L=${L} in directory: $job_directory"

for sbatch_file in "$job_directory"jobarray_${model}_L=${L}_U=*.sbatch; do

    batch_size=$(echo "$sbatch_file" | sed -E 's/.*_batch_size=([0-9]+)_batch_id.*/\1/')

    echo "--------------------------------------------------------"
    echo "Checking sbatch file: $sbatch_file"

    sbatch_filename="${sbatch_file##*/}"
    err_filename="e_${sbatch_filename#jobarray_}"
    err_filename="${err_filename%.sbatch}.err"
    err_file="${project_root}/cluster/${err_filename}"

    echo "Checking error file path: $err_file"

    if [ -f "$err_file" ] && [ ! -s "$err_file" ]; then
      echo "✅ Skipping $err_file — corresponding .err file exists and is **empty** (job completed successfully)."
      continue # Skip the rest of the loop and move to the next sbatch_file.
    fi
    
    # --- 2. Failed or Incomplete Job (PROCEED to Rerun) ---
    if [ -f "$err_file" ]; then
        echo "❌ NOT EMPTY : $err_file (Job FAILED, size > 0). Rerunning this job."
        cat "$err_file"
    else
        echo "❓ MISSING FILE : $err_file (Job likely Incomplete or Failed to start). Rerunning this job."
    fi

    # --- 3. Slot Check and Submission Logic ---
    # This entire block executes ONLY if the 'continue' was NOT hit above.
    # It handles both "job failed" and "file missing" cases.
    
    # waits until entire batch can run
    while true; do
        running=$(count_current_jobs)
        free_slots=$((N_max_jobs - running))

        if [ "$free_slots" -ge "$batch_size" ]; then
            echo "Enough slots free ($free_slots ≥ $batch_size). Submitting job."
            break
        else
            echo "Not enough slots free ($free_slots < $batch_size). Waiting..."
            sleep 10
            # Assuming 'sacct' works with your cluster for job status check
            sacct --format="user%10,jobid%15,jobname%30,state,ncpu,start,cputime,elapsed" | tail -n 50
        fi
    done
    # submit job
    echo "Running job: $sbatch_file"
    sbatch "$sbatch_file"
    sleep 10
done
echo "Finished simulations."
