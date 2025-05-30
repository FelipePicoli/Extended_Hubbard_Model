#!/bin/bash
SBATCH_SCRIPT="./single_job_EHM_Itensor_Phase_Diagram.sbatch"
MAX_JOBS=20

U0=-6.0
Uf=6.0
V0=-4.0
Vf=4.0
Npoints=50

count_current_jobs() {
    # sacct --format="user%10,jobid%10,jobname%30,state,ncpu,start,cputime,elapsed" | grep "$(whoami)" | grep "RUNNING\|PENDING" | wc -l
    squeue -u "$(whoami)" -t RUNNING,PENDING | tail -n +2 | wc -l
}

echo "Generating U and V values..."

U_VALUES=$(julia -e "U0=$U0; Uf=$Uf; Npoints=$Npoints; println.(range(U0, stop=Uf, length=Npoints))" | tr '\n' ' ')
V_VALUES=$(julia -e "V0=$V0; Vf=$Vf; Npoints=$Npoints; println.(range(V0, stop=Vf, length=Npoints))" | tr '\n' ' ')

echo "U values: $U_VALUES"
echo "V values: $V_VALUES"

for U_val in $U_VALUES; do
    for V_val in $V_VALUES; do
        while [ "$(count_current_jobs)" -ge "$MAX_JOBS" ]; do
            echo "Maximum job limit reached ($MAX_JOBS). Waiting..."
            sleep 60
            sacct --format="user%10,jobid%10,jobname%30,state,ncpu,start,cputime,elapsed" | tail -n 10 
        done

        # Construct job-specific names for --job-name, --output, --error
        JOB_NAME="EHM_L32_U${U_val}_V${V_val}"
        OUTPUT_FILE="o_L32_U${U_val}_V${V_val}.sbatch"
        ERROR_FILE="e_L32_U${U_val}_V${V_val}.sbatch"

        echo "Submitting job for U=${U_val}, V=${V_val}"
        sbatch "${SBATCH_SCRIPT}" "${U_val}" "${V_val}" "${Npoints}"
        sleep 50 
    done
done
echo "All jobs submitted."
