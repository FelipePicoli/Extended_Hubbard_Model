#!/usr/bin/python3
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description="Generate jobarrays in python with the given parameters.")
# System's parameters
parser.add_argument("-L", type=int, help="System size")
parser.add_argument("-U0", type=float, help="Initial U")
parser.add_argument("-Uf", type=float, help="Final U")
parser.add_argument("-V0", type=float, help="Initial V")
parser.add_argument("-Vf", type=float, help="Final V")

# Cluster and grid size parameters
parser.add_argument("-N_points", type=int, help="Number of points for square grid.")
parser.add_argument("-N_max_jobs", type=int, help="Max number of simultaneous jobs allowed.")

# DMRG parameters
parser.add_argument("-nsweeps", type=int, help="Number of sweeps.")
parser.add_argument("-m", type=int, help="Max linkdims.")

# Paths parameters
parser.add_argument("-results_path", type=str, help="Path where to store results.")
parser.add_argument("-project_path", type=str, help="Project path.")
parser.add_argument("-code_path", type=str, help="Name of executable.")
parser.add_argument("-code_simulation", type=str, help="Name of executable.")
parser.add_argument("-model", type=str, help="Name of simulation.")

args = parser.parse_args()

model = args.model
path_results = args.results_path
path_project = args.project_path
path_code = args.code_path
code_name = args.code_simulation

nsweeps = args.nsweeps
m = args.m

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def generate_job_arrays(L, U, V_0, batch, batch_id, batch_size, step_size,
                        model, pathresults, pathexec_dir,
                        execfile, nsweeps, m, partition, nodesjob, memory, total_time):
    """
        Write sbatch jobarray scripts for given U and V batch.
    """


    job_name = f"{model}_L={L}_U={U:.2f}_batch_size={batch_size}_batch_id={batch_id}"
    main_job_file = f"jobarray_{job_name}.sbatch"

    with open(main_job_file, "w") as f:

        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --partition={partition}\n")
        f.write(f"#SBATCH --nodes={nodesjob}\n")
        f.write(f"#SBATCH --mem={memory}\n")
        f.write(f"#SBATCH --cpus-per-task={cpustask}\n")
        f.write(f"#SBATCH --time={total_time}\n")
        f.write(f"#SBATCH --job-name={job_name}\n")
        f.write(f"#SBATCH --output=o_{job_name}.out\n")
        f.write(f"#SBATCH --error=e_{job_name}.err\n")

        # Array indices (batch defines which V indices to run)
        id_jobs_array = ",".join(str(i) for i in range(len(batch)))
        f.write(f"#SBATCH --array={id_jobs_array}\n\n")
        # ---- Environment ----
        f.write("module load julia\n\n")

        f.write('subjob_id=$(printf "%03d" $SLURM_ARRAY_TASK_ID)\n')

        v_values = [f"{val:.2f}" for val in batch]
        f.write("V_ARRAY=(" + " ".join(v_values) + ")\n")

        f.write('V0_val=${V_ARRAY[$SLURM_ARRAY_TASK_ID]}\n')
        f.write('V0_str=$(printf "%.2f" "$V0_val")\n\n')

        f.write(f'info_run="{job_name}_V0=$V0_str"\n')
        f.write('pathresults="' + str(pathresults) + '"\n\n')

        # jobfolder/jobfile names (as you had in bash)
        f.write('jobfile=job_$info_run\n')
        f.write('jobfolder=dirjob_$info_run\n')
        f.write('jobfolderDONE=dirjob_${info_run}_DONE\n')
        f.write('\n')

        # pathexec variable (folder where the original exec lives)
        f.write('pathexec="' + pathexec_dir + '"\n\n')

        # create per-job folder, copy exec, cd, run, and cleanup
        f.write('mkdir -p ${pathexec}${jobfolder}\n')
        f.write('cp "' + execfile + '" "${pathexec}${jobfolder}/"\n\n')

        f.write('cd ${pathexec}${jobfolder}\n')

        cmd = (f'julia "{execfile}" --results "{pathresults}" --model "{model}" '
               f'-L {L} --U0 {U:.6f} --V0 $V0_str --nsweeps {nsweeps} --m {m}')
        f.write(cmd + '\n\n')

        # after finishing, rename the jobfolder and remove it
        f.write('cd ..\n')
        f.write('mv ${jobfolder} ${jobfolderDONE}\n')
        f.write('rm -r ${jobfolderDONE}\n')
    print(f"[OK] Wrote {main_job_file}  (array indices: {id_jobs_array})")

N_points = args.N_points
N_max_jobs = args.N_max_jobs

L = args.L
U_0 = args.U0
U_f = args.Uf
V_0 = args.V0
V_f = args.Vf

U_values = np.linspace(U_0, U_f, N_points)

# step size needed to iterate in job array.
step_size = ((V_f - V_0)/ N_points)

partition = "SP2"
nodesjob = 1
mem = 8000
cpustask = 1
total_time = "00:30:00"

for i, U in enumerate(U_values):
    ''' 
        Create batches for jobarrays in V values.
    '''
    batches = list(chunks(np.linspace(V_0, V_f, N_points), N_max_jobs))
    for batch_id, batch in enumerate(batches):
        generate_job_arrays(L, U, V_0, batch, batch_id, len(batch), step_size, model,
                            path_results, path_code, code_name, nsweeps, m,
                            partition, nodesjob, mem, total_time)

