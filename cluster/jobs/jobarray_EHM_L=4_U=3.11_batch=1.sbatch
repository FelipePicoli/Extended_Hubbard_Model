#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --job-name=EHM_L=4_U=3.11_batch=1
#SBATCH --output=o_EHM_L=4_U=3.11_batch=1.out
#SBATCH --error=e_EHM_L=4_U=3.11_batch=1.err
#SBATCH --array=0,1,2

module load julia

subjob_id=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
V_ARRAY=(1.33 1.78 2.22)
V0_val=${V_ARRAY[$SLURM_ARRAY_TASK_ID]}
V0_str=$(printf "%.2f" "$V0_val")

info_run="EHM_L=4_U=3.11_batch=1_V0=$V0_str"
pathresults="../results/"

jobfile=job_$info_run
jobfolder=dirjob_$info_run
jobfolderDONE=dirjob_${info_run}_DONE

pathexec="../src/"

mkdir -p ${pathexec}${jobfolder}
cp "../src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" "${pathexec}${jobfolder}/"

cd ${pathexec}${jobfolder}
julia "../src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" --results "../results/" --model "EHM" --L 4 --U0 3.111111 --V0 $V0_str --nsweeps 6 --m 10

cd ..
mv ${jobfolder} ${jobfolderDONE}
rm -r ${jobfolderDONE}
