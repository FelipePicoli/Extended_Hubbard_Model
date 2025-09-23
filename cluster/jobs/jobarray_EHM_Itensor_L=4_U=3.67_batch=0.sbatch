#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --job-name=EHM_Itensor_L=4_U=3.67_batch=0
#SBATCH --output=o_EHM_Itensor_L=4_U=3.67_batch=0.out
#SBATCH --error=e_EHM_Itensor_L=4_U=3.67_batch=0.err
#SBATCH --array=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19

module load julia

subjob_id=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
V_ARRAY=(-6.00 -5.76 -5.51 -5.27 -5.02 -4.78 -4.53 -4.29 -4.04 -3.80 -3.55 -3.31 -3.06 -2.82 -2.57 -2.33 -2.08 -1.84 -1.59 -1.35)
V0_val=${V_ARRAY[$SLURM_ARRAY_TASK_ID]}
V0_str=$(printf "%.2f" "$V0_val")

info_run="EHM_Itensor_L=4_U=3.67_batch=0_V0=$V0_str"
pathresults="/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram"

jobfile=job_$info_run
jobfolder=dirjob_$info_run
jobfolderDONE=dirjob_${info_run}_DONE

pathexec="/home/jefter/Workspace/ehm_dmrg/src/"

mkdir -p ${pathexec}${jobfolder}
cp "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" "${pathexec}${jobfolder}/"

cd ${pathexec}${jobfolder}
julia "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" --results "/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram" --model "EHM_Itensor" --L 4 --U0 3.673469 --V0 $V0_str --nsweeps 32 --m 30

cd ..
mv ${jobfolder} ${jobfolderDONE}
rm -r ${jobfolderDONE}
