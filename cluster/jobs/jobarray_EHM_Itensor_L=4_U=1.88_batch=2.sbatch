#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --job-name=EHM_Itensor_L=4_U=1.88_batch=2
#SBATCH --output=o_EHM_Itensor_L=4_U=1.88_batch=2.out
#SBATCH --error=e_EHM_Itensor_L=4_U=1.88_batch=2.err
#SBATCH --array=0,1,2,3,4,5,6,7,8,9

module load julia

subjob_id=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
V_ARRAY=(3.80 4.04 4.29 4.53 4.78 5.02 5.27 5.51 5.76 6.00)
V0_val=${V_ARRAY[$SLURM_ARRAY_TASK_ID]}
V0_str=$(printf "%.2f" "$V0_val")

info_run="EHM_Itensor_L=4_U=1.88_batch=2_V0=$V0_str"
pathresults="/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram"

jobfile=job_$info_run
jobfolder=dirjob_$info_run
jobfolderDONE=dirjob_${info_run}_DONE

pathexec="/home/jefter/Workspace/ehm_dmrg/src/"

mkdir -p ${pathexec}${jobfolder}
cp "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" "${pathexec}${jobfolder}/"

cd ${pathexec}${jobfolder}
julia "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" --results "/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram" --model "EHM_Itensor" --L 4 --U0 1.877551 --V0 $V0_str --nsweeps 32 --m 30

cd ..
mv ${jobfolder} ${jobfolderDONE}
rm -r ${jobfolderDONE}
