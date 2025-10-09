#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=8000
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --job-name=EHM_Itensor_L=5_U=-4.04_batch_size=20_batch_id=1
#SBATCH --output=o_EHM_Itensor_L=5_U=-4.04_batch_size=20_batch_id=1.out
#SBATCH --error=e_EHM_Itensor_L=5_U=-4.04_batch_size=20_batch_id=1.err
#SBATCH --array=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19

module load julia

subjob_id=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
V_ARRAY=(-0.73 -0.57 -0.41 -0.24 -0.08 0.08 0.24 0.41 0.57 0.73 0.90 1.06 1.22 1.39 1.55 1.71 1.88 2.04 2.20 2.37)
V0_val=${V_ARRAY[$SLURM_ARRAY_TASK_ID]}
V0_str=$(printf "%.2f" "$V0_val")

info_run="EHM_Itensor_L=5_U=-4.04_batch_size=20_batch_id=1_V0=$V0_str"
pathresults="/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram"

jobfile=job_$info_run
jobfolder=dirjob_$info_run
jobfolderDONE=dirjob_${info_run}_DONE

pathexec="/home/jefter/Workspace/ehm_dmrg/src/"

mkdir -p ${pathexec}${jobfolder}
cp "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" "${pathexec}${jobfolder}/"

cd ${pathexec}${jobfolder}
julia "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" --results "/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram" --model "EHM_Itensor" -L 5 --U0 -4.040816 --V0 $V0_str --nsweeps 10 --m 10

cd ..
mv ${jobfolder} ${jobfolderDONE}
rm -r ${jobfolderDONE}
