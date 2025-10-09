#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --nodes=1
#SBATCH --mem=8000
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --job-name=EHM_Itensor_L=5_U=-2.82_batch_size=10_batch_id=2
#SBATCH --output=o_EHM_Itensor_L=5_U=-2.82_batch_size=10_batch_id=2.out
#SBATCH --error=e_EHM_Itensor_L=5_U=-2.82_batch_size=10_batch_id=2.err
#SBATCH --array=0,1,2,3,4,5,6,7,8,9

module load julia

subjob_id=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
V_ARRAY=(2.53 2.69 2.86 3.02 3.18 3.35 3.51 3.67 3.84 4.00)
V0_val=${V_ARRAY[$SLURM_ARRAY_TASK_ID]}
V0_str=$(printf "%.2f" "$V0_val")

info_run="EHM_Itensor_L=5_U=-2.82_batch_size=10_batch_id=2_V0=$V0_str"
pathresults="/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram"

jobfile=job_$info_run
jobfolder=dirjob_$info_run
jobfolderDONE=dirjob_${info_run}_DONE

pathexec="/home/jefter/Workspace/ehm_dmrg/src/"

mkdir -p ${pathexec}${jobfolder}
cp "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" "${pathexec}${jobfolder}/"

cd ${pathexec}${jobfolder}
julia "/home/jefter/Workspace/ehm_dmrg/src/code_EHM_Itensor_Phase_Diagram_Single_Point.jl" --results "/home/jefter/Workspace/ehm_dmrg/results/tmp_EHM_Itensor_Phase_Diagram" --model "EHM_Itensor" -L 5 --U0 -2.816327 --V0 $V0_str --nsweeps 10 --m 10

cd ..
mv ${jobfolder} ${jobfolderDONE}
rm -r ${jobfolderDONE}
