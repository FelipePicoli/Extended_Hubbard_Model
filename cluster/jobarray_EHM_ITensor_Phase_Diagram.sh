#!/bin/bash

# Parameters
#
U_0=0
V_0=0

L=$1
U_f=$2
V_f=$3

model='EHM_'${quench}


# Get the directory of this script
script_dir=$(dirname "$(realpath "$0")")
# Define the project root
project_root=$(realpath "$script_dir/../")

path_scratch_local_env=$(realpath "$script_dir/../../")'/quspin'

# Files and directories#
pathresults=${project_root}'/results/'${model}
job_directory=${project_root}'/cluster/jobs/'

pathexec=${project_root}'/src/'
execfile=$pathexec'code_EHM_Single_Quench_finite_time.py'

echo $model
echo $script_dir
echo $project_root
echo $path_scratch_local
echo $path_scratch_local_env
echo $pathexec
echo $execfile
echo $pathresults

# mkdir $pathresults


# discretization of time.
nhs=$(echo "print(int(abs((${tau_max}-${tau_min})/${dt})))" | python3)
# variáveis específicas da rodada do cluster

mem=32000
tempo=12
nodesjob=1
partition='SP2'
cpustask=4

# identificacao arquivo de output para conferir quais parametros ainda precisam rodar
outputmain_str='magnetization'

# speficic job info
main_info_run=$(echo "print('%s_L=%d_U0=%.1f_Uf=%.1f_V0=%.1f_Vf=%.1f' % ('${model}', $L, $U_0, $U_f, $V_0, $V_f))" | python3)
main_job_file=$(echo "print('jobarray_%s_tau_min=%.3f_tau_max=%.1f_dt=%.1f.sbatch' % ('${main_info_run}', ${tau_min}, ${tau_max}, ${dt}))" |python3)

echo "main_info_run = " ${main_info_run}
echo "main_job_file = " ${main_job_file}

echo "# Potentials = " ${nhs}

# loop sobre parâmetros
id_jobs_array=""
for (( j=0; j<$nhs; j++))
do
	tau_val=$(echo "print(${tau_min} + ${j} * ${dt})" |python3)
	info_run=$(echo "print('%s_tau_max=%.1f' % ('${main_info_run}', $tau_val))" |python3)

	outputfile_main=${outputmain_str}'_'${info_run}'.txt'
	echo $outputfile_main
	# checa se o output existe, caso contrário coloca na lista de parâmetros
	if [ ! -f $pathresults"/"$outputfile_main ]; then
		strid="${j},"
		if [ $(bc <<< "$j == $nhs") -eq 1 ]
		then
			strid="${j}"
		fi
			id_jobs_array+=$strid

		else
			echo "already calculated!"
	fi
done

echo "id jobs array = " $id_jobs_array

# começa a escrever no .sbatch que irá ser enviado para o sistema de filas
echo "#!/bin/bash" > $main_job_file
echo "# "$main_job_file >> $main_job_file

echo "#SBATCH --partition="$partition >> $main_job_file
echo "#SBATCH --nodes="$nodesjob >> $main_job_file
echo "#SBATCH --mem="$mem >> $main_job_file
echo "#SBATCH --cpus-per-task="$cpustask >> $main_job_file
echo "#SBATCH --time="$tempo":00:00" >> $main_job_file
echo "" >> $main_job_file
echo "#SBATCH --job-name="$main_job_file >> $main_job_file
echo "#SBATCH --output=o_"$main_job_file >> $main_job_file
echo "#SBATCH --error=e_"$main_job_file >> $main_job_file
echo "#SBATCH --array="$id_jobs_array >>$main_job_file

echo "" >> $main_job_file
echo "module load Anaconda" >> $main_job_file
echo "source activate ${path_scratch_local_env}" >> $main_job_file

echo "" >> $main_job_file
# cada job id é associado a um parâmetro
echo 'subjob_id='\`'printf %03d $SLURM_ARRAY_TASK_ID'\`>> $main_job_file
echo "" >> $main_job_file

# cada variável auxiliar necessária para executar um parâmetro é enviada ao jobfile
V0_min=$(echo "print(float('%.2f' % $tau_min))" |python3)
echo 'tau_val=$(echo "print(${SLURM_ARRAY_TASK_ID} *'$dt' + '$tau_min')" |python3)' >> $main_job_file
echo 'tau_str=$(echo "print('\''%.3f'\' '% ${tau_val})" |python3)' >> $main_job_file
echo 'info_run='${main_info_run}'_tau_max=$tau_str' >> $main_job_file

echo 'pathresults='$pathresults >> $main_job_file
echo 'outputfile_main='${outputmain_str}'_${info_run}.txt' >> $main_job_file

# Write the check for output file existence to the job file
echo '# Check if the output file already exists' >> $main_job_file
echo 'if [ -f ${pathresults}/${outputfile_main} ]; then' >> $main_job_file
echo '    echo "Output file ${outputfile_main} already exists. Skipping job."' >> $main_job_file
echo '    exit 0' >> $main_job_file
echo 'else' >> $main_job_file
echo '    echo "Output file ${outputfile_main} does not exist. Running the job."' >> $main_job_file
echo 'fi' >> $main_job_file

echo 'jobfile=job_$info_run' >> $main_job_file
echo 'jobfolder=dirjob_$info_run' >> $main_job_file
echo 'jobfolderDONE=dirjob_${info_run}_DONE' >> $main_job_file
#
echo "" >> $main_job_file

# repete parte da sequência de comandos do .sbath para um parâmetro único

echo "mkdir $pathexec\${jobfolder}" >> $main_job_file
echo "cp $execfile $pathexec\${jobfolder}/" >> $main_job_file

echo "" >> $main_job_file
echo "cd $pathexec\${jobfolder}" >> $main_job_file

cmdexec=$(echo "print('julia %s -L %d -U0 %f -Uf %f -V0 %f -Vf %f --results_dir %s --quench %s --model %s' % ('${execfile}', $L, $U_0, $U_f, $V_0, $V_f, '${pathresults}', '${quench}', '${model}'))" | python3)
cmdexec=$cmdexec" --tau \$tau_str"

echo "$cmdexec" >> $main_job_file
echo "" >> $main_job_file

# ao concluir, muda o nome da pasta
echo "cd ..">>$main_job_file
echo "mv \${jobfolder} \${jobfolderDONE}" >>$main_job_file
echo "rm -r  \${jobfolderDONE}" >>$main_job_file

echo "" >> $main_job_file
echo "" >> $main_job_file
