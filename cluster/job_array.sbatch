#!/bin/bash
model='QuenchC'
# Get the directory of this script
script_dir=$(dirname "$(realpath "$0")")
# Define the project root
project_root=$(realpath "$script_dir/../")

path_scratch_local_env=$(realpath "$script_dir/../../")'/quspin'

# Parameters
pathresults=${project_root}'/results/'${model}'_finite_time'

execfile='code_EHM_'${model}'_finite_time.py'
pathexec=${project_root}'src/'

job_directory=${project_root}'/cluster/jobs/'

# variaveis especificas do projeto
L=$1
Ni=$2
T=$3
tau=$4

U=-5.0

V0_min=-0.5
V0_max=-10.0
dV0=-1.0 # step size
nhs=$(echo "print(int(abs((${V0_max}-${V0_min})/${dV0})))" | python3)


# variaveis especificas da rodada do cluster
mem=32000
tempo=12
nodesjob=1
partition='SP2'
cpustask=4

# identificacao arquivo de output para conferir quais parametros ainda precisam rodar
outputmain_str='W1'

# speficic job info
main_info_run=$(echo "print('%s_L=%d_U=%.1f_T=%.1f_tau=%.4f_Ni=%d' % ('${model}', $L, $U, $T, $tau, $Ni))" | python3)
main_job_file=$(echo "print('jobarray_%s_V0_min=%.1f_V0_max=%.1f_dV0=%.1f.sbatch' % ('${main_info_run}', ${V0_min}, ${V0_max}, ${dV0}))" |python3)

echo "main_info_run = " ${main_info_run}
echo "main_job_file = " ${main_job_file}

echo "# Potentials = " ${nhs}

# loop sobre parâmetros
id_jobs_array=""

# Loop through the array using indices
for (( j=0; j<$nhs; j++))
do
	V0_val=$(echo "print(${h_min} + ${j} * ${dV0})" |python3)
	info_run=$(echo "print('%s_V0_max=%.1f' % ('${main_info_run}', $V0_val))" |python3)
	info_run=$(echo "print('%s_V0_max=%.1f' % ('${main_info_run}', $V0_val))" |python3)

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
V0_min=$(echo "print(float('%.2f' % $V0_min))" |python3)
echo 'V0_val=$(echo "print(${SLURM_ARRAY_TASK_ID} *'$dV0' + '$V0_min')" |python3)' >>$main_job_file
echo 'V0_str=$(echo "print('\''%.3f'\' '% ${V0_val})" |python3)' >>$main_job_file
echo 'info_run='${main_info_run}'_V0_max=$V0_str' >> $main_job_file

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
echo "cp $pathexec$execfile $pathexec\${jobfolder}/" >> $main_job_file

echo "" >> $main_job_file
echo "cd $pathexec\${jobfolder}" >> $main_job_file
cmdexec=$(echo "print('python3 %s --dirRes \'%s\' --permutations \'%s\' --L %d -tau %f -T %f --Ni %d' % ('${pathexec}${execfile}', '${pathresults}', '${path_permutation}', $L, $tau, $T, $Ni))" | python3)
cmdexec=$cmdexec" --V0 \$V0_str --Vf \$V0_str"

echo "$cmdexec" >> $main_job_file
echo "" >> $main_job_file

# ao concluir, muda o nome da pasta
echo "cd ..">>$main_job_file
echo "mv \${jobfolder} \${jobfolderDONE}" >>$main_job_file
echo "rm -r  \${jobfolderDONE}" >>$main_job_file

echo "" >> $main_job_file
echo "" >> $main_job_file
