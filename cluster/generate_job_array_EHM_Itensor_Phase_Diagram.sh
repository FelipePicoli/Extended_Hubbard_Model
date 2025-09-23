#!/bin/sh
#       This generates a list of jobarrays for computing the entire phase diagram of the EHM.
#       Input parameters of interest are : U_min, U_max, V_min, V_max, N_points, N_max_jobs 
#
code='Phase_Diagram'
# Get the directory of this script
script_dir=$(dirname "$(realpath "$0")")
# Define the project root
project_root=$(realpath "$script_dir/../")

path_scratch_local_env=$(realpath "$script_dir/../../")'/quspin'

# Paths 
pathresults=${project_root}'/results/EHM_Itensor_'${code}
execfile='code_EHM_Itensor_'${code}'_Single_Point.py'
pathexec=${project_root}'src/'
job_directory=${project_root}'/cluster/jobs/'

# System's variables 
L=$1
U_min=$2
U_max=$3
V_min=$4
V_max=$5

# DMRG parameters
Npoints=$6
nsweeps=$7
m=$8














