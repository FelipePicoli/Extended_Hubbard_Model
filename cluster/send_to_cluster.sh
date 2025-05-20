#!/bin/bash
# Get the directory of this script
script_dir=$(dirname "$(realpath "$0")")
# Define the project root
project_root=$(realpath "$script_dir/..") 

echo "PROJECT ROOT: "$project_root



source_code=$project_root"/src/"
scripts=$project_root"/scripts/"
cluster_scripts=$project_root"/cluster/"


scp -r $results_folder $source_code $scripts $cluster_scripts 12559016@shark.hpc.usp.br:/home/12559016/EHM_ED

