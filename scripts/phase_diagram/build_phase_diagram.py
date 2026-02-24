#!/usr/bin/python3
import numpy as np
import subprocess
import os

model = "EHM_ITensor"
output_name = "results"

simulation_results_folder_name = "Phase_Diagram"
# Configuration
N_points = 50
Ls = [5]
output_name = "results"
# Physics Parameters
U0, Uf = -6.0, 6.0
V0, Vf = -6.0, 6.0

nsweeps=6
m=15

# Generate ranges
U_values = np.linspace(U0, Uf, N_points)
V_values = np.linspace(V0, Vf, N_points)

U_values = [f'{x:.2f}' for x in U_values]
V_values = [f'{x:.2f}' for x in V_values]

# Paths Setup
exec_name = "EHM_ITensor_Phase_Diagram_Single_Point"
script_dir =    os.path.dirname(os.path.realpath(__file__))
project_root =  os.path.realpath(os.path.join(script_dir, "../../"))
dir_folder = f"{project_root}/"
path_results = os.path.join(dir_folder, "results", simulation_results_folder_name)
path_preprocessed_operators = os.path.join(dir_folder, "preprocessing")

execfile = os.path.join(dir_folder, "src", f"code_{exec_name}.jl")

print("script dir = ", script_dir)
print("project folder = ", dir_folder)
print("path results = ", path_results)
print("path preprocessing = ", path_preprocessed_operators)
print("execfile = ", execfile)
print("output_name = ", output_name)

# Ensure results directory exists
if not os.path.exists(path_results):
    os.makedirs(path_results)
for L in Ls:
    print(f"\n### Running for L = {L} ###")
    for U in U_values:
        for V in V_values:
            # Format to 2 decimal places to match filename convention
            result_file_name = f'{output_name}_{model}_L={L}_U={U}_V={V}_NPoints={N_points}.jld2'
            info_input = f'{path_results}/{result_file_name}'

            # Check if simulation point already exists
            if os.path.isfile(info_input):
                print(f"Skipping {file_name} (already exists)")
                continue
            print(f"Running: L={L}, U={U}, V={V}")
            try:
                subprocess.run([
                    "./run_single_point.sh",
                    str(L),
                    str(U),
                    str(V),
                    str(nsweeps),
                    str(m),
                    path_results,
                    model,
                    execfile,
                    path_preprocessed_operators], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running point L={L}, U={U}: {e}")
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
