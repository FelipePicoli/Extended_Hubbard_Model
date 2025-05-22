#!/bin/bash

N_POINTS=50
N_JOBS_MAX=20

L=$1

U0=$2
Uf=$3
V0=$4
Vf=$5

JULIA_SCRIPT="quench_UV_grid.jl"
OUTPUT_DIR="jobs"

mkdir -p $OUTPUT_DIR

TOTAL_POINTS=$((N_POINTS * N_POINTS))

# ceiling division
POINTS_PER_JOB=$(( (TOTAL_POINTS + N_JOBS_MAX - 1) / N_JOBS_MAX ))  
# ceiling again
N_CHUNKS=$(( (TOTAL_POINTS + POINTS_PER_JOB - 1) / POINTS_PER_JOB )) 

echo "Total points: $TOTAL_POINTS"
echo "Max jobs allowed: $N_JOBS_MAX"
echo "Splitting into $N_CHUNKS chunks of up to $POINTS_PER_JOB points each"


# ========== Create jobarray files ==========
for chunk in $(seq 0 $((N_CHUNKS-1))); do
    jobarray_file="$OUTPUT_DIR/jobarray_EHM_UV_chunk_$(printf "%02d" $chunk).sbatch"

    cat <<EOF > $jobarray_file
#!/bin/bash
#SBATCH --job-name=EHM_UV_L${L}_chunk${chunk}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --array=0-$((POINTS_PER_JOB - 1))
#SBATCH --output=slurm_outputs/output_%A_%a.out
#SBATCH --error=slurm_outputs/error_%A_%a.err

# Grid and chunk metadata
N_POINTS=$N_POINTS
U_MIN=$U_MIN
U_MAX=$U_MAX
V_MIN=$V_MIN
V_MAX=$V_MAX
L=$L

# Calculate the global index in the full N_POINTS^2 grid
CHUNK_ID=$chunk
LOCAL_ID=\$SLURM_ARRAY_TASK_ID
GLOBAL_INDEX=\$((CHUNK_ID * $POINTS_PER_JOB + LOCAL_ID))

# Bounds check: skip if beyond actual grid
if [ \$GLOBAL_INDEX -ge $TOTAL_POINTS ]; then
    echo "Skipping GLOBAL_INDEX=\$GLOBAL_INDEX (outside bounds)"
    exit 0
fi

# Map global index to (i, j)
U_IDX=\$(( GLOBAL_INDEX / N_POINTS ))
V_IDX=\$(( GLOBAL_INDEX % N_POINTS ))

# Compute Uf and Vf
Uf=\$(julia -e "println($U_MIN + \$U_IDX * ($U_MAX - $U_MIN) / ($N_POINTS - 1))")
Vf=\$(julia -e "println($V_MIN + \$V_IDX * ($V_MAX - $V_MIN) / ($N_POINTS - 1))")

# Call the Julia simulation
julia $JULIA_SCRIPT --L \$L --U0 \$U0 --V0 \$V0 --Uf \$Uf --Vf \$Vf
EOF
    echo "Generated: $jobarray_file"
done

