#!/bin/bash
# Simulation parameters
L=$1
U=$2
V=$3

# DMRG parameters
nsweeps=$4
m=$5

# Paths parameters
path_results=$6
model=$7
exec_file=$8

mkdir -p "$path_results"

julia $exec_file --L $L --U0 $U --V0 $V --results $path_results --model $model --nsweeps $nsweeps --m $m
