#!/bin/bash

run_mpi() {
    local file="$1"
    for j in {1..3}; do
        mpirun -np 1 ./exe "$file"
    done
    for i in 2 4 8 16 32 64; do
        for j in {1..5}; do
            mpirun -np $i ./exe "$file"
        done
    done
}

# Compilation
make clean
make

# Create results folder if it does not exist
mkdir -p results
cd results

folder_name=$(date +"%Y-%m-%d_%H-%M-%S")
mkdir "$folder_name"
cd "$folder_name"

for file in ../../data/*.txt; do
    echo "Running $(basename "$file")"
    run_mpi "$file"
done