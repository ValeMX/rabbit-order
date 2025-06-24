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

# Create results folder if it does not exist
mkdir -p results
cd results

folder_name=$(date +"%Y-%m-%d_%H-%M-%S")
mkdir "$folder_name"
cd "$folder_name"

cp ../../bin/exe .

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 file1.txt [file2.txt ...]"
    exit 1
fi

for file in "$@"; do
    full_path="../../data/$file"
    echo "Running $(basename "$full_path")"
    run_mpi "$full_path"
done
