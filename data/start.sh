#!/bin/bash

run_mpi() {
    local nomeFile="$1"
    for j in {1..3}; do
        mpirun -np 1 ./exe "$nomeFile"
    done
    for i in 2 4 8 16 32 64; do
        for j in {1..10}; do
            mpirun -np $i ./exe "$nomeFile"
        done
    done
}

echo "Running email-Eu-core.txt"
run_mpi "email-Eu-core.txt"

echo "Running 1138_bus.txt"
run_mpi "1138_bus.txt"

echo "Running facebook_combined.txt"
run_mpi "facebook_combined.txt"

echo "Running Wiki-Vote.txt"
run_mpi "Wiki-Vote.txt"

echo "Running Slashdot0902.txt"
run_mpi "Slashdot0902.txt"
