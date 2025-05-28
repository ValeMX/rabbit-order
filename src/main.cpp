#include <mpi.h>
#include <omp.h>

#include <ctime>
#include <sstream>
#include <string>

#include "graph.h"

using namespace std;

char* fileName = NULL;

void parseArguments(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-f" && i + 1 < argc) {
            fileName = argv[++i];
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize OpenMP with the number of threads equal to the number of MPI processes
    omp_set_num_threads(size);

    // Start timing
    double start_time = MPI_Wtime();

    string s(fileName);
    s += "_" + to_string(rank) + ".txt";

    Graph g(s);

    // End timing
    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Total execution time: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}