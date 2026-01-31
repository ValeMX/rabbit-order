# Distributed Rabbit Order

**Distributed Rabbit Order** is a distributed parallel scalable graph reordering algorithm for High-Performance Computing. The objective of this project is to implement the Rabbit Order graph reordering algorithm for distributed memory systems using MPI.

## Repository Structure

```

rabbit-order/
├── data/                           # Input datasets 
├── results/                        # Output results and graphs
├── scripts/                        # Results visualization scripts (Python)
├── src/                            # Source code (C++)
├── Distributed Rabbit Order.pdf    # Thesis
├── run.sh                          # Execution without compilation
├── start.sh                        # Compilation and execution
├── Makefile                        # Build configuration
└── README.md                       # Project overview and instructions

```

---

## Prerequisites

Before you begin, make sure you have the following installed:

- **GCC-9.1.0**
- **MPICH-3.2.1** 
- **Python 3.8+** (only for results visualization)
- **Make**

## Build Instructions

### 1. Clone the Repository

```bash
git clone https://github.com/ValeMX/rabbit-order.git
cd rabbit-order
```

### 2. Build

```bash
make
```

## Execution Instructions

### 1. Manual Execution

You can manually execute the program by typing the following:

```bash
mpirun -np <number_of_processes> <executable> <input_file>
```

Results will be generated in an ```out.csv``` file in the same folder.

### 2. Batch Execution

Assuming that the executable is ```/bin/exe``` file, you can execute the program with multiple inputs by typing the following:

```bash
chmod +x run.sh
./run.sh file1.txt [file2.txt ...]
```
This will run the algorithm on every input file in the following way: 
- 3 times with 1 process 
- 5 times with 2, 4, 8, 16, 32, 64 processes each

Results will be saved in a folder with the current datetime in the ```results``` directory.

### 3. Build and run on whole dataset

In order to execute all the tests, it is possible to type the following:

```bash
chmod +x start.sh
./start.sh
```

This will build and run the algorithm on all the input ```*.txt``` files in the ```data``` directory. Results will be saved in a subfolder in the ```results``` directory.

## Results visualization

In the ```scripts``` folder there are some scripts to visualize results. Assuming python is provided with:

- pandas
- numpy
- matplotlib

It is possible to execute the following scripts by running them and inserting the ```.csv``` results file:

- **bars.py** shows the execution times of different computation steps for each number of processes
- **speedup.py** computes and shows the strong scaling and efficiency for each input graph and for each number of processess
- **times.py** shows the total and specific (single computation step) execution times of all the input graphs by the number of their edges and their nodes for each number of processes