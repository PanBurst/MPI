#!/bin/sh
#SBATCH -p short
#SBATCH -n64
#SBATCH -C alpha
#SBATCH -o MPI_EQUAL_RESULTS
#SBATCH -e MPI_EQUAL_ERRORS
mpic++ -pthread -o MPI_EQUAL MPI_EQUAL.cpp
mpirun MPI_EQUAL
