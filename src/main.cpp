#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Config.h"
#include "Simulation.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace simulation;

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        // print usage with specifying number of atoms
        printf("Usage: %s <number of atoms> <num_runs>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int number_of_atoms = atoi(argv[1]);
    int number_of_runs = atoi(argv[2]);
    bool use_sdl = (argc > 3 && strcmp(argv[3], "sdl") == 0);

#ifdef USE_OPENMP
    // print number of threads
    printf("Number of threads: %d\n", omp_get_num_procs());
    // Set the number of threads to the number of cores
    omp_set_num_threads(omp_get_num_procs());
#endif

#ifdef USE_MPI
    // Initialize MPI variables for ranks and size
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

#ifdef USE_SDL
    printf("Using SDL\n");
#endif

#ifdef USE_MPI
    if (rank == 0)
    {
        //number of mpi prcesses
        printf("Number of MPI processes: %d\n", size);
        printf("From rank %d: Starting simulation with the following config...\n", rank);

        printf("==============================================================\n");
        printf("\tQuadtree Config:\n");
        printf("\t   QUADTREE_MAX_DEPTH: %d\n", QUADTREE_MAX_DEPTH);
        printf("\t   QUADTREE_MAX_VALUES: %d\n", QUADTREE_MAX_VALUES);
        printf("\t   QUADTREE_THRESHOLD: %d\n", QUADTREE_THRESHOLD);
        printf("\tSimualtion configuration:\n");
        printf("\t   Number of atoms: %d\n", number_of_atoms);
        printf("\t   Number of runs: %d\n", number_of_runs);
        printf("==============================================================\n");
        if (use_sdl)
        {
            printf("\t   SDL Visualization: ENABLED\n");
        }
        else
        {
            printf("\t   SDL Visualization: DISABLED\n");
        }
        printf("==============================================================\n");
    }
#else
    printf("Starting simulation with the following config...\n");

    printf("==============================================================\n");
    printf("\tQuadtree Config:\n");
    printf("\t   QUADTREE_MAX_DEPTH: %d\n", QUADTREE_MAX_DEPTH);
    printf("\t   QUADTREE_MAX_VALUES: %d\n", QUADTREE_MAX_VALUES);
    printf("\t   QUADTREE_THRESHOLD: %d\n", QUADTREE_THRESHOLD);
    printf("\tSimualtion configuration:\n");
    printf("\t   Number of atoms: %d\n", number_of_atoms);
    printf("\t   Number of runs: %d\n", number_of_runs);
    printf("==============================================================\n");
    if (use_sdl)
    {
        printf("\t   SDL Visualization: ENABLED\n");
    }
    else
    {
        printf("\t   SDL Visualization: DISABLED\n");
    }
    printf("==============================================================\n");
#endif
    Simulation sim(number_of_atoms, number_of_runs, use_sdl);

#ifdef USE_MPI
    if (rank == 0)
    {
        sim.setup();
        sim.start_timer();
    }
#else
    sim.setup();
    sim.start_timer();
#endif

#ifdef USE_MPI
    sim.init_mpi();
    printf("From rank %d: Starting MPI simulation...\n", rank);
    sim.mpi_start(rank, size);
#else
    sim.start();
#endif

#ifdef USE_MPI
    if (rank == 0)
    {
        sim.end_timer();
        sim.print_time();
    }
    // clean up MPI
    MPI_Finalize();
#else
    sim.end_timer();

    sim.print_time();
#endif

    return EXIT_SUCCESS;
}