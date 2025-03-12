#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include "Config.h"
#include "Simulation.h"

using namespace simulation;

int main(int argc, char **argv)
{
    if(argc < 3)
    {
        //print usage with specifying number of atoms
        printf("Usage: %s <number of atoms> <num_runs>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int number_of_atoms = atoi(argv[1]);
    int number_of_runs = atoi(argv[2]);
    bool use_sdl = (argc > 3 && strcmp(argv[3], "sdl") == 0);

    printf("Starting simulation with the following config...\n");

    printf ("==============================================================\n");
    printf ("\tQuadtree Config:\n");
    printf ("\t   QUADTREE_MAX_DEPTH: %d\n", QUADTREE_MAX_DEPTH);
    printf ("\t   QUADTREE_MAX_VALUES: %d\n", QUADTREE_MAX_VALUES);
    printf ("\t   QUADTREE_THRESHOLD: %d\n", QUADTREE_THRESHOLD);
    printf ("\tSimualtion configuration:\n");
    printf ("\t   Number of atoms: %d\n", number_of_atoms);
    printf ("\t   Number of runs: %d\n", number_of_runs);
    printf ("==============================================================\n");
    if (use_sdl)
    {
        printf("\t   SDL Visualization: ENABLED\n");
    }
    else
    {
        printf("\t   SDL Visualization: DISABLED\n");
    }
    printf("==============================================================\n");
    
    Simulation sim(number_of_atoms, number_of_runs, use_sdl);

    sim.setup();
    
    sim.start_timer();
    sim.start();
    sim.end_timer();

    sim.print_time();

    return EXIT_SUCCESS;
}