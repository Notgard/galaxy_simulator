#pragma once

#include "Simulation.h"

#include <mpi.h>

namespace simulation
{
    class SimulationMPI : public Simulation
    {
    public:
        SimulationMPI(int num_particles, int num_runs = -1)
            : Simulation(num_particles, num_runs, false)
        {
        }

        SimulationMPI(int num_particles, int num_runs, bool graphical)
            : Simulation(num_particles, num_runs, graphical)
        {
        }

        ~SimulationMPI()
        {

            MPI_Type_free(&mpi_particle_data_type);
        }

        void mpi_start(int rank, int size);
        void mpi_step(double dtime, int rank, int size);

        void init_mpi();

        void distribute_subtrees();

    private:
        // MPI specific variables
        MPI_Datatype mpi_particle_data_type;
    };
}