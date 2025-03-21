#pragma once

#include "Simulation.h"

#include <mpi.h>

namespace simulation
{
    class SimulationMPI : public Simulation
    {
    public:
        SimulationMPI(int num_particles, int num_runs, int rank, int)
            : Simulation(num_particles, num_runs), rank(rank), size(size)
        {
        }

        ~SimulationMPI()
        {
        }

        void mpi_start();
        void mpi_step(double dtime);

        void init_mpi();

        void distribute_subtrees();
        void distribute_particles(std::vector<ParticleData> &local_particles, int &startpos, int &endpos, double dtime);
        void gather_particles(std::vector<Vector2<double>> &local_accelerations);

        void compute_workload(std::vector<ParticleData> &local_particles, std::vector<Vector2<double>> &local_accelerations, double dtime, int startpos, int endpos);

        void clean_mpi();

    private:
        int rank;
        int size;
        // MPI specific variables
        MPI_Datatype mpi_particle_data_type;
    };
}