#pragma once

#include "Simulation.h"

#include <mpi.h>

namespace simulation
{
    class SimulationMPI : public virtual Simulation
    {
    public:
        SimulationMPI(int num_particles, int num_runs, int rank, int size)
            : Simulation(num_particles, num_runs), rank(rank), size(size)
        {
            bounding_box = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
            nb_particles = num_particles + CELESTIAL_BODY_COUNT;
        }

        ~SimulationMPI()
        {
        }

        virtual void mpi_start();
        void mpi_step(double dtime);

        void init_mpi();

        void distribute_subtrees();
        void distribute_particles(std::vector<ParticleData> &local_particles, double dtime);
        void gather_particles(std::vector<ParticleData> &local_particles, double dtime);
        void gather_particles(std::vector<double> &local_accelerations, double dtime);
        //void gather_particles(std::vector<ParticleData> &local_particles, std::vector<Vector2<double>> &local_accelerations, double dtime);

        void compute_workload(std::vector<ParticleData> &local_particles, std::vector<double> &local_accelerations, double dtime);
        //void compute_workload(std::vector<ParticleData> &local_particles, std::vector<Vector2<double>> &local_accelerations, double dtime, int startpos, int endpos);

        void clean_mpi();

    protected:
        int rank;
        int size;

        size_t nb_particles;

        int startpos, endpos;

        quadtree::Box<float> bounding_box;

        // MPI specific variables
        MPI_Datatype mpi_particle_data_type;
        MPI_Datatype mpi_particle_type;
    };
}