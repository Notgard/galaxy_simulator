#pragma once

#include "SimulationMPI.h"
#include "SimulationSDL.h"

namespace simulation
{
    class SimulationMPISDL : public SimulationMPI, public SimulationSDL
    {
    public:
        SimulationMPISDL(int num_particles, int num_runs,int mpi_version, int rank, int size, bool graphical = false)
            : Simulation(num_particles, num_runs), 
            SimulationMPI(num_particles, num_runs,mpi_version, rank, size),  
            SimulationSDL(num_particles, num_runs, graphical)
        {
        }

        ~SimulationMPISDL()
        {
        }

        void mpi_start() override;
    };
}