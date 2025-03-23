#include "SimulationMPI.h"

#include <mpi.h>

using namespace quadtree;

void simulation::SimulationMPI::init_mpi()
{
    const int num_fields = 8;
    MPI_Aint displacements[num_fields] = {
        offsetof(ParticleData, id), offsetof(ParticleData, x), offsetof(ParticleData, y),
        offsetof(ParticleData, ax), offsetof(ParticleData, ay),
        offsetof(ParticleData, vx), offsetof(ParticleData, vy),
        offsetof(ParticleData, mass)};

    int block_lengths[num_fields] = {1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype dtypes[num_fields] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_create_struct(num_fields, block_lengths, displacements, dtypes, &mpi_particle_data_type);
    MPI_Type_commit(&mpi_particle_data_type);

    std::cout << "Created MPI particle data type" << std::endl;
    std::cout << "Now calculating chunk size..." << std::endl;

    if (rank != 0)
    {
        // Compute chunk size only for non-root processes over the entire array
        int workers = size - 1; // Exclude rank 0 from the division
        int chunk_size = nb_particles / workers;
        int remainder = nb_particles % workers;

        // Assign start and end positions for each process
        startpos = (rank - 1) * chunk_size + std::min(rank - 1, remainder);
        endpos = startpos + chunk_size + ((rank - 1) < remainder ? 1 : 0);
    }

    std::cout << "MPI initialized, sending particles..." << std::endl;
}

void simulation::SimulationMPI::mpi_start()
{
    // same as start but this time to divide the computation of the mpi through multiple mpi processes
    is_running = true;

    while (is_running && /* (t < t_end) &&  */ (run_count < num_runs))
    {
        // only the compute step is parallelized throughout the mpi processes
        mpi_step(t);
        t += delta_time;
        run_count++;
    }

    // MPI_Barrier(MPI_COMM_WORLD);
}

// the step function is parallelized throughout the mpi processes
// typically in the sequential version, we create the quadtree, update the masses, and compute the acceleration of each particle based on the applied forces
// then we update the position of the particles based on the computed acceleration
// in the mpi version, we first create the tree and update the masses in the root process
// then we broadcast the 4 children of the root to the other processes
// each process will recursivly send other subtrees to the other processes
// then each process computes the acceleration of the particles in the subtree it has
// then we gather the computed accelerations from all the processes and update the position of the particles
// this way we can parallelize the computation of the forces and the position of the particles
// this is done in the mpi_step function

// the simplest parallelization method would be with a strict amount of 4 processes (1 for each child of the root)
// this way we can divide the computation of the forces and the position of the particles in 4 processes
// we will try both methods and compare the results

// for now we will implement a basic method that simply chunks the particles and distributes them to the processes
//  each process will compute the forces of the particles in the chunk it has with the entire barnes tree in memory
void simulation::SimulationMPI::mpi_step(double dtime)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // std::cout << "Rank " << rank << " is running" << std::endl;
    std::vector<ParticleData> local_particles;

    // distribute the particles to each process
    distribute_particles(local_particles, dtime);

    //if (rank != 0)
    //    std::cout << "Rank " << rank << " has to process " << endpos - startpos << " particles for " << dtime << std::endl;

    // std::cout << "Rank " << rank << " has " << local_particles.size() << " particles" << std::endl;
    // std::cout << "Rank " << rank << " now computing workload" << std::endl;

    // MPI_Barrier(MPI_COMM_WORLD);

    // save the updated particle accelerations in a buffer to send back
    std::vector<double> local_particles_acceleration /*  = std::vector<double>(nb_particles * 2) */;

    // half step
    if (rank == 0)
    {
        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i]->update_position(worldBounds, dtime * 0.5);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // compute the workload for each process
    if (rank != 0)
        compute_workload(local_particles, local_particles_acceleration, dtime);

    // MPI_Barrier(MPI_COMM_WORLD);

    // if (rank == 0)
    //     std::cout << "Rank " << rank << " now gathering particles" << std::endl;

    gather_particles(local_particles, local_particles_acceleration, dtime);
    // gather_particles(local_particles);

    // keep it simple for now and just update positions directly without integrating
    if (rank == 0)
    {
        // #pragma omp parallel for
        for (size_t i = 0; i < particles.size(); i++)
        {
            // assuming velocity is update by each other process along with acceleration
            particles[i]->update_position(worldBounds, dtime * 0.5);
        }
/*         std::cout << "Particle 1: " << particles[1]->position.x << ", " << particles[1]->position.y << std::endl;
        std::cout << "Particle 1 velocity: " << particles[1]->velocity.x << ", " << particles[1]->velocity.y << std::endl;
        std::cout << "Particle 1 acceleration: " << particles[1]->acceleration.x << ", " << particles[1]->acceleration.y << std::endl;
        std::cout << "Particle 1 mass: " << particles[1]->mass << std::endl; */
    }
}

void simulation::SimulationMPI::distribute_subtrees()
{
}

void simulation::SimulationMPI::distribute_particles(std::vector<ParticleData> &local_particles, double dtime)
{
    // serialize the particles
    // size_t all_particles = 0;

    if (rank == 0)
    {
        // all_particles = particles.size();
        //  set data for root process
        for (int idx = 0; idx < nb_particles; idx++)
        {
            local_particles.push_back(particles[idx]->serialize()); // Convert before sending
            // std::cout << "Rank " << rank << " sending Particle " << idx << " with acceleration " << local_particles[idx].ax << ", " << local_particles[idx].ay << std::endl;
        }
    }

    // Broadcast `all_particles` so all ranks know the correct number of particles
    // MPI_Bcast(&all_particles, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    local_particles.resize(nb_particles);

    // send all the particles to each process
    MPI_Bcast(local_particles.data(), nb_particles, mpi_particle_data_type, 0, MPI_COMM_WORLD);

    // send the current delta time to each process
    // MPI_Bcast(&dtime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void simulation::SimulationMPI::gather_particles(std::vector<ParticleData> &local_particles, std::vector<double> &local_accelerations, double dtime)
{
    // send each particle buffer one by one to the root process
    // this is done to avoid the need to allocate a large buffer to store all the particles
    // this way we can send the particles one by one and update the particles in the root process
    if (rank != 0)
    {
        MPI_Send(local_particles.data(), local_particles.size(), mpi_particle_data_type, 0, 0, MPI_COMM_WORLD);
    }
    else
    { // root process
        std::vector<ParticleData> global_particles;
        global_particles.resize(particles.size());
        for (int i = 0; i < size; i++)
        {
            if (i != 0)
            {
                MPI_Recv(global_particles.data(), global_particles.size(), mpi_particle_data_type, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // std::cout << "Rank " << i << " sent particles" << std::endl;
                for (int idx = 0; idx < global_particles.size(); idx++)
                {
                    particles[idx]->accumulate(global_particles[idx]);
                }
                for (int idx = 0; idx < global_particles.size(); idx++)
                {
                    // std::cout << "Rank " << i << " Particle " << idx << " at " << global_particles[idx].x << ", " << global_particles[idx].y << std::endl;
                    particles[idx]->update_velocity(dtime);
                }
            }
        }
    }
}

void simulation::SimulationMPI::compute_workload(std::vector<ParticleData> &local_particles, std::vector<double> &local_accelerations, double dtime)
{
    // this function ecompasses the work done by each process
    // this contains creating the subtree with the satrt and end positions locally and computing the forces on the particles in the subtree
    // this function will be called by each process to compute the forces on the particles in the subtree it has

    // create the tree with the particles in the subtree
    std::unique_ptr<TREE_TYPE> quadtree = std::make_unique<TREE_TYPE>(bounding_box, getBoxFunc);

    // allocate the particles as dynamic memory
    std::vector<std::unique_ptr<Particle>> all_particles = std::vector<std::unique_ptr<Particle>>(nb_particles);
    for (int idx = 0; idx < nb_particles; idx++)
    {
        all_particles[idx] = std::make_unique<Particle>();
        all_particles[idx]->deserialize(local_particles[idx]);
        // std::cout << "Rank " << rank << " Particle " << idx << " has starting acceleration " << local_particles[idx].ax << ", " << local_particles[idx].ay << std::endl;
    }

    // add the subtree particles
    for (int i = startpos; i < endpos; i++)
    {
        quadtree->add(all_particles[i].get());
    }

    // update the masses of the nodes
    quadtree->update_tree_masses();

    //if (run_count == 0)
    //    quadtree->printTree();

    // compute the forces on the particles in the subtree
    // quadtree->update_barnes_hut_forces(dtime);
    // #pragma omp parallel for
    for (int idx = 0; idx < nb_particles; idx++)
    {
        quadtree->update_individual_barnes_hut_forces(all_particles[idx].get(), dtime);
        //quadtree->update_barnes_hut_forces(dtime);
        /*         local_accelerations[idx * 2] = all_particles[idx]->acceleration.x;
                local_accelerations[idx * 2 + 1] = all_particles[idx]->acceleration.y; */
        local_particles[idx] = all_particles[idx]->serialize();
        //std::cout << "Rank " << rank << " Particle " << idx << " has after acceleration " << all_particles[idx]->acceleration.x << ", " << all_particles[idx]->acceleration.y << std::endl;
    }
}

void simulation::SimulationMPI::clean_mpi()
{
    MPI_Type_free(&mpi_particle_data_type);
}