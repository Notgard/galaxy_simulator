#include "SimulationMPI.h"

#include <mpi.h>

void simulation::SimulationMPI::init_mpi()
{

    MPI_Aint displacements[8] = {
        offsetof(ParticleData, id), offsetof(ParticleData, x), offsetof(ParticleData, y),
        offsetof(ParticleData, ax), offsetof(ParticleData, ay), 
        offsetof(ParticleData, vx), offsetof(ParticleData, vy),
        offsetof(ParticleData, mass)};

    int block_lengths[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype dtypes[8] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_create_struct(8, block_lengths, displacements, dtypes, &mpi_particle_data_type);
    MPI_Type_commit(&mpi_particle_data_type);

    std::cout << "MPI initialized, sending particles..." << std::endl;

    // Broadcast particle data
}

void simulation::SimulationMPI::mpi_start(int rank, int size)
{
    // same as start but this time to divide the computation of the mpi through multiple mpi processes
    is_running = true;

    while (is_running && (t < t_end) && (run_count < num_runs))
    {
#ifdef USE_SDL
        if (rank == 0 && graphical)
        {
            ticks = SDL_GetTicks();
            perf_ticks = SDL_GetPerformanceCounter();
        }
#endif
#ifdef USE_SDL
        if (rank == 0)
            process_inputs();
#endif
        // only the compute step is parallelized throughout the mpi processes
        mpi_step(t, rank, size);
#ifdef USE_SDL
        if (rank == 0)
            render();
#endif
        t += delta_time;
        run_count++;
    }
#ifdef USE_SDL
    if (rank == 0)
        clean_sdl();
#endif
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
void simulation::SimulationMPI::mpi_step(double dtime, int rank, int size)
{

    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Rank " << rank << " is running" << std::endl;

    // each process create the tree and update the masses
    // recreate_tree();
    quadtree::Box<float> bounding_box = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
    std::unique_ptr<TREE_TYPE> quadtree = std::make_unique<TREE_TYPE>(bounding_box, getBoxFunc);

    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        std::cout << "Rank " << rank << " has added particle " << i << " at " << particles[i]->position.x << ", " << particles[i]->position.y << std::endl;
        qt->add(particles[i].get());
    }

    // update the masses of the nodes
    qt->update_tree_masses();

    // recreate_tree(quadtree, bounding_box);

    std::cout << "Rank " << rank << " has created the tree" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // Split the particle vector into chunks per process
    int chunk_size = num_particles / size;
    int remainder = num_particles % size;

    // Adjust last process's chunk to take any remainder
    int startpos = rank * chunk_size + std::min(rank, remainder);
    int endpos = (rank + 1) * chunk_size + std::min(rank + 1, remainder);

    // compute the forces for the particles in the chunk (leapfrog integration)
    for (int idx = startpos; idx < endpos; idx++)
    {
        particles[idx]->update_position(worldBounds, dtime * 0.5);
    }

    for (int idx = startpos; idx < endpos; idx++)
    {
        quadtree->update_individual_barnes_hut_forces(particles[idx].get(), dtime);
    }

    for (int idx = startpos; idx < endpos; idx++)
    {
        particles[idx]->update_velocity(dtime);
    }

    for (int idx = startpos; idx < endpos; idx++)
    {
        particles[idx]->update_position(worldBounds, dtime * 0.5);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // share the updated particles with the other processes
    std::vector<ParticleData> local_particles;
    for (int idx = startpos; idx < endpos; idx++)
    {
        local_particles.push_back(particles[idx]->serialize()); // Convert before sending
    }

    std::vector<int> recv_counts(size);
    std::vector<int> displacements(size);

    int offset = 0;
    for (int i = 0; i < size; i++)
    {
        recv_counts[i] = (num_particles / size) + (i < remainder ? 1 : 0);
        displacements[i] = offset;
        offset += recv_counts[i];
    }

    std::vector<ParticleData> all_particles(num_particles); // Buffer for all processes' data

    MPI_Allgatherv(local_particles.data(), local_particles.size(), mpi_particle_data_type,
                   all_particles.data(), recv_counts.data(), displacements.data(), mpi_particle_data_type,
                   MPI_COMM_WORLD);

    // Update each Particle with received data
    for (int i = 0; i < num_particles; i++)
    {
        if (particles[i]->id == all_particles[i].id)
            particles[i]->deserialize(all_particles[i]);
    }
}

void simulation::SimulationMPI::distribute_subtrees()
{
}