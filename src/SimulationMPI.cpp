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
}

void simulation::SimulationMPI::mpi_start()
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
        mpi_step(t);
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
void simulation::SimulationMPI::mpi_step(double dtime)
{

    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Rank " << rank << " is running" << std::endl;
    std::vector<ParticleData> local_particles;
    int startpos, endpos;

    // distribute the particles to each process
    distribute_particles(local_particles, startpos, endpos, dtime);

    std::cout << "Rank " << rank << " has " << local_particles.size() << " particles" << std::endl;
    std::cout << "Rank " << rank << " now computing workload" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // compute the workload for each process
    compute_workload(local_particles, dtime, startpos, endpos);

    if (rank == 0)
        std::cout << "Rank " << rank << " now gathering particles" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    gather_particles(local_particles);
    /*
        // keep it simple for now and just update positions directly without integrating
        if (rank == 0)
        {
    #pragma omp parallel for
            for (size_t i = 0; i < particles.size(); i++)
            {
                // assuming velocity is update by each other process along with acceleration
                particles[i]->update_position(worldBounds, dtime);
            }
        } */
}

void simulation::SimulationMPI::distribute_subtrees()
{
}

void simulation::SimulationMPI::distribute_particles(std::vector<ParticleData> &local_particles, int &startpos, int &endpos, double dtime)
{
    // serialize the particles
    size_t all_particles = 0;

    if (rank == 0)
    {
        all_particles = particles.size();
        // set data for root process
        for (int idx = 0; idx < all_particles; idx++)
        {
            local_particles.push_back(particles[idx]->serialize()); // Convert before sending
        }
    }

    // Broadcast `all_particles` so all ranks know the correct number of particles
    MPI_Bcast(&all_particles, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    local_particles.resize(all_particles);

    // send all the particles to each process
    MPI_Bcast(local_particles.data(), all_particles, mpi_particle_data_type, 0, MPI_COMM_WORLD);

    // calculate chunk of particles for each process so that each process has acess to all particles and can compute the forces on a local subtree
    // the given chunk will be the number of particles in the subtree for each process
    int chunk_size = all_particles / size;
    int remainder = all_particles % size;

    // Adjust last process's chunk to take any remainder
    startpos = rank * chunk_size + std::min(rank, remainder);
    endpos = (rank + 1) * chunk_size + std::min(rank + 1, remainder);

    // send the current delta time to each process
    MPI_Bcast(&dtime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void simulation::SimulationMPI::gather_particles(std::vector<ParticleData> &local_particles)
{
    // gather the particles from all the processes from the root process
    std::vector<ParticleData> global_particles;
    if (rank == 0)
    {
        std::cout << "New size: " << particles.size() << std::endl;
        std::cout << "New size: " << local_particles.size() << std::endl;
        std::cout << "New size: " << global_particles.size() << std::endl;
        global_particles.resize(particles.size());
    }

    // gather the particles from all the processes
    MPI_Gather(local_particles.data(), local_particles.size(), mpi_particle_data_type,
               global_particles.data(), global_particles.size(), mpi_particle_data_type, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        // Update each Particle with received data
        for (int i = 0; i < particles.size(); i++)
        {
            if (particles[i]->id == global_particles[i].id)
                particles[i]->deserialize(global_particles[i]);
        }
    }
}

void simulation::SimulationMPI::compute_workload(std::vector<ParticleData> &local_particles, double dtime, int startpos, int endpos)
{
    // this function ecompasses the work done by each process
    // this contains creating the subtree with the satrt and end positions locally and computing the forces on the particles in the subtree
    // this function will be called by each process to compute the forces on the particles in the subtree it has

    // create the tree with the particles in the subtree
    quadtree::Box<float> bounding_box = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
    std::unique_ptr<TREE_TYPE> quadtree = std::make_unique<TREE_TYPE>(bounding_box, getBoxFunc);

    size_t buffer_size = local_particles.size();

    // allocate the particles as dynamic memory
    std::vector<std::unique_ptr<Particle>> all_particles = std::vector<std::unique_ptr<Particle>>(buffer_size);
    for (int idx = 0; idx < buffer_size; idx++)
    {
        all_particles[idx] = std::make_unique<Particle>();
        all_particles[idx]->deserialize(local_particles[idx]);
    }

    // add the subtree particles
    for (int i = startpos; i < endpos; i++)
    {
        quadtree->add(all_particles[i].get());
    }

    // update the masses of the nodes
    quadtree->update_tree_masses();

    // compute the forces on the particles in the subtree
    quadtree->update_barnes_hut_forces(dtime);

// update the velocity of the particles
#pragma omp parallel for
    for (int idx = 0; idx < all_particles.size(); idx++)
    {
        all_particles[idx]->update_velocity(dtime);
    }

    // serialize the particles to send back
    for (int idx = 0; idx < buffer_size; idx++)
    {
        local_particles[idx] = all_particles[idx]->serialize();
    }
}

void simulation::SimulationMPI::clean_mpi()
{
    MPI_Type_free(&mpi_particle_data_type);
}