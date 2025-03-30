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

        std::cout << "Rank " << rank << " processing " << endpos - startpos << " particles" << std::endl;
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

    // process particles depending on the method
    if (mpi_version == SUB_TREE_PROCESSING)
    {
        sub_tree_processing(local_particles, dtime);
    }
    else if (mpi_version == CHUNKED_PROCESSING)
    {
        chunked_processing(local_particles, dtime);
    }
}

void simulation::SimulationMPI::distribute_subtrees()
{
}

void simulation::SimulationMPI::distribute_particles(std::vector<ParticleData> &local_particles, double dtime)
{
    // serialize the particles
    if (rank == 0)
    {
// all_particles = particles.size();
//  set data for root process
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (size_t idx = 0; idx < nb_particles; idx++)
        {
            local_particles.push_back(particles[idx]->serialize()); // Convert before sending
        }
    }

    // Broadcast `all_particles` so all ranks know the correct number of particles
    local_particles.resize(nb_particles);

    // send all the particles to each process
    MPI_Bcast(local_particles.data(), nb_particles, mpi_particle_data_type, 0, MPI_COMM_WORLD);
}

void simulation::SimulationMPI::gather_particles(std::vector<ParticleData> &local_particles, double dtime)
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
                for (size_t idx = 0; idx < global_particles.size(); idx++)
                {
                    particles[idx]->accumulate(global_particles[idx]);
                }
                for (size_t idx = 0; idx < global_particles.size(); idx++)
                {
                    particles[idx]->update_velocity(dtime);
                }
            }
        }
    }
}

void simulation::SimulationMPI::gather_particles(std::vector<double> &local_accelerations, double dtime)
{
    // Reduce local accelerations to sum onto process 0
    std::vector<double> global_accelerations;
    if (rank == 0)
    {
        global_accelerations.resize(local_accelerations.size(), 0.0);
    }

    MPI_Reduce(local_accelerations.data(), global_accelerations.data(), local_accelerations.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
// Apply accumulated accelerations to particles
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i]->acceleration.x = global_accelerations[2 * i];
            particles[i]->acceleration.y = global_accelerations[2 * i + 1];
            particles[i]->update_velocity(dtime);
            particles[i]->update_position(worldBounds, dtime);
        }
    }
}

void simulation::SimulationMPI::gather_chunk_particles(std::vector<ParticleData> &updated_chunk, double dtime)
{
    int local_size = updated_chunk.size();
    std::vector<int> recv_counts(size);
    std::vector<int> displs(size);

    // Root process collects the sizes of each chunk
    MPI_Gather(&local_size, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<ParticleData> gathered_particles; // Only used on rank 0

    if (rank == 0)
    {
        int total_particles = 0;
        displs[0] = 0;

        for (int i = 0; i < size; i++)
        {
            displs[i] = total_particles;
            total_particles += recv_counts[i];
        }

        gathered_particles.resize(total_particles);
    }

    // Gather all particle data back to rank 0
    MPI_Gatherv(updated_chunk.data(), local_size, mpi_particle_data_type,
                gathered_particles.data(), recv_counts.data(), displs.data(),
                mpi_particle_data_type, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < gathered_particles.size(); i++)
        {
            particles[gathered_particles[i].id]->update(gathered_particles[i]);
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
    //std::vector<std::unique_ptr<Particle>> all_particles = std::vector<std::unique_ptr<Particle>>(nb_particles);
    std::vector<Particle> all_particles = std::vector<Particle>(nb_particles);
    for (size_t idx = 0; idx < nb_particles; idx++)
    {
        all_particles[idx].id = local_particles[idx].id;
        all_particles[idx].position = {local_particles[idx].x, local_particles[idx].y};
        all_particles[idx].acceleration = {local_particles[idx].ax, local_particles[idx].ay};
        all_particles[idx].velocity = {local_particles[idx].vx, local_particles[idx].vy};
        all_particles[idx].mass = local_particles[idx].mass;
        //all_particles[idx] = std::make_unique<Particle>();
        //all_particles[idx]->deserialize(local_particles[idx]);
    }

    // add the subtree particles
    for (int i = startpos; i < endpos; i++)
    {
        //quadtree->add(all_particles[i].get());
        quadtree->add(&all_particles[i]);
    }

    // update the masses of the nodes
    quadtree->update_tree_masses();

// if (run_count == 0)
//     quadtree->printTree();

// compute the forces on the particles in the subtree
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t idx = 0; idx < nb_particles; idx++)
    {
        quadtree->update_individual_barnes_hut_forces(&all_particles[idx], dtime);
        local_accelerations[idx * 2] = all_particles[idx].acceleration.x;
        local_accelerations[idx * 2 + 1] = all_particles[idx].acceleration.y;
        //quadtree->update_individual_barnes_hut_forces(all_particles[idx].get(), dtime);
        //local_accelerations[idx * 2] = all_particles[idx]->acceleration.x;
        //local_accelerations[idx * 2 + 1] = all_particles[idx]->acceleration.y;
        // local_particles[idx] = all_particles[idx]->serialize();
    }
}

void simulation::SimulationMPI::compute_entire_workload(std::vector<ParticleData> &local_particles, std::vector<ParticleData> &computed_chunk, double dtime)
{
    // this function ecompasses the work done by each process
    // this contains creating the subtree with the satrt and end positions locally and computing the forces on the particles in the subtree
    // this function will be called by each process to compute the forces on the particles in the subtree it has

    // create the tree with the particles in the subtree
    std::unique_ptr<TREE_TYPE> quadtree = std::make_unique<TREE_TYPE>(bounding_box, getBoxFunc);

    // allocate the particles as dynamic memory
    //std::vector<std::unique_ptr<Particle>> all_particles = std::vector<std::unique_ptr<Particle>>(nb_particles);
    std::vector<Particle> all_particles = std::vector<Particle>(nb_particles);
    for (size_t idx = 0; idx < nb_particles; idx++)
    {
        all_particles[idx].id = local_particles[idx].id;
        all_particles[idx].position = {local_particles[idx].x, local_particles[idx].y};
        all_particles[idx].acceleration = {local_particles[idx].ax, local_particles[idx].ay};
        all_particles[idx].velocity = {local_particles[idx].vx, local_particles[idx].vy};
        all_particles[idx].mass = local_particles[idx].mass;
        //all_particles[idx] = std::make_unique<Particle>();
        //all_particles[idx]->deserialize(local_particles[idx]);
    }

    // add the entire tree of particles
    for (size_t i = 0; i < nb_particles; i++)
    {
        //quadtree->add(all_particles[i].get());
        quadtree->add(&all_particles[i]);
    }

    // update the masses of the nodes
    quadtree->update_tree_masses();
/*
    if (run_count == 0)
        quadtree->printTree(); */

// compute the forces on the particles in the entire tree for a portion of the particles
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int idx = startpos; idx < endpos; idx++)
    {
        quadtree->update_individual_barnes_hut_forces(&all_particles[idx], dtime);
        //quadtree->update_individual_barnes_hut_forces(all_particles[idx].get(), dtime);
        //all_particles[idx]->update_velocity(dtime);
        all_particles[idx].update_velocity(dtime);
        //ParticleData p_data = all_particles[idx]->updated_position(bounding_box, dtime);
        ParticleData p_data = all_particles[idx].updated_position(bounding_box, dtime);
        computed_chunk.push_back(p_data);
        // all_particles[idx]->update_position(bounding_box, dtime);
        // computed_chunk.push_back(all_particles[idx]->serialize());
    }
}

void simulation::SimulationMPI::chunked_processing(std::vector<ParticleData> &local_particles, double dtime)
{
    // save the updated particle info in a buffer to send back
    std::vector<ParticleData> updated_chunk = std::vector<ParticleData>();

    // compute the workload for each process
    if (rank != 0)
    {
        updated_chunk.reserve(endpos - startpos);
        compute_entire_workload(local_particles, updated_chunk, dtime);
    }

    gather_chunk_particles(updated_chunk, dtime);
}

void simulation::SimulationMPI::sub_tree_processing(std::vector<ParticleData> &local_particles, double dtime)
{
    // save the updated particle accelerations in a buffer to send back
    std::vector<double> local_particles_acceleration = std::vector<double>(nb_particles * 2);

    // compute the workload for each process
    if (rank != 0)
    {
        compute_workload(local_particles, local_particles_acceleration, dtime);
    }

    gather_particles(local_particles_acceleration, dtime);
    // gather_particles(local_particles, dtime);
}

void simulation::SimulationMPI::clean_mpi()
{
    MPI_Type_free(&mpi_particle_data_type);
}