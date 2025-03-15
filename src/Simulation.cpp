#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <thread>
#include <atomic>

#include "Quadtree.h"
#include "Simulation.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace quadtree;

simulation::Simulation::~Simulation()
{
#ifdef USE_SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
#endif
    particles.clear();
}

void simulation::Simulation::draw_circle(int centerX, int centerY, SDL_Color color, float radius)
{
#ifdef USE_SDL
    SDL_SetRenderDrawColor(this->renderer, color.r, color.g, color.b, color.a);

    // Loop through angles 0 to 2*PI, drawing points along the circle's circumference
    for (double angle = 0; angle < 2 * M_PI; angle += 0.05) // step size is small to smooth out the circle
    {
        int x = centerX + static_cast<int>(radius * cos(angle));
        int y = centerY + static_cast<int>(radius * sin(angle));
        SDL_RenderDrawPoint(this->renderer, x, y);
    }
#endif
}

void simulation::Simulation::draw_box(const Box<float> &box, SDL_Color color)
{
#ifdef USE_SDL
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);

    SDL_Rect rect = {
        static_cast<int>(box.getCenter().x - box.getSize().x / 2),
        static_cast<int>(box.getCenter().y - box.getSize().y / 2),
        static_cast<int>(box.getSize().x),
        static_cast<int>(box.getSize().y)};
    SDL_RenderDrawRect(renderer, &rect);
#endif
}

void simulation::Simulation::draw_cross(int x, int y, int size)
{
#ifdef USE_SDL
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    // Draw horizontal line
    SDL_RenderDrawLine(renderer, x - size / 2, y, x + size / 2, y);
    // Draw vertical line
    SDL_RenderDrawLine(renderer, x, y - size / 2, x, y + size / 2);
#endif
}

void simulation::Simulation::draw_text(const std::string &text, int x, int y)
{
#ifdef USE_SDL
    SDL_Color color = {255, 255, 255, 255};
    SDL_Surface *surface = TTF_RenderText_Solid(font, text.c_str(), color);
    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);

    int texW = 0;
    int texH = 0;
    SDL_QueryTexture(texture, NULL, NULL, &texW, &texH);
    SDL_Rect dstrect = {x, y, texW, texH};

    SDL_RenderCopy(renderer, texture, NULL, &dstrect);

    SDL_DestroyTexture(texture);
    SDL_FreeSurface(surface);
#endif
}

void simulation::Simulation::draw_particle(Particle *particle)
{

    SDL_Color color = {static_cast<Uint8>(particle->color.r * 255),
                       static_cast<Uint8>(particle->color.g * 255),
                       static_cast<Uint8>(particle->color.b * 255),
                       static_cast<Uint8>(particle->color.a * 255)};
    draw_circle(particle->position.x, particle->position.y, color, particle->radius);
}

void simulation::Simulation::draw_particle_path(Particle *particle)
{
    // Loop through the position history and draw lines between points
    for (size_t i = 1; i < particle->position_history.size(); ++i)
    {
        const glm::vec2 &p1 = particle->position_history[i - 1];
        const glm::vec2 &p2 = particle->position_history[i];

        SDL_SetRenderDrawColor(this->renderer,
                               static_cast<Uint8>(particle->color.r * 255),
                               static_cast<Uint8>(particle->color.g * 255),
                               static_cast<Uint8>(particle->color.b * 255),
                               static_cast<Uint8>(particle->color.a * 255));

        int x1 = static_cast<int>(p1.x);
        int y1 = static_cast<int>(p1.y);
        int x2 = static_cast<int>(p2.x);
        int y2 = static_cast<int>(p2.y);

        // Draw the line between the points
        SDL_RenderDrawLine(this->renderer, x1, y1, x2, y2);
    }
}

void simulation::Simulation::draw_quadtree(const TREE_TYPE::Node *node, const quadtree::Box<float> &box, int depth = 0)
{
    if (!node)
        return;

    SDL_Color color = {0, 0, 255};

    draw_box(box, color);

    // Recursively draw child nodes
    for (int i = 0; i < QUADTREE_MAX_VALUES; ++i)
    {
        if (node->children[i])
        {
            draw_quadtree(node->children[i].get(), qt->computeBox(box, i), depth + 1);
        }
    }

    // Draw center of mass
    if (render_center_of_mass_flag)
    {
        if (!qt->isLeaf(node))
            draw_cross(node->centerOfMass.x, node->centerOfMass.y, 5);
    }
}

void simulation::Simulation::draw_galaxy()
{
    // #pragma omp parallel for
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        draw_particle(particles[i].get());
        // draw_particle_path(particles[i].get());
    }
}

void simulation::Simulation::render_tree()
{
    draw_quadtree(qt->getRoot(), qt->getBox());
}

void simulation::Simulation::recreate_tree()
{
    // reallocate the quadtree
    worldBounds = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
    qt = std::make_unique<TREE_TYPE>(worldBounds, getBoxFunc);

    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        qt->add(particles[i].get());
    }

    // update the masses of the nodes
    qt->update_tree_masses();

    // dump the tree in debug mode
    if (counter == 0 && DEBUG)
    {
        qt->printTree();
    }
}

void simulation::Simulation::recreate_tree(std::unique_ptr<TREE_TYPE>& quadtree, quadtree::Box<float> box)
{
    box = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
    quadtree = std::make_unique<TREE_TYPE>(box, getBoxFunc);

    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        qt->add(particles[i].get());
    }

    // update the masses of the nodes
    qt->update_tree_masses();
}

void simulation::Simulation::setup()
{
    std::cout << "Setting up particle simulation..." << std::endl;

    worldBounds = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
    Vector2<float> c = worldBounds.getCenter();
    glm::vec2 center = {c.x, c.y};
    std::cout << "Center of the world: " << center.x << ", " << center.y << std::endl;

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> radiusDist(10.0f, worldBounds.width / 4.0f); // Ensure particles are within bounds
    std::uniform_real_distribution<float> angleDist(0.0f, 2.0f * M_PI);
    // std::uniform_real_distribution<float> xDist(worldBounds.left, worldBounds.getRight());
    // std::uniform_real_distribution<float> yDist(worldBounds.top, worldBounds.getBottom());

    // generate random color values
    std::uniform_real_distribution<float> colorDist(0.0f, 1.0f);

    // generate random particle velocities
    std::uniform_real_distribution<float> velocityDist(-1.0f, 1.0f);

    // we start at 1 because the first particle is the sun
    for (int i = CELESTIAL_BODY_COUNT; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        float angle = angleDist(gen);
        float radius = radiusDist(gen);

        glm::vec2 position = center + glm::vec2(radius * cos(angle), radius * sin(angle));
        std::cout << "Particle " << i << " at " << position.x << ", " << position.y << std::endl;
        glm::vec4 color = {colorDist(gen), colorDist(gen), colorDist(gen), 1.0f};

        particles[i] = std::make_unique<Particle>();

        particles[i]->id = i;
        particles[i]->position = position;
        particles[i]->mass = PARTICLE_MASS;
        particles[i]->velocity = {velocityDist(gen), velocityDist(gen)};

        // Compute circular velocity
        float r = glm::distance(position, center);
        if (r > 0)
        {
            float v_circular = sqrt(G * particles[0]->mass / r);
            particles[i]->velocity.x = -v_circular * sin(angle);
            particles[i]->velocity.y = v_circular * cos(angle);
        }

        particles[i]->color = color;

        particle_count++;
    }

    recreate_tree();

#ifdef USE_SDL
    init_sdl();
#endif

    std::cout << "Particle simulation setup complete" << std::endl;
}

void simulation::Simulation::init_sdl()
{
#ifndef USE_SDL
    return;
#endif
    if (!graphical)
        return;

    std::cout << "Initializing SDL..." << std::endl;
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
        exit(1);
    }

    // ttf init
    if (TTF_Init() < 0)
    {
        std::cerr << "TTF could not initialize! TTF_Error: " << TTF_GetError() << std::endl;
        exit(1);
    }

    window = SDL_CreateWindow("Particle Simulation",
                              SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                              WINDOW_WIDTH, WINDOW_HEIGHT,
                              SDL_WINDOW_SHOWN);
    if (!window)
    {
        std::cerr << "Window could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        exit(1);
    }

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer)
    {
        std::cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        exit(1);
    }
    SDL_SetHint(SDL_HINT_RENDER_DRIVER, "openGL");

    font = TTF_OpenFont("JetBrainsMono-ExtraLight.ttf", 21);
    if (!font)
    {
        std::cerr << "Error loading font: " << TTF_GetError() << std::endl;
        exit(1);
    }
    std::cout << "SDL initialized" << std::endl;
}

// first order leapfrog integration
void simulation::Simulation::leapfrog(double dtime)
{
// move by half step, update forces, update velocities and move by half step again for the full step
#pragma omp parallel for
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_position(worldBounds, dtime * 0.5);
    }

    qt->update_barnes_hut_forces(dtime);

#pragma omp parallel for
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_velocity(dtime);
    }

#pragma omp parallel for
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_position(worldBounds, dtime * 0.5);
    }
}

void simulation::Simulation::brute_force(double dtime)
{
#pragma omp parallel for
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        glm::dvec2 a_g(0.0);

        for (int j = 0; j < num_particles + CELESTIAL_BODY_COUNT; j++)
        {
            if (particles[i]->id != particles[j]->id && !particles[i]->is_sun)
            {
                glm::dvec2 r = particles[i]->position - particles[j]->position;
                double r_mag = glm::length(r);

                if (r_mag > 0)
                {
                    double acceleration = -G * particles[j]->mass / (r_mag * r_mag);
                    glm::dvec2 r_unit_vector = r / r_mag;
                    a_g += acceleration * r_unit_vector;
                }
            }
        }

        particles[i]->velocity += a_g * dtime;
    }

#pragma omp parallel for
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_position(worldBounds, dtime);
    }
}

void simulation::Simulation::step(double dtime)
{
    counter++;

    leapfrog(dtime);
    // qt->update_barnes_hut_forces(dtime);
    recreate_tree();
    // brute_force(dtime);
}

void simulation::Simulation::process_inputs()
{
#ifndef USE_SDL
    return;
#endif
    if (!graphical)
        return;
    SDL_Event event;
    while (SDL_PollEvent(&event))
    {

        switch (event.type)
        {
        case SDL_QUIT:
            is_running = false;
            break;
        case SDL_KEYDOWN:
            if (event.key.keysym.sym == SDLK_ESCAPE)
            {
                std::cout << "Exiting simulation" << std::endl;
                is_running = false;
            }
            else if (event.key.keysym.sym == SDLK_c)
            {
                std::cout << "Toggling center of mass rendering" << std::endl;
                render_center_of_mass_flag = !render_center_of_mass_flag;
            }
            else if (event.key.keysym.sym == SDLK_r)
            {
                std::cout << "Toggling quadtree rendering" << std::endl;
                render_tree_flag = !render_tree_flag;
            }
            else if (event.key.keysym.sym == SDLK_p)
            {
                std::cout << "Pausing simulation" << std::endl;
                render_pause_flag = !render_pause_flag;
            }
            break;
        }
    }
}

void simulation::Simulation::render()
{
#ifndef USE_SDL
    return;
#endif
    if (!graphical)
        return;
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    if (render_tree_flag)
    {
        render_tree();
    }

    draw_galaxy();

    // End frame timing
    Uint32 endTicks = SDL_GetTicks();
    Uint64 endPerf = SDL_GetPerformanceCounter();
    Uint64 framePerf = endPerf - perf_ticks;
    float frameTime = (endTicks - ticks) / 1000.0f;
    total_frame_ticks += endTicks - ticks;

    // Strings to display
    std::string fps = "Current FPS: " + std::to_string(1.0f / frameTime);
    std::string avg = "Average FPS: " + std::to_string(1000.0f / ((float)total_frame_ticks / run_count));
    std::string time = "Time (dt): " + std::to_string(t);
    std::string runs = "Runs: " + std::to_string(run_count);

    draw_text(fps, 10, 10);
    draw_text(avg, 10, 30);
    draw_text(time, 10, 50);
    draw_text(runs, 10, 70);

    SDL_RenderPresent(renderer);
}

void simulation::Simulation::clean_sdl()
{
#ifndef USE_SDL
    return;
#endif
    if (!graphical)
        return;

    // cleanup sdl
    TTF_CloseFont(font);
    font = nullptr;
    SDL_DestroyRenderer(renderer);
    renderer = nullptr;
    SDL_DestroyWindow(window);
    window = nullptr;
    TTF_Quit();
    SDL_Quit();
}

void simulation::Simulation::start()
{
    is_running = true;

    while (is_running && (t < t_end) && (run_count < num_runs))
    {
#ifdef USE_SDL
        if (graphical)
        {
            if (render_pause_flag)
                continue;
            ticks = SDL_GetTicks();
            perf_ticks = SDL_GetPerformanceCounter();
        }
#endif
        process_inputs();

        step(t);

        render();

        t += delta_time;
        run_count++;
    }

    clean_sdl();
}

void simulation::Simulation::init_mpi()
{
    MPI_Aint displacements[6] = {
        offsetof(ParticleData, id), offsetof(ParticleData, x), offsetof(ParticleData, y),
        offsetof(ParticleData, vx), offsetof(ParticleData, vy), offsetof(ParticleData, mass)};

    int block_lengths[6] = {1, 1, 1, 1, 1, 1};
    MPI_Datatype dtypes[6] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_create_struct(6, block_lengths, displacements, dtypes, &mpi_particle_data_type);
    MPI_Type_commit(&mpi_particle_data_type);
}

void simulation::Simulation::mpi_start(int rank, int size)
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

        if(rank == 0)
            process_inputs();

        // only the compute step is parallelized throughout the mpi processes
        mpi_step(t, rank, size);

        if (rank == 0)
            render();

        t += delta_time;
        run_count++;
    }

    if (rank == 0)
        clean_sdl();
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
void simulation::Simulation::mpi_step(double dtime, int rank, int size)
{
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Rank " << rank << " is running" << std::endl;

    // each process create the tree and update the masses
    //recreate_tree();
    quadtree::Box<float> bounding_box;
    std::unique_ptr<TREE_TYPE> quadtree;

    recreate_tree(quadtree, bounding_box);

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
        if(particles[i]->id == all_particles[i].id)
            particles[i]->deserialize(all_particles[i]);
    }
#endif
}

void simulation::Simulation::distribute_subtrees()
{
#ifdef USE_MPI

#endif
}