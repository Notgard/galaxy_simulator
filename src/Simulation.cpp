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

using namespace quadtree;

simulation::Simulation::~Simulation()
{
    particles.clear();
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

void simulation::Simulation::setup()
{
    std::cout << "Setting up particle simulation..." << std::endl;

    worldBounds = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
    Vector2<float> c = worldBounds.getCenter();
    Vector2<double> center = {c.x, c.y};

    std::cout << "Center of the world: " << center.x << ", " << center.y << std::endl;

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(5);
    std::uniform_real_distribution<double> radiusDist(10.0f, worldBounds.width / 4.0f); // Ensure particles are within bounds
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
        double angle = angleDist(gen);
        double radius = radiusDist(gen);

        Vector2<double> position = center + Vector2<double>(radius * cos(angle), radius * sin(angle));
        if (DEBUG)
            std::cout << "Particle " << i << " at " << position.x << ", " << position.y << std::endl;
        float color[4] = {colorDist(gen), colorDist(gen), colorDist(gen), 1.0f};

        particles[i] = std::make_unique<Particle>();

        particles[i]->id = i;
        particles[i]->position = position;
        particles[i]->mass = PARTICLE_MASS;
        particles[i]->velocity = {velocityDist(gen), velocityDist(gen)};

        // Compute circular velocity
        double r = Vector2<double>::distance(position, center);
        if (r > 0)
        {
            double v_circular = sqrt(G * particles[0]->mass / r);
            particles[i]->velocity.x = -v_circular * sin(angle);
            particles[i]->velocity.y = v_circular * cos(angle);
        }

        particles[i]->color[0] = color[0];
        particles[i]->color[1] = color[1];
        particles[i]->color[2] = color[2];
        particles[i]->color[3] = color[3];

        particle_count++;
    }

    for (int i = 1; i < CELESTIAL_BODY_COUNT; i++)
    {
        double dx = particles[i]->position.x - 410; // Distance from center
        double dy = particles[i]->position.y - 410;
        double r = sqrt(dx * dx + dy * dy);

        if (r > 0)
        {
            double v_circular = sqrt(G * SUN_MASS / r);        // Orbital velocity formula
            particles[i]->velocity.x = -v_circular * (dy / r); // Perpendicular to radius
            particles[i]->velocity.y = v_circular * (dx / r);
        }
    }

    recreate_tree();

    std::cout << "Particle simulation setup complete" << std::endl;
}

// first order leapfrog integration
void simulation::Simulation::leapfrog(double dtime)
{
    // move by half step, update forces, update velocities and move by half step again for the full step
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_position(worldBounds, dtime * 0.5);
    }

    qt->update_barnes_hut_forces(dtime);

    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_velocity(dtime);
        particles[i]->update_position(worldBounds, dtime * 0.5);
    }
}

void simulation::Simulation::brute_force(double dtime)
{
#pragma omp parallel for
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        Vector2<double> a_g(0.0);

        for (int j = 0; j < num_particles + CELESTIAL_BODY_COUNT; j++)
        {
            if (particles[i]->id != particles[j]->id && !particles[i]->is_sun)
            {
                Vector2<double> r = particles[i]->position - particles[j]->position;
                double r_mag = Vector2<double>::length(r);

                if (r_mag > 0)
                {
                    double acceleration = -G * particles[j]->mass / (r_mag * r_mag);
                    Vector2<double> r_unit_vector = r / r_mag;
                    a_g += r_unit_vector * acceleration;
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
    // brute_force(dtime);
    recreate_tree();
}

void simulation::Simulation::start()
{
    is_running = true;

    while (is_running /* && (t < t_end) */ && (run_count < num_runs))
    {
        step(t);
        t += delta_time;
        run_count++;
    }
}