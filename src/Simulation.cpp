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

// Function to scale solar system positions to fit inside the simulation box
static std::vector<std::array<double, 8>> scale_solar_system(const std::vector<std::array<double, 8>> &real_solar_system)
{
    std::vector<std::array<double, 8>> scaled_solar_system;

    double box_center_x = BOX_LEFT + BOX_WIDTH / 2;
    double box_center_y = BOX_TOP + BOX_HEIGHT / 2;

    // Scale factor based on max planet distance
    double scale_factor = (BOX_WIDTH - 20.0) / (2.0 * simulation::MAX_REAL_DISTANCE); // Keep a margin of 10px

    for (const auto &body : real_solar_system)
    {
        double x_scaled = box_center_x + (body[0] * scale_factor); // Centered around the Sun
        double y_scaled = box_center_y + (body[1] * scale_factor); // Should be 0 initially (2D view)

        scaled_solar_system.push_back({
            x_scaled, y_scaled,       // x, y positions
            body[2], body[3],         // Properly scaled velocities
            body[4],                  // Properly scaled mass
            body[5], body[6], body[7] // r, g, b colors unchanged
        });
    }

    return scaled_solar_system;
}

simulation::Simulation::~Simulation()
{
    particles.clear();
}

void simulation::Simulation::recreate_tree()
{
    // reallocate the quadtree    
    qt.reset(); 

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
    particles = std::vector<std::unique_ptr<Particle>>(num_particles + CELESTIAL_BODY_COUNT);

    std::cout << "Setting up simulation..." << std::endl;

    std::vector<std::array<double, 8>> scaled_solar_system = scale_solar_system(solar_system);

    // Assign planets based on scaled positions
    for (size_t i = 0; i < CELESTIAL_BODY_COUNT; i++)
    {
        particles[i] = std::make_unique<Particle>();
        particles[i]->id = i;
        particles[i]->position = {scaled_solar_system[i][0], scaled_solar_system[i][1]};
        particles[i]->velocity = {scaled_solar_system[i][2], scaled_solar_system[i][3]};
        particles[i]->mass = scaled_solar_system[i][4];

        float r = scaled_solar_system[i][5];
        float g = scaled_solar_system[i][6];
        float b = scaled_solar_system[i][7];

        particles[i]->color[0] = r;
        particles[i]->color[1] = g;
        particles[i]->color[2] = b;
        particles[i]->color[3] = 1.0f;
        particles[i]->radius = PLANET_RADIUS;
        if (i == 0)
        {
            std::cout << "Particle 0 is the SUN" << std::endl;
            particles[i]->radius = SUN_RADIUS;
            particles[i]->is_sun = true;
        }
        std::cout << "Particle " << i << " at " << particles[i]->position.x << ", " << particles[i]->position.y << std::endl;
        std::cout << " with velocity " << particles[i]->velocity.x << ", " << particles[i]->velocity.y << std::endl;
        std::cout << "with mass " << particles[i]->mass << std::endl;
    }
    std::cout << "Solar system setup complete" << std::endl;

    std::cout << "Setting up particle simulation..." << std::endl;
    
    Vector2<double> c = worldBounds.getCenter();
    Vector2<double> center = {c.x, c.y};

    std::cout << "Center of the world: " << center.x << ", " << center.y << std::endl;

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> radiusDist(10.0, worldBounds.width / 4.0); // Ensure particles are within bounds
    std::uniform_real_distribution<double> angleDist(0.0, 2.0 * M_PI);
    // std::uniform_real_distribution<float> xDist(worldBounds.left, worldBounds.getRight());
    // std::uniform_real_distribution<float> yDist(worldBounds.top, worldBounds.getBottom());

    // generate random color values
    std::uniform_real_distribution<float> colorDist(0.0f, 1.0f);

    // generate random particle velocities
    std::uniform_real_distribution<double> velocityDist(-1.0, 1.0);

    // we start at 1 because the first particle is the sun
    for (size_t i = CELESTIAL_BODY_COUNT; i < num_particles + CELESTIAL_BODY_COUNT; i++)
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

    for (size_t i = 1; i < CELESTIAL_BODY_COUNT; i++)
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
    for (size_t i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        Vector2<double> a_g(0.0);

        for (size_t j = 0; j < num_particles + CELESTIAL_BODY_COUNT; j++)
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
    for (size_t i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_position(worldBounds, dtime);
    }
}

void simulation::Simulation::step(double dtime)
{
    counter++;

    // leapfrog(dtime);
    qt->update_barnes_hut_forces(dtime);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        particles[i]->update_velocity(dtime);
        particles[i]->update_position(worldBounds, dtime * 0.5);
    }
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