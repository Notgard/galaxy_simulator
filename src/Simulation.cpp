#include <cassert>
#include <chrono>
#include <iostream>
#include <random>

#include "Quadtree.h"
#include "Simulation.h"

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

void simulation::Simulation::recreate_tree()
{
    worldBounds = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
    qt = std::make_unique<TREE_TYPE>(worldBounds, getBoxFunc);

    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        qt->add(particles[i].get());
        /*         qt->printTree();
                std::cout << std::endl; */
    }
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
    // if (!qt->isLeaf(node))
    //    draw_cross(node->centerOfMass.x, node->centerOfMass.y, 5);
}

void simulation::Simulation::render_tree()
{
    draw_quadtree(qt->getRoot(), qt->getBox());
}

void simulation::Simulation::setup()
{
    std::cout << "Setting up particle simulation..." << std::endl;

    worldBounds = quadtree::Box<float>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(6);
    std::uniform_real_distribution<float> xDist(worldBounds.left, worldBounds.getRight());
    std::uniform_real_distribution<float> yDist(worldBounds.top, worldBounds.getBottom());

    // generate random color values
    std::uniform_real_distribution<float> colorDist(0.0f, 1.0f);

    // generate random particle velocities
    std::uniform_real_distribution<float> velocityDist(-1.0f, 1.0f);

    // we start at 1 because the first particle is the sun
    for (int i = CELESTIAL_BODY_COUNT; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        glm::vec2 position = {xDist(gen), yDist(gen)};
        std::cout << "Particle " << i << " at " << position.x << ", " << position.y << std::endl;
        glm::vec4 color = {colorDist(gen), colorDist(gen), colorDist(gen), 1.0f};

        particles[i] = std::make_unique<Particle>();

        particles[i]->id = i;
        particles[i]->position = position;
        particles[i]->mass = PARTICLE_MASS;
        particles[i]->velocity = {velocityDist(gen), velocityDist(gen)};
        particles[i]->color = color;

        particle_count++;
    }

    for (size_t i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        double dx = particles[i]->position.x - particles[0]->position.x; // Distance from Sun center
        double dy = particles[i]->position.y - particles[0]->position.y;
        double r = sqrt(dx * dx + dy * dy);

        if (r > 0)
        {
            double v_circular = sqrt(G * particles[0]->mass / r); // Orbital velocity
            particles[i]->velocity.x = -v_circular * (dy / r);    // Perpendicular to radius
            particles[i]->velocity.y = v_circular * (dx / r);
        }
        std::cout << "Particle " << particles[i]->id << " velocity " << particles[i]->velocity.x << ", " << particles[i]->velocity.y << std::endl;
    }

    recreate_tree();

    if (graphical)
    {
        init_sdl();
    }
}

void simulation::Simulation::init_sdl()
{
#ifdef USE_SDL
    std::cout << "Initializing SDL..." << std::endl;
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
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

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer)
    {
        std::cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        exit(1);
    }
#endif
}

void simulation::Simulation::leapfrog(double dtime)
{
    // move by half step, update forces, update velocities and move by half step again for the full step
    for (auto &particle : particles)
    {
        particle->update_position(worldBounds, dtime * 0.5f);
    }

    qt->update_barnes_hut_forces(dtime);

    for (auto &particle : particles)
    {
        particle->update_velocity(dtime);
    }

    for (auto &particle : particles)
    {
        particle->update_position(worldBounds, dtime * 0.5f);
    }
}

void simulation::Simulation::brute_force(double dtime)
{
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        glm::vec<2, double> a_g = {0.0f, 0.0f};
        for (int j = 0; j < num_particles + CELESTIAL_BODY_COUNT; j++)
        {
            if (particles[i]->id != particles[j]->id && !particles[i]->is_sun)
            {
                glm::vec<2, double> r;
                r = particles[i]->position - particles[j]->position;
                double r_mag = sqrt(r.x * r.x + r.y * r.y);
                double acceleration = -1.0 * G * (particles[j]->mass) / (r_mag * r_mag);
                glm::vec<2, double> r_unit_vector = {r.x / r_mag, r.y / r_mag};
                a_g += acceleration * r_unit_vector;
            }
        }
        particles[i]->velocity += a_g * dtime;
    }

    for (auto &particle : particles)
    {
        particle->update_position(worldBounds, dtime);
    }
}

void simulation::Simulation::step(double dtime)
{

    if (counter == 0)
    {
        qt->printTree();
    }

    counter++;

#ifdef USE_SDL
    render_tree();

    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        draw_particle(particles[i].get());
        // draw_particle_path(particles[i].get());
    }
#endif

    leapfrog(dtime);
    // qt->update_barnes_hut_forces(dtime);
    //brute_force(dtime);

    recreate_tree();
}

void simulation::Simulation::start()
{
#ifdef USE_SDL
    SDL_Event event;
#endif
    is_running = true;
    while (is_running && (t < t_end) && (run_count < num_runs))
    {
#ifdef USE_SDL
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                is_running = false;
            }
        }
        // Clear the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
#endif
        step(t);

        //  increase delta time step
        t += dt;

#ifdef USE_SDL
        // Render the frame
        SDL_RenderPresent(renderer);
        SDL_Delay(16); // ~60 FPS
#endif
        run_count++;
    }
#ifdef USE_SDL
    // Cleanup
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
#endif
}