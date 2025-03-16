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
#ifdef USE_SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
#endif
    particles.clear();
}
#ifdef USE_SDL
void simulation::Simulation::draw_circle(int centerX, int centerY, SDL_Color color, float radius)
{

    SDL_SetRenderDrawColor(this->renderer, color.r, color.g, color.b, color.a);

    // Loop through angles 0 to 2*PI, drawing points along the circle's circumference
    for (double angle = 0; angle < 2 * M_PI; angle += 0.05) // step size is small to smooth out the circle
    {
        int x = centerX + static_cast<int>(radius * cos(angle));
        int y = centerY + static_cast<int>(radius * sin(angle));
        SDL_RenderDrawPoint(this->renderer, x, y);
    }
}

void simulation::Simulation::draw_box(const Box<float> &box, SDL_Color color)
{
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);

    SDL_Rect rect = {
        static_cast<int>(box.getCenter().x - box.getSize().x / 2),
        static_cast<int>(box.getCenter().y - box.getSize().y / 2),
        static_cast<int>(box.getSize().x),
        static_cast<int>(box.getSize().y)};
    SDL_RenderDrawRect(renderer, &rect);
}

void simulation::Simulation::draw_cross(int x, int y, int size)
{
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    // Draw horizontal line
    SDL_RenderDrawLine(renderer, x - size / 2, y, x + size / 2, y);
    // Draw vertical line
    SDL_RenderDrawLine(renderer, x, y - size / 2, x, y + size / 2);
}

void simulation::Simulation::draw_text(const std::string &text, int x, int y)
{
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
#endif

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
        if(DEBUG)
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

#ifdef USE_SDL
    init_sdl();
#endif

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
    // brute_force(dtime);
    recreate_tree();
}

void simulation::Simulation::start()
{
    is_running = true;

    while (is_running && (t < t_end) && (run_count < num_runs))
    {
#ifdef USE_SDL
        if (graphical)
        {
            // if (render_pause_flag)
            //     continue;
            ticks = SDL_GetTicks();
            perf_ticks = SDL_GetPerformanceCounter();
        }

        process_inputs();
#endif
        step(t);
#ifdef USE_SDL
        render();
#endif
        t += delta_time;
        run_count++;
    }
#ifdef USE_SDL
    clean_sdl();
#endif
}