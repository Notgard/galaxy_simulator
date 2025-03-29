#include "SimulationSDL.h"

void simulation::SimulationSDL::draw_quadtree(const TREE_TYPE::Node *node, const quadtree::Box<double> &box, int depth)
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

void simulation::SimulationSDL::draw_galaxy()
{
    for (int i = 0; i < num_particles + CELESTIAL_BODY_COUNT; i++)
    {
        draw_particle(particles[i].get());
        // draw_particle_path(particles[i].get());
    }
}

void simulation::SimulationSDL::render_tree()
{
    draw_quadtree(qt->getRoot(), qt->getBox());
}

void simulation::SimulationSDL::process_inputs()
{
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

void simulation::SimulationSDL::render()
{
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


void simulation::SimulationSDL::step(double dtime)
{
    counter++;

    leapfrog(dtime);
    recreate_tree();
}

void simulation::SimulationSDL::start()
{
    is_running = true;

    init_sdl();
    while (is_running /* && (t < t_end)  */&& (run_count < num_runs))
    {
        if (graphical)
        {
            // if (render_pause_flag)
            //     continue;
            ticks = SDL_GetTicks();
            perf_ticks = SDL_GetPerformanceCounter();
        }

        process_inputs();

        step(t);

        render();

        t += delta_time;
        run_count++;
    }
    clean_sdl();
}