#include "SimulationMPISDL.h"

void simulation::SimulationMPISDL::mpi_start()
{
    is_running = true;

    if (rank == 0)
        init_sdl();

    while (is_running && (run_count < num_runs))
    {
        if (rank == 0)
        {
            ticks = SDL_GetTicks();
            perf_ticks = SDL_GetPerformanceCounter();
            process_inputs();
        }

        // only the compute step is parallelized throughout the mpi processes
        mpi_step(t);

        if (rank == 0)
        {
            recreate_tree();
            render();
        }

        t += delta_time;

        run_count++;
    }

    if (rank == 0)
        clean_sdl();
}