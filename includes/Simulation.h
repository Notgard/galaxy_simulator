#pragma once

#include "Box.h"
#include "Quadtree.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <chrono>
#include <iostream>
#include <thread>
#include <atomic>

#define BOX_GETTER_TYPE std::function<quadtree::Box<double>(Particle *)>
#define TREE_TYPE quadtree::Quadtree<Particle *, BOX_GETTER_TYPE>

namespace simulation
{
    static const double DIST_SCALE = 2.67e-9; // Distance scaling factor
    static const double TIME_SCALE = 1e5;     // Time step in seconds
    static const double VEL_SCALE = DIST_SCALE * TIME_SCALE;

    // x, y, vx, vy, mass, r, g, b
    static const std::vector<std::array<double, 8>> solar_system = {
        {0.0, 0.0, 0.0, 0.0, SUN_MASS /* 1.989e30 */, 1.0, 1.0, 0.0},                      // Sun (Yellow)
        {57.909e9, 0.0, 0.0, 47.36e-5 * VEL_SCALE, 100 /* 0.33011e24 */, 0.6, 0.6, 0.6},   // Mercury (Gray)
        {108.209e9, 0.0, 0.0, 35.02e-5 * VEL_SCALE, 600 /* 4.8675e24 */, 0.9, 0.7, 0.3},   // Venus (Pale Yellow)
        {149.596e9, 0.0, 0.0, 29.78e-5 * VEL_SCALE, 650 /* 5.9724e24 */, 0.0, 0.5, 1.0},   // Earth (Blue-Green)
        {227.923e9, 0.0, 0.0, 24.07e-5 * VEL_SCALE, 200 /* 0.64171e24 */, 0.8, 0.2, 0.2},  // Mars (Red)
        {778.570e9, 0.0, 0.0, 13e-5 * VEL_SCALE, 40000 /* 1898.19e24 */, 0.9, 0.6, 0.4},   // Jupiter (Orange-White)
        {1433.529e9, 0.0, 0.0, 9.68e-5 * VEL_SCALE, 10000 /* 568.34e24 */, 0.8, 0.7, 0.5}, // Saturn (Pale Yellow)
        {2872.463e9, 0.0, 0.0, 6.80e-5 * VEL_SCALE, 4000 /* 86.813e24 */, 0.5, 0.8, 1.0},  // Uranus (Cyan)
        {4495.060e9, 0.0, 0.0, 5.43e-5 * VEL_SCALE, 5500 /* 102.413e24 */, 0.3, 0.3, 1.0}  // Neptune (Deep Blue)
    };

    // Max real distance (Neptune's distance in meters)
    constexpr double MAX_REAL_DISTANCE = 4.495060e12;

    class Simulation
    {
    public:
        Simulation(int num_particles, int num_runs = -1)
            : num_particles(num_particles), num_runs(num_runs),
              getBoxFunc([](Particle *p)
                         { return quadtree::Box<double>(
                               quadtree::Vector2<double>(p->position.x, p->position.y),
                               quadtree::Vector2<double>(0.0, 0.0)); }),
              qt(std::make_unique<TREE_TYPE>(
                  worldBounds, getBoxFunc))
        {
        }

        ~Simulation();

        void recreate_tree();

        void setup();
        virtual void start();
        virtual void step(double dtime);
        void stop() { is_running = false; }

        void start_timer() { start_time = std::chrono::steady_clock::now(); }
        void end_timer() { end_time = std::chrono::steady_clock::now(); }
        void print_time()
        {
            std::cout << "Time: (seconds) " << std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count() << "s" << std::endl;
            std::cout << "Time: (milliseconds) " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;
            std::cout << "Time: (microseconds) " << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "us" << std::endl;
        }

    protected:
        BOX_GETTER_TYPE getBoxFunc;

        quadtree::Box<double> worldBounds = quadtree::Box<double>(BOX_LEFT, BOX_TOP, BOX_WIDTH, BOX_HEIGHT);
        std::unique_ptr<TREE_TYPE> qt;

        std::vector<std::unique_ptr<Particle>> particles;

        void leapfrog(double dtime);
        void brute_force(double dtime);

        unsigned long int num_runs; // Number of runs (-1 for infite number of runs till manually stopped)
        int num_particles;
        int counter = 0;
        int particle_count = 0;

        bool is_running = false;

        double t_0 = 0.0;
        double t = t_0;
        unsigned long int run_count = 0;

        // keep the start and and time of the simulation from std::chrono
        std::chrono::time_point<std::chrono::steady_clock> start_time;
        std::chrono::time_point<std::chrono::steady_clock> end_time;
    };
}
