#pragma once

#include "Simulation.h"
#include "SDLUtils.h"

namespace simulation
{
    class SimulationSDL : virtual public Simulation, public SDLUtils
    {
    public:
        SimulationSDL(int num_particles, int num_runs = -1, bool graphical = false)
            : Simulation(num_particles, num_runs)
        {
            this->graphical = graphical;
        }

        ~SimulationSDL()
        {
        }

        void start() override;
        void step(double dtime) override;
    protected:
        
        void draw_quadtree(const TREE_TYPE::Node *node, const quadtree::Box<double> &box, int depth = 0);
        void draw_galaxy();
        void render_tree();

        void process_inputs();
        void render();
    };
}
