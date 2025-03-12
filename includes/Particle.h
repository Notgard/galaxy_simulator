#pragma once

#include "Config.h"
#include "Box.h"

#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtx/string_cast.hpp>

#include <vector>
#include <iostream>

struct Particle
{
    int id;
    glm::vec<2, double> position;
    glm::vec4 color;

    glm::vec<2, double> speed = {0.0, 0.0};
    glm::vec<2, double> acceleration = {0.0, 0.0};
    glm::vec<2, double> velocity = {0.0, 0.0};
    double mass = PARTICLE_MASS;
    std::vector<glm::vec<2, double>> position_history = std::vector<glm::vec<2, double>>(PARTICLE_POSITION_HISTORY_SIZE);
    float radius = PARTICLE_RADIUS;
    bool is_sun = false;

    void update_position(quadtree::Box<float> worldBounds, double dtime)
    {
        position_history.push_back(position);
        if (position_history.size() > PARTICLE_POSITION_HISTORY_SIZE) {
            position_history.erase(position_history.begin());
        }

        position += velocity * dtime;
        acceleration = {0.0f, 0.0f};

        //position before wrapping around
        //std::cout << "Particle " << id << ", Position before wrapping around: " << position.x << ", " << position.y << std::endl;

        // Wrap around to the other side if out of bounds
        float right = worldBounds.getRight();
        float bottom = worldBounds.getBottom();

        //print world bounds
        //std::cout << "World bounds: " << worldBounds.left << ", " << worldBounds.top << ", " << worldBounds.width << ", " << worldBounds.height << std::endl;
    
        float worldWidth = worldBounds.width;
        float worldHeight = worldBounds.height;
    
        // Correct wrapping using modulo to handle large out-of-bounds values
        if (position.x < worldBounds.left)
            position.x = right - std::fmod(worldBounds.left - position.x, worldWidth);
        else if (position.x > right)
            position.x = worldBounds.left + std::fmod(position.x - right, worldWidth);
    
        if (position.y < worldBounds.top)
            position.y = bottom - std::fmod(worldBounds.top - position.y, worldHeight);
        else if (position.y > bottom)
            position.y = worldBounds.top + std::fmod(position.y - bottom, worldHeight);

        //position after wrapping around
        //std::cout << "Particle " << id << ", Position after wrapping around: " << position.x << ", " << position.y << std::endl;
    }

    void update_velocity(double dtime)
    {
        velocity += acceleration * dtime;
    }

    void update_acceleration(glm::vec2 acceleration)
    {
        this->acceleration += acceleration;
    }
};