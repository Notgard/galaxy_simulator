#pragma once

#include "Config.h"
#include "Box.h"

#ifdef USE_GLM
#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/ext/vector_double1_precision.hpp>
#include <glm/ext/vector_double2_precision.hpp>
#include <glm/ext/vector_double3_precision.hpp>
#endif

#include "Vector2d.h"

#include <vector>
#include <iostream>

using namespace quadtree;

struct ParticleData
{
    int id;
    double x, y, ax, ay, vx, vy, mass;
};

struct Particle
{
    int id;
    float color[4];
    Vector2<double> position;
    Vector2<double> acceleration = {0.0, 0.0};
    Vector2<double> velocity = {0.0, 0.0};

    double mass = PARTICLE_MASS;
    float radius = PARTICLE_RADIUS;
    bool is_sun = false;

    ParticleData serialize()
    {
        ParticleData data;
        data.id = id;
        data.x = position.x;
        data.y = position.y;
        data.ax = acceleration.x;
        data.ay = acceleration.y;
        data.vx = velocity.x;
        data.vy = velocity.y;
        data.mass = mass;
        return data;
    }

    void deserialize(ParticleData data)
    {
        id = data.id;
        position = {data.x, data.y};
        acceleration = {data.ax, data.ay};
        velocity = {data.vx, data.vy};
        mass = data.mass;
    }

    void accumulate(ParticleData data) {
        acceleration.x += data.ax;
        acceleration.y += data.ay;
    }

    void update_position(quadtree::Box<float> worldBounds, double dtime)
    {
        /*         position_history.push_back(position);
                if (position_history.size() > PARTICLE_POSITION_HISTORY_SIZE) {
                    position_history.erase(position_history.begin());
                } */

        position += velocity * dtime;

        // position before wrapping around
        // std::cout << "Particle " << id << ", Position before wrapping around: " << position.x << ", " << position.y << std::endl;

        // Wrap around to the other side if out of bounds
        float right = worldBounds.getRight();
        float bottom = worldBounds.getBottom();

        // print world bounds
        // std::cout << "World bounds: " << worldBounds.left << ", " << worldBounds.top << ", " << worldBounds.width << ", " << worldBounds.height << std::endl;

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

        // position after wrapping around
        // std::cout << "Particle " << id << ", Position after wrapping around: " << position.x << ", " << position.y << std::endl;
    }

    void update_velocity(double dtime)
    {
        velocity += acceleration * dtime;
        acceleration = {0.0, 0.0};
    }

    void update_acceleration(Vector2<double> acceleration)
    {
        this->acceleration += acceleration;
    }
};