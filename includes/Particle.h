#pragma once

#include "Config.h"
#include "Box.h"

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

    void update(ParticleData data) {
        position = {data.x, data.y};
        velocity = {data.vx, data.vy};
    }

    void update_position(quadtree::Box<double> worldBounds, double dtime)
    {
        Vector2<double> pos;
        pos.x = this->position.x;
        pos.y = this->position.y;
        pos += velocity * dtime;

        // position before wrapping around
        // std::cout << "Particle " << id << ", Position before wrapping around: " << position.x << ", " << position.y << std::endl;

        // Wrap around to the other side if out of bounds
        double right = worldBounds.getRight();
        double bottom = worldBounds.getBottom();

        // print world bounds
        // std::cout << "World bounds: " << worldBounds.left << ", " << worldBounds.top << ", " << worldBounds.width << ", " << worldBounds.height << std::endl;

        double worldWidth = worldBounds.width;
        double worldHeight = worldBounds.height;

        // Correct wrapping using modulo to handle large out-of-bounds values
        if (pos.x < worldBounds.left)
            pos.x = right - std::fmod(worldBounds.left - pos.x, worldWidth);
        else if (pos.x > right)
            pos.x = worldBounds.left + std::fmod(pos.x - right, worldWidth);

        if (pos.y < worldBounds.top)
            pos.y = bottom - std::fmod(worldBounds.top - pos.y, worldHeight);
        else if (pos.y > bottom)
            pos.y = worldBounds.top + std::fmod(pos.y - bottom, worldHeight);

        this->position.x = pos.x;
        this->position.y = pos.y;

        // position after wrapping around
        // std::cout << "Particle " << id << ", Position after wrapping around: " << position.x << ", " << position.y << std::endl;
    }

    ParticleData updated_position(quadtree::Box<double> worldBounds, double dtime)
    {
        ParticleData data;

        Vector2<double> pos;
        pos.x = this->position.x;
        pos.y = this->position.y;

        pos += velocity * dtime;

        // position before wrapping around
        // std::cout << "Particle " << id << ", Position before wrapping around: " << position.x << ", " << position.y << std::endl;

        // Wrap around to the other side if out of bounds
        double right = worldBounds.getRight();
        double bottom = worldBounds.getBottom();

        // print world bounds
        // std::cout << "World bounds: " << worldBounds.left << ", " << worldBounds.top << ", " << worldBounds.width << ", " << worldBounds.height << std::endl;

        double worldWidth = worldBounds.width;
        double worldHeight = worldBounds.height;

        // Correct wrapping using modulo to handle large out-of-bounds values
        if (pos.x < worldBounds.left)
            pos.x = right - std::fmod(worldBounds.left - pos.x, worldWidth);
        else if (pos.x > right)
            pos.x = worldBounds.left + std::fmod(pos.x - right, worldWidth);

        if (pos.y < worldBounds.top)
            pos.y = bottom - std::fmod(worldBounds.top - pos.y, worldHeight);
        else if (pos.y > bottom)
            pos.y = worldBounds.top + std::fmod(pos.y - bottom, worldHeight);

        this->position.x = pos.x;
        this->position.y = pos.y;

        data.id = id;
        data.x = position.x;
        data.y = position.y;
        data.ax = acceleration.x;
        data.ay = acceleration.y;
        data.vx = velocity.x;
        data.vy = velocity.y;

        return data;
    }

    void update_velocity(double dtime)
    {
        velocity += acceleration * dtime;
        acceleration = {0.0, 0.0};
    }

    void update_acceleration(Vector2<double> accel)
    {
        this->acceleration += accel;
    }
};