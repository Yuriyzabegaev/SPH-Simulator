#pragma once
#include "simulation.hpp"
#include <memory>

std::shared_ptr<Simulation> initialize_simulation() {
    vec3<size_t> grid_dims = {1, 50, 60};
    vec3<double> grid_limits = {
        static_cast<double>(grid_dims.z),
        static_cast<double>(grid_dims.y),
        static_cast<double>(grid_dims.x),
    };
    grid_limits *= PARTICLE_RADIUS;
    Grid grid(std::move(grid_dims), std::move(grid_limits));

    std::vector<std::unique_ptr<Particle>> particles;
    for (auto j = 5; j < grid_dims.y - 5; j += 1) {
        for (auto i = 15; i < grid_dims.x - 15; i += 1) {
            particles.emplace_back(std::make_unique<Particle>(make_particle(
                {
                    grid_limits.z / 2,
                    PARTICLE_RADIUS / 2 + PARTICLE_RADIUS * j,
                    PARTICLE_RADIUS / 2 + PARTICLE_RADIUS * i,
                },
                RHO_0)));
        }
    }
    return std::make_shared<Simulation>(std::move(grid), std::move(particles));
}