#include "../core/simulation.hpp"
#include <memory>

std::shared_ptr<Simulation> initialize_simulation() {
    vec3<size_t> grid_cell_size = {1, 60, 80};
    vec3<double> grid_limits = {
        static_cast<double>(grid_cell_size.z),
        static_cast<double>(grid_cell_size.y),
        static_cast<double>(grid_cell_size.x),
    };
    grid_limits *= PARTICLE_RADIUS;
    Grid grid(std::move(grid_cell_size), std::move(grid_limits));

    std::vector<std::unique_ptr<Particle>> particles;
    for (auto j = 1; j < (grid_cell_size.y / 2); ++j) {
        for (auto i = 1; i < grid_cell_size.x - 1; ++i) {
            particles.emplace_back(std::make_unique<Particle>(make_particle(
                {
                    grid_limits.z / 2,
                    PARTICLE_RADIUS * j,
                    PARTICLE_RADIUS * i,
                },
                1000)));
        }
    }
    auto res = std::make_shared<Simulation>(std::move(grid), std::move(particles));
    return res;
}