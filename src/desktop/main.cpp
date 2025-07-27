#include "../core/grid.hpp"
#include "../core/renderer.hpp"
#include "../core/simulation.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

int main() {
    vec3<size_t> grid_cell_size = {1, 60, 80};
    vec3<double> grid_limits = {
        static_cast<double>(grid_cell_size.z),
        static_cast<double>(grid_cell_size.y),
        static_cast<double>(grid_cell_size.x),
    };
    grid_limits *= PARTICLE_RADIUS;
    Grid grid(std::move(grid_cell_size), std::move(grid_limits));

    std::vector<std::unique_ptr<Particle>> particles;
    for (auto i = 12; i < 42; ++i) {
        for (auto j = 12; j < 42; ++j) {
            particles.emplace_back(std::make_unique<Particle>(
                make_particle({PARTICLE_RADIUS / 2, .12 * j, .16 * i}, 1000)));
        }
    }
    Particle p = make_particle({PARTICLE_RADIUS / 2, 3., 0.}, 1000);
    p.velocity.x += 1;
    particles.emplace_back(std::make_unique<Particle>(std::move(p)));

    Simulation simulation(std::move(grid), std::move(particles));
    SFMLRenderer renderer(std::move(simulation));
    renderer.run_until_complete();
}
