#include "grid.hpp"
#include "renderer.hpp"
#include "simulation.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

int main() {
    Grid grid(vec3<size_t>{1, 30, 40}, vec3<double>{0.5, 6., 8.});

    std::vector<std::unique_ptr<Particle>> particles;
    for (auto i = 12; i < 42; ++i) {
        for (auto j = 12; j < 42; ++j) {
            particles.emplace_back(std::make_unique<Particle>(
                Particle{{0.25, .12 * j, .16 * i}, 0.1, 4.19}));
        }
    }
    Particle p({0.25, 3., 0.}, 0.1, 4.19);
    p.velocity.x += 1;
    particles.emplace_back(std::make_unique<Particle>(std::move(p)));

    Simulation simulation(std::move(grid), std::move(particles));
    SFMLRenderer renderer(std::move(simulation));
    renderer.run_until_complete();
}
