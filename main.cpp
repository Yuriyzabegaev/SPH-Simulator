#include "grid.hpp"
#include "materials.hpp"
#include "renderer.hpp"
#include "simulation.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

int main() {
    Grid grid(vec3<size_t>{1, 4, 4}, vec3<double>{0.5, 1, 1});

    std::vector<Particle> particles;
    for (auto i = 0; i < 50; ++i) {
        particles.emplace_back(Particle{{0.25, .25, .25}});
    }

    Simulation simulation(std::move(grid), std::move(particles));
    SFMLRenderer renderer(std::move(simulation));
    renderer.run_until_complete();
}
