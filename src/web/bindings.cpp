#include "../core/grid.hpp"
#include "../core/simulation.hpp"
#include <emscripten/bind.h>
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
    for (auto i = 12; i < 42; ++i) {
        for (auto j = 12; j < 42; ++j) {
            particles.emplace_back(std::make_unique<Particle>(
                make_particle({PARTICLE_RADIUS / 2, .12 * j, .16 * i}, 1000)));
        }
    }
    Particle p = make_particle({PARTICLE_RADIUS / 2, 3., 0.}, 1000);
    p.velocity.x += 1;
    particles.emplace_back(std::make_unique<Particle>(std::move(p)));

    return std::make_shared<Simulation>(std::move(grid), std::move(particles));
}

using vec3d = vec3<double>;
using vec3u = vec3<size_t>;

EMSCRIPTEN_BINDINGS(vec3_bindings) {
    emscripten::value_object<vec3d>("vec3d")
        .field("z", &vec3d::z)
        .field("y", &vec3d::y)
        .field("x", &vec3d::x);

    emscripten::value_object<vec3u>("vec3u")
        .field("z", &vec3u::z)
        .field("y", &vec3u::y)
        .field("x", &vec3u::x);
}

// Particle binding
EMSCRIPTEN_BINDINGS(particle_bindings) {
    emscripten::class_<Particle>("Particle")
        .constructor<vec3d, double>()
        .property("position", &Particle::position)
        .property("velocity", &Particle::velocity)
        .property("mass", &Particle::mass);
}

EMSCRIPTEN_BINDINGS(particle_vector) {
    emscripten::register_vector<Particle *>("VectorParticlePtr");
}

// Grid binding
EMSCRIPTEN_BINDINGS(grid_bindings) {
    emscripten::class_<Grid>("Grid")
        .constructor<vec3u, vec3d>()
        .property("domain_limits_", &Grid::domain_limits_)
        .property("grid_dims_", &Grid::grid_dims_)
        .property("grid_cell_size_", &Grid::grid_cell_size_);
}

// Simulation binding
EMSCRIPTEN_BINDINGS(simulation_bindings) {
    emscripten::class_<Simulation>("Simulation")
        .smart_ptr<std::shared_ptr<Simulation>>("shared_ptr<Simulation>")
        .function("add_particle", &Simulation::add_particle)
        .function("update", &Simulation::update)
        .function("get_particles", &Simulation::get_particles_raw);
    emscripten::function("initialize_simulation", &initialize_simulation);
}
