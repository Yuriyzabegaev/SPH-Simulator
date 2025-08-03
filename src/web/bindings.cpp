#include "../core/initialize_simulation.hpp"
#include <emscripten/bind.h>

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
    emscripten::register_vector<Particle *>("VectorParticlePtr");
    emscripten::function("make_particle", &make_particle);
}

// Simulation binding
EMSCRIPTEN_BINDINGS(simulation_bindings) {
    emscripten::class_<Simulation>("Simulation")
        .smart_ptr<std::shared_ptr<Simulation>>("shared_ptr<Simulation>")
        .function("get_domain_limits", &Simulation::get_domain_limits)
        .function("apply_central_force", &Simulation::apply_central_force)
        .function("set_gravity", &Simulation::set_gravity)
        .function("set_specific_volume", &Simulation::set_specific_volume)
        .function("add_particle", &Simulation::add_particle)
        .function("update", &Simulation::update)
        .function("set_viscosity", &Simulation::set_viscosity)
        .function("get_particles", &Simulation::get_particles_raw);
    emscripten::function("initialize_simulation", &initialize_simulation);
}
